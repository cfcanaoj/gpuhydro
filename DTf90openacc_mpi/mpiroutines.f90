module mpimod
  use mpi
  implicit none
  integer, parameter :: mreq  = 300
  integer :: stat(MPI_STATUS_SIZE,mreq)                     
  integer :: req(mreq)
  
  integer :: ierr,myid_w, nprocs_w
  integer :: mpi_comm_hyd,myid_hyd, nprocs_hyd
  integer :: comm3d,myid, nprocs
  logical :: periodic(3)
  integer :: ntiles(3), coords(3)
  logical :: reorder
  integer :: n1m, n1p, n2m, n2p, n3m, n3p
  integer :: nreq, nsub
  integer ::   gpuid, ngpus
!$acc declare create(myid_w)
contains
subroutine InitializeMPI
  use openacc
  implicit none
  integer::key,color
  integer::np_hyd

! Initialize MPI
  call MPI_INIT( ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs_w, ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid_w  , ierr )
  
  ntiles(1)=1
  ntiles(2)=2
  ntiles(3)=2
  periodic(1)=.true.
  periodic(2)=.true.
  periodic(3)=.true.
  if(myid_w == 0) then
     print *, "MPI process=",nprocs_w
     print *, "decomposition=",ntiles(1),ntiles(2),ntiles(3)
  endif

  call MPI_BCAST(ntiles,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(periodic,3,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

! Making 3D strucure
  np_hyd = ntiles(1)*ntiles(2)*ntiles(3)
  color = int(myid_w/np_hyd)
  key   = myid_w   
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,mpi_comm_hyd,ierr)
  call MPI_COMM_SIZE( mpi_comm_hyd, nprocs_hyd, ierr )
  call MPI_COMM_RANK( mpi_comm_hyd, myid_hyd , ierr )     
  
! Create a virtual Cartesian topology for the domain decomposition.
!
  call MPI_CART_CREATE( mpi_comm_hyd, 3, ntiles, periodic &
       &                    , reorder, comm3d, ierr )
  call MPI_COMM_RANK( comm3d, myid,     ierr )
  call MPI_COMM_SIZE( comm3d, nprocs,   ierr )
!
! Find the ranks of my neighbors; find my virtual Cartesian coords.
!
  call MPI_CART_SHIFT( comm3d, 0, 1, n1m, n1p, ierr )
  call MPI_CART_SHIFT( comm3d, 1, 1, n2m, n2p, ierr )
  call MPI_CART_SHIFT( comm3d, 2, 1, n3m, n3p, ierr )
  !
  call MPI_CART_COORDS( comm3d, myid, 3, coords, ierr )

  ngpus = acc_get_num_devices(acc_device_nvidia)
  if(myid_w == 0) then
     print *, "num of GPUs = ", ngpus
  end if

  gpuid = mod(myid_w, ngpus)
  if(ngpus == 0) gpuid = -1
  if(gpuid >= 0) then
     call acc_set_device_num(gpuid, acc_device_nvidia)
  end if
  
!$acc update device (myid_w)
  return
end subroutine InitializeMPI

subroutine FinalizeMPI
  implicit none
  call MPI_FINALIZE(ierr)
end subroutine FinalizeMPI

subroutine MPIminfind(bufinp, bufout)
  implicit none
  real(8),intent(in) :: bufinp(2)
  real(8),intent(out):: bufout(2)

!$acc host_data use_device(bufinp,bufout)
       call MPI_ALLREDUCE( bufinp(1), bufout(1), 1 &
     &                   , MPI_2DOUBLE_PRECISION   &
     &                   , MPI_MINLOC, comm3d, ierr)      
!$acc end host_data

end subroutine MPIminfind

subroutine MPImaxfind(bufinp, bufout)
  implicit none
  real(8),intent(in) :: bufinp(2)
  real(8),intent(out):: bufout(2)
!$acc host_data use_device(bufinp,bufout)
       call MPI_ALLREDUCE( bufinp(1), bufout(1), 1 &
     &                   , MPI_2DOUBLE_PRECISION   &
     &                   , MPI_MAXLOC, comm3d, ierr)
!$acc end host_data
end subroutine MPImaxfind

end module mpimod


!==================================================
! DATA IO
!==================================================

module mpiiomod
  implicit none
  private
  integer::SAG1D,SAG2D,SAG3D,SAD3D
  
  integer,dimension(3):: ntotal
  integer,dimension(3):: npart
  integer:: nvars,nvarg
  character(len= 2),parameter :: id ="DT"
  character(len=10),parameter :: datadir="bindata/"
    
  real(8),dimension(:,:),allocatable,save :: gridX, gridY, gridZ
  real(8),dimension(:,:,:,:),allocatable,save :: data3D
  public ntotal,npart
  public nvars,nvarg
  
  public gridX,gridY,gridZ,data3D
  public MPIOutputBindary
contains  
  subroutine MPIOutputBindary(timeid)
    use mpimod
    implicit NONE
    integer,intent(in) :: timeid
    integer :: i, j, k, l, m, n

    integer::iss,jss,kss
    integer::iee,jee,kee
    integer::itot,jtot,ktot
      
    integer strtoi 
    character(len=15) :: unffile
    character(len=40)::usrfile
    character(len=30) :: fpathbin,fpathunf
    integer, parameter :: unitunf=560
    integer,save:: unitd3d,unitg1d, unitg2d, unitg3d, unitg0d
    data unitd3d / 512 /
    data unitg1d / 513 /
    data unitg2d / 514 /
    data unitg3d / 515 /
    data unitg0d / 516 /

    logical :: fileflag

    logical,save :: is_inited
    data is_inited / .false. /
    
    integer,dimension(4)::Asize,Ssize,Start
    integer(kind=MPI_OFFSET_KIND) idisp
    data idisp / 0 /


!   print *, "p1"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D GRID PREPARE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    init1D: if(.not. is_inited )then

       Asize(1) = nvarg
       Ssize(1) = nvarg
       Start(1) = 0
       
       Asize(2) = ntotal(1)+1 ! total izones + edge
       Ssize(2) = npart(1) ! izones in 1 process
       if(coords(1) .eq. ntiles(1)-1)Ssize(2)=Ssize(2)+1  ! + edge
       Start(2) = npart(1) * coords(1)
              
       call MPI_TYPE_CREATE_SUBARRAY( &
     & 2, & ! dimension of array
     & Asize,Ssize,Start, &
     & MPI_ORDER_FORTRAN, &
     & MPI_DOUBLE_PRECISION,& 
     & SAG1D, &! Data type of Subarray for Grid 1D
     & ierr)
       
       call MPI_TYPE_COMMIT(SAG1D,ierr)

      write(usrfile,"(a3,a2)")'g1d',id
      fpathbin = trim(datadir)//usrfile
      call MPI_FILE_OPEN(MPI_COMM_WORLD, &
     &                         fpathbin, &  ! file path
     &  MPI_MODE_WRONLY+MPI_MODE_CREATE, &
     &            MPI_INFO_NULL,unitg1d,ierr)

      call MPI_FILE_SET_VIEW( &
     &   unitg1d, &! file path
     &     idisp, &! 
     & MPI_DOUBLE_PRECISION, & 
     &     SAG1D, &! data type
     & 'NATIVE', MPI_INFO_NULL,ierr)

      call MPI_FILE_WRITE_ALL( &
     &   unitg1d, &  ! file path
     &     gridX, &  ! the data
     & Ssize(2)*nvarg, &! total data number
     & MPI_DOUBLE_PRECISION, & 
     & mpi_status_ignore, &
     & ierr)
      call MPI_FILE_CLOSE(unitg1d,ierr)
      
   endif init1D
!   print *, "p2"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D GRID PREPARE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   init2D: if(.not. is_inited )then
      Asize(1) = nvarg
      Ssize(1) = nvarg
      Start(1) = 0
      
      Asize(2) = ntotal(2)+1 ! total jzones + edge
      Ssize(2) = npart(2)    ! jzones in 1 process
      if(coords(2) .eq. ntiles(2)-1)Ssize(2)=Ssize(2)+1  ! + edge
      Start(2) = npart(2) * coords(2)
      call MPI_TYPE_CREATE_SUBARRAY(&
     & 2, & ! dimension of array
     & Asize,Ssize,Start,&
     & MPI_ORDER_FORTRAN,&
     & MPI_DOUBLE_PRECISION,&
     & SAG2D, &! Data type of Subarray for Grid 2D
     & ierr)
      call MPI_TYPE_COMMIT(SAG2D,ierr)
         
      write(usrfile,"(a3,a2)")'g2d',id
      fpathbin = trim(datadir)//usrfile
 
      call MPI_FILE_OPEN(MPI_COMM_WORLD, &
     &                         fpathbin, & ! file path
     &  MPI_MODE_WRONLY+MPI_MODE_CREATE, &
     &            MPI_INFO_NULL,unitg2d,ierr)
      call MPI_FILE_SET_VIEW(&
     &   unitg2d, & ! file ID
     &     idisp, & ! 
     & MPI_DOUBLE_PRECISION,&
     &     SAG2D, & ! data type
     & 'NATIVE', MPI_INFO_NULL,ierr)

      call MPI_FILE_WRITE_ALL( &
     &   unitg2d,  &! file ID
     &     gridY,  &! the data
     & Ssize(2)*nvarg,& ! total data number
     & MPI_DOUBLE_PRECISION,&
     & mpi_status_ignore,&
     & ierr)
      call MPI_FILE_CLOSE(unitg2d,ierr)

      endif init2D
      
      init3D: if(.not. is_inited )then
         Asize(1) = nvarg
         Ssize(1) = nvarg
         Start(1) = 0
         Asize(2) = ntotal(3)+1  ! total kzones+edge
         Ssize(2) = npart(3) ! kzones in 1 process
         if(coords(3) .eq. ntiles(3)-1)Ssize(2)=Ssize(2)+1  ! + edge
         Start(2) = npart(3) * coords(3)
         call MPI_TYPE_CREATE_SUBARRAY(&
     & 2, & ! dimension of array
     & Asize,Ssize,Start, &
     & MPI_ORDER_FORTRAN, &
     & MPI_DOUBLE_PRECISION,&
     & SAG3D, &! Data type of Subarray for Grid 3D
     & ierr)
         call MPI_TYPE_COMMIT(SAG3D,ierr)

      write(usrfile,"(a3,a2)")'g3d',id
      fpathbin = trim(datadir)//usrfile
      
      call MPI_FILE_OPEN(MPI_COMM_WORLD,&
     &                         fpathbin,&  ! file path
     &  MPI_MODE_WRONLY+MPI_MODE_CREATE,&
     &            MPI_INFO_NULL,unitg3d,ierr)
      call MPI_FILE_SET_VIEW(&
     &  unitg3d, &  ! file path
     &    idisp, &  ! 
     & MPI_DOUBLE_PRECISION, &
     &     SAG3D, & ! data type
     & 'NATIVE', MPI_INFO_NULL,ierr)

      call MPI_FILE_WRITE_ALL(&
     &  unitg3d, & ! file path
     &    gridZ, & ! the data
     & Ssize(2)*nvarg,& ! total data number
     & MPI_DOUBLE_PRECISION, &
     & mpi_status_ignore, &
     & ierr)
      call MPI_FILE_CLOSE(unitg3d,ierr)
      
   endif init3D
          
    initdata: if(.not. is_inited )then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DATA PREPARE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       Asize(1) = nvars
       Ssize(1) = nvars
       Start(1) = 0 
       Asize(2) = ntotal(1) ! total zones for 1D 
       Asize(3) = ntotal(2) ! total zones for 2D
       Asize(4) = ntotal(3) ! total zones for 3D
       Ssize(2) =  npart(1) ! partial zones in 1 process 
       Ssize(3) =  npart(2) ! partial zones in 1 process 
       Ssize(4) =  npart(3) ! partial zones in 1 process 
       Start(2) =  npart(1) * coords(1)
       Start(3) =  npart(2) * coords(2)
       Start(4) =  npart(3) * coords(3)

       call MPI_TYPE_CREATE_SUBARRAY(&
     & 4, &! dimension of array
     & Asize,Ssize,Start,&
     & MPI_ORDER_FORTRAN,&
     & MPI_DOUBLE_PRECISION,&
     & SAD3D,& ! Data type of Subarray for data 3D
     & ierr)
       call MPI_TYPE_COMMIT(SAD3D,ierr)

    endif initdata
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DATA WRITE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(usrfile,"(a3,a2,a1,i5.5)")'d3d',id,'.',timeid
    fpathbin = trim(datadir)//usrfile
    
      call MPI_FILE_OPEN(MPI_COMM_WORLD, &
     &                         fpathbin, & ! file path
     &  MPI_MODE_WRONLY+MPI_MODE_CREATE, &
     &            MPI_INFO_NULL,unitd3d,ierr)
      call MPI_FILE_SET_VIEW(&
     &  unitd3d,  &! file path
     &     idisp, & ! 
     & MPI_DOUBLE_PRECISION,& 
     &     SAD3D, & ! data type
     & 'NATIVE', MPI_INFO_NULL,ierr)

      call MPI_FILE_WRITE_ALL(&
     &   unitd3d,  &! file path
     &    data3D,  &! the data
     & npart(1)*npart(2)*npart(3)*nvars,& ! total data number
     & MPI_DOUBLE_PRECISION,&  
     & mpi_status_ignore,&
     &      ierr)
      call MPI_FILE_CLOSE(unitd3d,ierr)
      
      is_inited = .true.

      return
    end subroutine MPIOutputBindary
  end module mpiiomod
