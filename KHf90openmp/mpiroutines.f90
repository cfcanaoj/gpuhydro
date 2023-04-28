module mpipara
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

contains
subroutine InitializeMPI
  implicit none
  integer::key,color
  integer::np_hyd

! Initialize MPI
  call MPI_INIT( ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs_w, ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid_w  , ierr )
  
  ntiles(1)=2
  ntiles(2)=2
  ntiles(3)=1
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
      
  return
end subroutine InitializeMPI

subroutine FinalizeMPI
  implicit none
  call MPI_FINALIZE(ierr)
end subroutine FinalizeMPI

end module  mpipara
