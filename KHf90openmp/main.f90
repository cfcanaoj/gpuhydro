
module basicmod
  implicit none
  integer::nhy
  integer,parameter::nhymax=20000
  real(8)::time,dt
  data time / 0.0d0 /
  integer,parameter::ngrid=250
  integer,parameter::mgn=2
  integer,parameter::in=ngrid+2*mgn+1 &
 &                  ,jn=ngrid+2*mgn+1 &
 &                  ,kn=1
  integer,parameter::is=mgn+1 &
 &                  ,js=mgn+1 &
 &                  ,ks=1
  integer,parameter::ie=ngrid+mgn &
 &                  ,je=ngrid+mgn &
 &                  ,ke=1

  real(8),parameter:: x1min=-0.5d0,x1max=0.5d0
  real(8),parameter:: x2min=-0.5d0,x2max=0.5d0
  real(8),dimension(in)::x1a,x1b
  real(8),dimension(jn)::x2a,x2b
  real(8),dimension(kn)::x3a,x3b
  
  real(8),dimension(in,jn,kn)::d,et,mv1,mv2,mv3
  real(8),dimension(in,jn,kn)::p,ei,v1,v2,v3,cs
  
  real(8),parameter::gam=5.0d0/3.0d0
  
  integer,parameter::nbc=5
end module basicmod
      
module fluxmod
  use basicmod, only : in,jn,kn
  implicit none
  integer,parameter::nden=1,nve1=2,nve2=3,nve3=4,nene=5,npre=6
  integer,parameter::nhyd=6
  real(8),dimension(nhyd,in,jn,kn):: svc
  
  integer,parameter::mudn=1,muvu=2,muvv=3,muvw=4,muet=5   &
 &                  ,mfdn=6,mfvu=7,mfvv=8,mfvw=9,mfet=10  &
 &                  ,mcsp=11,mvel=12,mpre=13
  integer,parameter:: mflx=5,madd=3

  integer,parameter:: mden=1,mrv1=2,mrv2=3,mrv3=4,meto=5  &
 &                          ,mrvu=muvu,mrvv=muvv,mrvw=muvw
  real(8),dimension(mflx,in,jn,kn):: nflux1,nflux2,nflux3
      
end module fluxmod

      program main
      use basicmod
      use omp_lib
      use mpimod
      implicit none
      real(8)::time_begin,time_end
      logical,parameter::nooutput=.true.
      call InitializeMPI
      if(myid_w == 0) print *, "setup grids and fiels"
      call GenerateGrid
      call GenerateProblem
      call ConsvVariable
      if(myid_w == 0) write(6,*) "entering main loop"
! main loop
      time_begin = omp_get_wtime()
      do nhy=1,nhymax
         if(mod(nhy,100) .eq. 0 .and. .not. nooutput)write(6,*)nhy,time,dt
         call TimestepControl
         call BoundaryCondition
         call StateVevtor
         call NumericalFlux1
         call NumericalFlux2
         call UpdateConsv
         call PrimVariable
         time=time+dt
         if(.not. nooutput)call Output
      enddo
      time_end = omp_get_wtime()

      if(myid_w == 0) print *, "sim time [s]:", time_end-time_begin
      if(myid_w == 0) print *, "time/count/cell", (time_end-time_begin)/(ngrid**2)/nhymax

      call FinalizeMPI
      print *, "program has been finished"
      end program main

      subroutine GenerateGrid
      use basicmod
      implicit none
      real(8)::dx,dy
      integer::i,j,k
      dx=(x1max-x1min)/ngrid
      do i=1,in
         x1a(i) = dx*(i-(mgn+1))+x1min
      enddo
      do i=1,in-1
         x1b(i) = 0.5d0*(x1a(i+1)+x1a(i))
      enddo

      dy=(x2max-x2min)/ngrid
      do j=1,jn
         x2a(j) = dy*(j-(mgn+1))+x2min
      enddo
      do j=1,jn-1
         x2b(j) = 0.5d0*(x2a(j+1)+x2a(j))
      enddo
      
      return
      end subroutine GenerateGrid

      subroutine GenerateProblem
      use basicmod
      implicit none
      integer::i,j,k
      real(8) :: rho1,rho2,Lsm,u1,u2
      data rho1  /  1.0d0 /
      data rho2  /  2.0d0 /
      data u1    /  0.5d0 /
      data u2    / -0.5d0 /
      data Lsm   /  0.025d0 /

      real(8)::pi
      pi=acos(-1.0d0)

      d(:,:,:) = 1.0d0

      do k=ks,ke
      do j=js,je
      do i=is,ie

         if      ( x2b(j) .gt. 0.25d0 )then
            v1(i,j,k) =    u1 - (  u1-  u2)/2.0d0*exp(-( x2b(j)-0.25d0)/Lsm)
             d(i,j,k) =  rho1 - (rho1-rho2)/2.0d0*exp(-( x2b(j)-0.25d0)/Lsm)
         else if (x2b(j) .gt.  0.0d0 )then
            v1(i,j,k) =    u2 + (  u1-  u2)/2.0d0*exp(-( 0.25d0-x2b(j))/Lsm)
             d(i,j,k) =  rho2 + (rho1-rho2)/2.0d0*exp(-( 0.25d0-x2b(j))/Lsm)
         else if (x2b(j) .gt. -0.25d0)then
            v1(i,j,k) =    u2 + (  u1-  u2)/2.0d0*exp(-( x2b(j)+0.25d0)/Lsm)
             d(i,j,k) =  rho2 + (rho1-rho2)/2.0d0*exp(-( x2b(j)+0.25d0)/Lsm)
          else
            v1(i,j,k) =    u1 - (  u1-  u2)/2.0d0*exp(-(-0.25d0-x2b(j))/Lsm)
             d(i,j,k) =  rho1 - (rho1-rho2)/2.0d0*exp(-(-0.25d0-x2b(j))/Lsm)
         endif

          p(i,j,k) = 2.5d0
         v2(i,j,k) = 0.01d0*sin(4.0d0*pi*x1b(i))
         v3(i,j,k) = 0.0d0
      enddo
      enddo
      enddo


      do k=ks,ke
      do j=js,je
      do i=is,ie
          ei(i,j,k) = p(i,j,k)/(gam-1.0d0)
          cs(i,j,k) = sqrt(gam*p(i,j,k)/d(i,j,k))
      enddo
      enddo
      enddo
      
      call BoundaryCondition

      return
      end subroutine GenerateProblem

subroutine BoundaryCondition
  use basicmod
  implicit none
  real(8),dimension(mgn,jn,kn,nbc):: varsendXstt,varsendXend
  real(8),dimension(in,mgn,kn,nbc):: varsendYstt,varsendYend
  real(8),dimension(mgn,jn,kn,nbc):: varrecvXstt,varrecvXend
  real(8),dimension(in,mgn,kn,nbc):: varrecvYstt,varrecvYend
  integer::i,j,k

  k=ks
  do j=1,jn-1
  do i=1,mgn
     varsendXend(i,j,k,1) =  d(ie-mgn+i,j,k)
     varsendXend(i,j,k,2) = ei(ie-mgn+i,j,k)
     varsendXend(i,j,k,3) = v1(ie-mgn+i,j,k)
     varsendXend(i,j,k,4) = v2(ie-mgn+i,j,k)
     varsendXend(i,j,k,5) = v3(ie-mgn+i,j,k)

     varsendXstt(i,j,k,1) =  d(  is+i-1,j,k)
     varsendXstt(i,j,k,2) = ei(  is+i-1,j,k)
     varsendXstt(i,j,k,3) = v1(  is+i-1,j,k)
     varsendXstt(i,j,k,4) = v2(  is+i-1,j,k)
     varsendXstt(i,j,k,5) = v3(  is+i-1,j,k)
  enddo
  enddo

  k=ks
  do i=1,in-1
  do j=1,mgn
     varsendYend(i,j,k,1) =  d(i,je-mgn+j,k)
     varsendYend(i,j,k,2) = ei(i,je-mgn+j,k)
     varsendYend(i,j,k,3) = v1(i,je-mgn+j,k)
     varsendYend(i,j,k,4) = v2(i,je-mgn+j,k)
     varsendYend(i,j,k,5) = v3(i,je-mgn+j,k)

     varsendYstt(i,j,k,1) =  d(i,  js+j-1,k)
     varsendYstt(i,j,k,2) = ei(i,  js+j-1,k)
     varsendYstt(i,j,k,3) = v1(i,  js+j-1,k)
     varsendYstt(i,j,k,4) = v2(i,  js+j-1,k)
     varsendYstt(i,j,k,5) = v3(i,  js+j-1,k)
  enddo
  enddo


  call XbcSendRecv(varsendXstt,varsendXend,varrecvXstt,varrecvXend)
  k=ks
  do j=1,jn-1
  do i=1,mgn
      d(i,j,k) = varrecvXend(i,j,k,1)
     ei(i,j,k) = varrecvXend(i,j,k,2)
     v1(i,j,k) = varrecvXend(i,j,k,3)
     v2(i,j,k) = varrecvXend(i,j,k,4)
     v3(i,j,k) = varrecvXend(i,j,k,5)
  enddo
  enddo

  k=ks
  do j=1,jn-1
  do i=1,mgn
      d(ie+i,j,k) = varrecvXstt(i,j,k,1)
     ei(ie+i,j,k) = varrecvXstt(i,j,k,2)
     v1(ie+i,j,k) = varrecvXstt(i,j,k,3)
     v2(ie+i,j,k) = varrecvXstt(i,j,k,4)
     v3(ie+i,j,k) = varrecvXstt(i,j,k,5)
  enddo
  enddo

  call YbcSendRecv(varsendYstt,varsendYend,varrecvYstt,varrecvYend)

  k=ks
  do i=1,in-1
  do j=1,mgn
      d(i,j,k) = varrecvYstt(i,j,k,1)
     ei(i,j,k) = varrecvYstt(i,j,k,2)
     v1(i,j,k) = varrecvYstt(i,j,k,3)
     v2(i,j,k) = varrecvYstt(i,j,k,4)
     v3(i,j,k) = varrecvYstt(i,j,k,5)
  enddo
  enddo

  k=ks
  do i=1,in-1
  do j=1,mgn
      d(i,je+j,k) = varrecvYend(i,j,k,1)
     ei(i,je+j,k) = varrecvYend(i,j,k,2)
     v1(i,je+j,k) = varrecvYend(i,j,k,3)
     v2(i,je+j,k) = varrecvYend(i,j,k,4)
     v3(i,je+j,k) = varrecvYend(i,j,k,5)
  enddo
  enddo

  return
end subroutine BoundaryCondition

subroutine XbcSendRecv(varsendXstt,varsendXend,varrecvXstt,varrecvXend)
  use   mpimod
  use basicmod
  use mpi
  implicit none
  real(8),dimension(mgn,jn,kn,nbc),intent(in) ::varsendXstt,varsendXend
  real(8),dimension(mgn,jn,kn,nbc),intent(out)::varrecvXstt,varrecvXend
  
  if(ntiles(1) == 1) then
     varrecvXstt(:,:,:,:) = varsendXend(:,:,:,:)
     varrecvXend(:,:,:,:) = varsendXstt(:,:,:,:)
  else

     nreq = nreq + 1         
     call MPI_IRECV(varrecvXstt,mgn*jn*kn*nbc &
    & , MPI_DOUBLE &
    & , n1m,1100, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_ISEND(varsendXstt,mgn*jn*kn*nbc &
    & , MPI_DOUBLE &
    & , n1m, 1200, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_IRECV(varrecvXend,mgn*jn*kn*nbc &
    & , MPI_DOUBLE &
    & , n1p,1200, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_ISEND(varsendXend,mgn*jn*kn*nbc &
    & , MPI_DOUBLE &
    & , n1p, 1100, comm3d, req(nreq), ierr)

     if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
     nreq = 0

  endif

  return
end subroutine XbcSendRecv

subroutine YbcSendRecv(varsendYstt,varsendYend,varrecvYstt,varrecvYend)
  use   mpimod
  use basicmod
  use mpi
  implicit none
  real(8),dimension(in,mgn,kn,nbc),intent(in) ::varsendYstt,varsendYend
  real(8),dimension(in,mgn,kn,nbc),intent(out)::varrecvYstt,varrecvYend

  if(ntiles(2) == 1) then
     varrecvYstt(:,:,:,:) = varsendYend(:,:,:,:)
     varrecvYend(:,:,:,:) = varsendYstt(:,:,:,:)
  else

     nreq = nreq + 1         
     call MPI_IRECV(varrecvYstt,mgn*in*kn*nbc &
    & , MPI_DOUBLE &
    & , n2m, 2100, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_ISEND(varsendYstt,mgn*in*kn*nbc &
    & , MPI_DOUBLE &
    & , n2m, 2200, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_IRECV(varrecvYend,mgn*in*kn*nbc &
    & , MPI_DOUBLE &
    & , n2p,2200, comm3d, req(nreq), ierr)

     nreq = nreq + 1
     call MPI_ISEND(varsendYend,mgn*in*kn*nbc &
    & , MPI_DOUBLE &
    & , n2p, 2100, comm3d, req(nreq), ierr)

     if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
     nreq = 0
  endif

  return
end subroutine YbcSendRecv


      subroutine ConsvVariable
      use basicmod
      implicit none
      integer::i,j,k
!$omp parallel do collapse(3)
      do k=ks,ke
      do j=js,je
      do i=is,ie
          et(i,j,k) = 0.5d0*d(i,j,k)*(     &
     &                    +v1(i,j,k)**2    &
     &                    +v2(i,j,k)**2    &
     &                    +v3(i,j,k)**2)   &
     &                    +ei(i,j,k)   
          mv1(i,j,k) =d(i,j,k)*v1(i,j,k)
          mv2(i,j,k) =d(i,j,k)*v2(i,j,k)
          mv3(i,j,k) =d(i,j,k)*v3(i,j,k)
      enddo
      enddo
      enddo
!$omp end parallel do   
      
      return
      end subroutine Consvvariable

      subroutine PrimVariable
      use basicmod
      implicit none
      integer::i,j,k
      
!$omp parallel do collapse(3)
      do k=ks,ke
      do j=js,je
      do i=is,ie
          v1(i,j,k) = mv1(i,j,k)/d(i,j,k)
          v2(i,j,k) = mv2(i,j,k)/d(i,j,k)
          v3(i,j,k) = mv3(i,j,k)/d(i,j,k)

          ei(i,j,k) =  et(i,j,k)           &
     &          -0.5d0*d(i,j,k)*(          &
     &                    +v1(i,j,k)**2    &
     &                    +v2(i,j,k)**2    &
     &                    +v3(i,j,k)**2)

           p(i,j,k) =  ei(i,j,k)*(gam-1.0d0)
          cs(i,j,k) =  sqrt(gam*p(i,j,k)/d(i,j,k))
      enddo
      enddo
      enddo
!$end omp parallel do
      return
      end subroutine PrimVariable

      subroutine TimestepControl
      use basicmod
      implicit none
      real(8)::dtl1
      real(8)::dtl2
      real(8)::dtl3
      real(8)::dtlocal
      real(8)::dtmin
      integer::i,j,k
      dtmin=1.0d90

!$omp parallel do collapse(3) reduction(min:dtmin) private(dtl1,dtl2,dtl3,dtlocal)      
      do k=ks,ke
      do j=js,je
      do i=is,ie
         dtl1 =(x1a(i+1)-x1a(i))/(abs(v1(i,j,k)) +cs(i,j,k))
         dtl2 =(x2a(j+1)-x2a(j))/(abs(v2(i,j,k)) +cs(i,j,k))
!         dtl3 =(x1a(i+1)-x1a(i))/(abs(v1(i,j,k)) +cs(i,j,k))
         dtlocal = min (dtl1,dtl2)
         if(dtlocal .lt. dtmin) dtmin = dtlocal
      enddo
      enddo
      enddo
!$omp end parallel do
   
      dt = 0.05d0 * dtmin

      return
      end subroutine TimestepControl

      subroutine StateVevtor
      use basicmod
      use fluxmod
      implicit none
      integer::i,j,k
!$omp parallel do collapse(3)
      do k=ks,ke
      do j=1,jn-1
      do i=1,in-1
         svc(nden,i,j,k) =  d(i,j,k)
         svc(nve1,i,j,k) = v1(i,j,k)
         svc(nve2,i,j,k) = v2(i,j,k)
         svc(nve3,i,j,k) = v3(i,j,k)
         svc(nene,i,j,k) = ei(i,j,k)/d(i,j,k)
         svc(npre,i,j,k) = ei(i,j,k)*(gam-1.0d0)
      enddo
      enddo
      enddo
!$omp end parallel do
      
      return
      end subroutine StateVevtor

      subroutine minmod(a,b,d)
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::a,b
      real(8),dimension(nhyd),intent(out)::d
      integer:: n

      do n=1,nhyd
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n)) &
     &                                        ,sign(1.0d0,a(n))*b(n)))
      enddo

      return
      end subroutine minmod


      subroutine vanLeer(dvp,dvm,dv)
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::dvp,dvm
      real(8),dimension(nhyd),intent(out)::dv
      integer:: n

      do n=1,nhyd
         if(dvp(n)*dvm(n) .gt. 0.0d0)then
            dv(n) =2.0d0*dvp(n)*dvm(n)/(dvp(n)+dvm(n))
         else
            dv(n) = 0.0d0
         endif

      enddo

      return
      end subroutine vanLeer

      subroutine MClimiter(a,b,c,d)
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::a,b,c
      real(8),dimension(nhyd),intent(out)::d
      integer:: n

      do n=1,nhyd
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))         &
     &                                  ,sign(1.0d0,a(n))*b(n)   &
     &                                  ,sign(1.0d0,a(n))*c(n)))
      enddo

      return
      end subroutine MClimiter

      subroutine NumericalFlux1
      use basicmod, only: is,ie,in,js,je,jn,ks,ke,kn,gam
      use fluxmod
      implicit none
      integer::i,j,k,n
      real(8),dimension(nhyd):: dsvp,dsvm,dsv
      real(8),dimension(nhyd,in,jn,kn):: leftpr,rigtpr
      real(8),dimension(2*mflx+madd,in,jn,kn):: leftco,rigtco
      real(8),dimension(2*mflx+madd):: leftst,rigtst
      real(8),dimension(mflx):: nflux

!$omp parallel do collapse(3) private(dsv,dsvp,dsvm) 
      do k=ks,ke
      do j=js,je
      do i=is-1,ie+1
        do n=1,nhyd
         dsvp(n) = (svc(n,i+1,j,k) -svc(n,i,j,k)                 )
         dsvm(n) = (                svc(n,i,j,k) - svc(n,i-1,j,k))
        enddo

         call vanLeer(dsvp,dsvm,dsv)
!         call minmod(dsvp,dsvm,dsv)
         do n=1,nhyd
            leftpr(n,i+1,j,k) = svc(n,i,j,k) + 0.5d0*dsv(n)
            rigtpr(n,i  ,j,k) = svc(n,i,j,k) - 0.5d0*dsv(n)
         enddo
      enddo
      enddo
      enddo
!$omp end parallel do

!$omp parallel do collapse(3)
      do k=ks,ke
      do j=js,je
      do i=is,ie+1
         leftco(mudn,i,j,k)=leftpr(nden,i,j,k) ! rho
         leftco(muvu,i,j,k)=leftpr(nve1,i,j,k)*leftpr(nden,i,j,k)     ! rho v_x
         leftco(muvv,i,j,k)=leftpr(nve2,i,j,k)*leftpr(nden,i,j,k)     ! rho v_y
         leftco(muvw,i,j,k)=leftpr(nve3,i,j,k)*leftpr(nden,i,j,k)     ! rho v_z
         leftco(muet,i,j,k)=leftpr(nene,i,j,k)*leftpr(nden,i,j,k)   & ! e_i+ rho v^2/2
     &               +0.5d0*leftpr(nden,i,j,k)*(           &
     &                     +leftpr(nve1,i,j,k)**2         &
     &                     +leftpr(nve2,i,j,k)**2          &
     &                     +leftpr(nve3,i,j,k)**2)

         leftco(mfdn,i,j,k)=leftpr(nden,i,j,k)                   *leftpr(nve1,i,j,k)
         leftco(mfvu,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve1,i,j,k)*leftpr(nve1,i,j,k)  &
     &                     +leftpr(npre,i,j,k)
         leftco(mfvv,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve2,i,j,k)*leftpr(nve1,i,j,k)
         leftco(mfvw,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve3,i,j,k)*leftpr(nve1,i,j,k)
         leftco(mfet,i,j,k)=(leftpr(nene,i,j,k)*leftpr(nden,i,j,k) &
     &               +0.5d0*leftpr(nden,i,j,k)*(   &
     &                     +leftpr(nve1,i,j,k)**2  &
     &                     +leftpr(nve2,i,j,k)**2 &
     &                     +leftpr(nve3,i,j,k)**2) &
     &                     +leftpr(npre,i,j,k) &
     &                       )                                  *leftpr(nve1,i,j,k)

         leftco(mcsp,i,j,k)= sqrt(gam*(gam-1.0d0)*leftpr(nene,i,j,k))
         leftco(mvel,i,j,k)= leftpr(nve1,i,j,k)
         leftco(mpre,i,j,k)= leftpr(npre,i,j,k)


         rigtco(mudn,i,j,k)=rigtpr(nden,i,j,k)
         rigtco(muvu,i,j,k)=rigtpr(nve1,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvv,i,j,k)=rigtpr(nve2,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvw,i,j,k)=rigtpr(nve3,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muet,i,j,k)=rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k)   &
     &               +0.5d0*rigtpr(nden,i,j,k)*(   &
     &                     +rigtpr(nve1,i,j,k)**2 &
     &                     +rigtpr(nve2,i,j,k)**2 &
     &                     +rigtpr(nve3,i,j,k)**2)

         rigtco(mfdn,i,j,k)=rigtpr(nden,i,j,k)                   *rigtpr(nve1,i,j,k)
         rigtco(mfvu,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve1,i,j,k)*rigtpr(nve1,i,j,k)  &
     &                     +rigtpr(npre,i,j,k)
         rigtco(mfvv,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve2,i,j,k)*rigtpr(nve1,i,j,k)
         rigtco(mfvw,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve3,i,j,k)*rigtpr(nve1,i,j,k)
         rigtco(mfet,i,j,k)=(rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k)  &
     &               +0.5d0*rigtpr(nden,i,j,k)*(  &
     &                     +rigtpr(nve1,i,j,k)**2 &
     &                     +rigtpr(nve2,i,j,k)**2 &
     &                     +rigtpr(nve3,i,j,k)**2) &
     &                     +rigtpr(npre,i,j,k) &
     &                      )                                    *rigtpr(nve1,i,j,k)

         rigtco(mcsp,i,j,k)= sqrt(gam*(gam-1.0d0)*rigtpr(nene,i,j,k))
         rigtco(mvel,i,j,k)= rigtpr(nve1,i,j,k)
         rigtco(mpre,i,j,k)= rigtpr(npre,i,j,k)

      enddo
      enddo
      enddo 
!$omp end parallel do

!$omp parallel do collapse(3) private(leftst,rigtst,nflux) 
      do k=ks,ke
      do j=js,je
      do i=is,ie+1
         leftst(:)=leftco(:,i,j,k)
         rigtst(:)=rigtco(:,i,j,k)
         call HLLE(leftst,rigtst,nflux)
         nflux1(mden,i,j,k)=nflux(mden)
         nflux1(mrv1,i,j,k)=nflux(mrvu)
         nflux1(mrv2,i,j,k)=nflux(mrvv)
         nflux1(mrv3,i,j,k)=nflux(mrvw)
         nflux1(meto,i,j,k)=nflux(meto)
      enddo
      enddo

      enddo
!$omp end parallel do
      
      return
      end subroutine Numericalflux1

      subroutine NumericalFlux2
      use basicmod, only: is,ie,in,js,je,jn,ks,ke,kn,gam
      use fluxmod
      implicit none
      integer::i,j,k
      real(8),dimension(nhyd):: dsvp,dsvm,dsvc,dsv
      real(8),dimension(nhyd,in,jn,kn):: leftpr,rigtpr
      real(8),dimension(2*mflx+madd,in,jn,kn):: leftco,rigtco
      real(8),dimension(2*mflx+madd):: leftst,rigtst
      real(8),dimension(mflx):: nflux

!$omp parallel do collapse(3) private(dsv,dsvp,dsvm)
      do k=ks,ke
      do i=is,ie
      do j=js-1,je+1
         dsvp(:) = (svc(:,i,j+1,k) -svc(:,i,j,k)                 )
         dsvm(:) = (                svc(:,i,j,k) - svc(:,i,j-1,k))

         call vanLeer(dsvp,dsvm,dsv)
!         call minmod(dsvp,dsvm,dsv)
         leftpr(:,i,j+1,k) = svc(:,i,j,k) + 0.5d0*dsv(:)
         rigtpr(:,i,j  ,k) = svc(:,i,j,k) - 0.5d0*dsv(:)

!         leftpr(:,i,j,k) = svc(:,i,j-1,k)
!         rigtpr(:,i,j,k) = svc(:,i,j  ,k)

     enddo
     enddo
     enddo
!$omp end parallel do

!$omp parallel do collapse(3)
      do k=ks,ke
      do i=is,ie
      do j=js,je+1
         leftco(mudn,i,j,k)=leftpr(nden,i,j,k)
         leftco(muvw,i,j,k)=leftpr(nve1,i,j,k)*leftpr(nden,i,j,k)
         leftco(muvu,i,j,k)=leftpr(nve2,i,j,k)*leftpr(nden,i,j,k) ! rho v
         leftco(muvv,i,j,k)=leftpr(nve3,i,j,k)*leftpr(nden,i,j,k)
         leftco(muet,i,j,k)=leftpr(nene,i,j,k)*leftpr(nden,i,j,k)  &
     &               +0.5d0*leftpr(nden,i,j,k)*(                   &
     &                     +leftpr(nve1,i,j,k)**2                  &
     &                     +leftpr(nve2,i,j,k)**2                  &
     &                     +leftpr(nve3,i,j,k)**2)

         leftco(mfdn,i,j,k)=leftpr(nden,i,j,k)                   *leftpr(nve2,i,j,k) ! rho v
         leftco(mfvw,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve1,i,j,k)*leftpr(nve2,i,j,k)
         leftco(mfvu,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve2,i,j,k)*leftpr(nve2,i,j,k) &
     &                     +leftpr(npre,i,j,k)
         leftco(mfvv,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve3,i,j,k)*leftpr(nve2,i,j,k)
         leftco(mfet,i,j,k)=(leftpr(nene,i,j,k)*leftpr(nden,i,j,k) &
     &               +0.5d0*leftpr(nden,i,j,k)*(   &
     &                     +leftpr(nve1,i,j,k)**2  &
     &                     +leftpr(nve2,i,j,k)**2  &
     &                     +leftpr(nve3,i,j,k)**2) &
     &                     +leftpr(npre,i,j,k)     &
     &                                       )*leftpr(nve2,i,j,k)

         leftco(mcsp,i,j,k)= sqrt(gam*(gam-1.0d0)*leftpr(nene,i,j,k))
         leftco(mvel,i,j,k)= leftpr(nve2,i,j,k)
         leftco(mpre,i,j,k)= leftpr(npre,i,j,k)


         rigtco(mudn,i,j,k)=rigtpr(nden,i,j,k)
         rigtco(muvw,i,j,k)=rigtpr(nve1,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvu,i,j,k)=rigtpr(nve2,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvv,i,j,k)=rigtpr(nve3,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muet,i,j,k)=rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) &
     &               +0.5d0*rigtpr(nden,i,j,k)*(   &
     &                     +rigtpr(nve1,i,j,k)**2  &
     &                     +rigtpr(nve2,i,j,k)**2  &
     &                     +rigtpr(nve3,i,j,k)**2) 

         rigtco(mfdn,i,j,k)=rigtpr(nden,i,j,k)                   *rigtpr(nve2,i,j,k)
         rigtco(mfvw,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve1,i,j,k)*rigtpr(nve2,i,j,k)
         rigtco(mfvu,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve2,i,j,k)*rigtpr(nve2,i,j,k) &
     &                     +rigtpr(npre,i,j,k)
         rigtco(mfvv,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve3,i,j,k)*rigtpr(nve2,i,j,k)
         rigtco(mfet,i,j,k)=(rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) &
     &               +0.5d0*rigtpr(nden,i,j,k)*(    &
     &                     +rigtpr(nve1,i,j,k)**2   &
     &                     +rigtpr(nve2,i,j,k)**2   &
     &                     +rigtpr(nve3,i,j,k)**2)  &
     &                     +rigtpr(npre,i,j,k)      &
     &                                       )*rigtpr(nve2,i,j,k)

         rigtco(mcsp,i,j,k)= sqrt(gam*(gam-1.0d0)*rigtpr(nene,i,j,k))
         rigtco(mvel,i,j,k)= rigtpr(nve2,i,j,k)
         rigtco(mpre,i,j,k)= rigtpr(npre,i,j,k)

      enddo
      enddo
      enddo
!$omp end parallel do

!$omp parallel do collapse(3) private(leftst,rigtst,nflux)
      do k=ks,ke
      do i=is,ie
      do j=js,je+1
         leftst(:)=leftco(:,i,j,k)
         rigtst(:)=rigtco(:,i,j,k)
         call HLLE(leftst,rigtst,nflux)
         nflux2(mden,i,j,k)=nflux(mden)
         nflux2(mrv1,i,j,k)=nflux(mrvw)
         nflux2(mrv2,i,j,k)=nflux(mrvu) ! mrv2=3, mrvu=2
         nflux2(mrv3,i,j,k)=nflux(mrvv)
         nflux2(meto,i,j,k)=nflux(meto)
      enddo
      enddo
      enddo
!$omp end parallel do

      return
      end subroutine Numericalflux2

      subroutine HLLE(leftst,rigtst,nflux)
      use fluxmod
      implicit none
      real(8),dimension(2*mflx+madd),intent(in)::leftst,rigtst
      real(8),dimension(mflx),intent(out)::nflux
      real(8),dimension(mflx)::ul,ur,fl,fr
      real(8)::csl,csr
      real(8):: vl, vr
      real(8):: sl, sr

      ul(1:mflx) = leftst(1:mflx)
      fl(1:mflx) = leftst(mflx+1:2*mflx)
      ur(1:mflx) = rigtst(1:mflx)
      fr(1:mflx) = rigtst(mflx+1:2*mflx)
      csl=leftst(mcsp)
      csr=rigtst(mcsp)
       vl=leftst(mvel)
       vr=rigtst(mvel)

       sl = min(vl,vr) - max(csl,csr)
       sl = min(0.0d0,sl)
       sr = max(vl,vr) + max(csl,csr)
       sr = max(0.0d0,sr)

       nflux(:) = (sr*fl(:)-sl*fr(:) +sl*sr*(ur(:)-ul(:)))/(sr-sl)

      return
      end subroutine HLLE

      subroutine UpdateConsv
      use basicmod
      use fluxmod
      implicit none
      integer::i,j,k

!$omp parallel do collapse(3)
      do k=ks,ke
      do j=js,je
      do i=is,ie
         
         d(i,j,k) = d(i,j,k)  &
     & +dt*( &
     &  (- nflux1(mden,i+1,j,k) &
     &   + nflux1(mden,i  ,j,k))/(x1a(i+1)-x1a(i))  &
     & +(- nflux2(mden,i,j+1,k) &
     &   + nflux2(mden,i,j  ,k))/(x2a(j+1)-x2a(j))  &
     &      )

         mv1(i,j,k) = mv1(i,j,k)  &
     & +dt*( &
     &  (- nflux1(mrv1,i+1,j,k) &
     &   + nflux1(mrv1,i  ,j,k))/(x1a(i+1)-x1a(i))  &
     & +(- nflux2(mrv1,i,j+1,k) &
     &   + nflux2(mrv1,i,j  ,k))/(x2a(j+1)-x2a(j)) & 
     &      )

         mv2(i,j,k) = mv2(i,j,k)  &
     & +dt*( &
     &  (- nflux1(mrv2,i+1,j,k) &
     &   + nflux1(mrv2,i  ,j,k))/(x1a(i+1)-x1a(i)) & 
     & +(- nflux2(mrv2,i,j+1,k) &
     &   + nflux2(mrv2,i,j  ,k))/(x2a(j+1)-x2a(j))  &
     &      )

         mv3(i,j,k) = mv3(i,j,k)  &
     & +dt*( &
     &  (- nflux1(mrv3,i+1,j,k) &
     &   + nflux1(mrv3,i  ,j,k))/(x1a(i+1)-x1a(i))  &
     & +(- nflux2(mrv3,i,j+1,k) &
     &   + nflux2(mrv3,i,j  ,k))/(x2a(j+1)-x2a(j)) & 
     &      )

          et(i,j,k) = et(i,j,k)  &
     & +dt*( &
     &  (- nflux1(meto,i+1,j,k) &
     &   + nflux1(meto,i  ,j,k))/(x1a(i+1)-x1a(i)) & 
     & +(- nflux2(meto,i,j+1,k) &
     &   + nflux2(meto,i,j  ,k))/(x2a(j+1)-x2a(j)) & 
     &      )
      enddo
      enddo
      enddo
!$end omp parallel do

      return
      end subroutine UpdateConsv

      subroutine Output
      use basicmod
      implicit none
      integer::i,j,k
      character(20),parameter::dirname="snapshots/"
      character(40)::filename
      real(8),save::tout
      data tout / 0.0d0 / 
      real(8),parameter:: dtout=1.0d-2
      integer::nout
      data nout / 1 /
      integer,parameter::unitout=13
      logical,save::is_inited
      data is_inited / .false. /
      
      if(.not. is_inited)then
         call makedirs(dirname)
         is_inited = .true.
      endif
      if(time .lt. tout+dtout) return

      write(filename,'(a2,i5.5,a4)')"Sc",nout,".xss"
      filename = trim(dirname)//filename
      open(unitout,file=filename,status='replace',form='formatted') 

      write(unitout,*) "# ",time
      k=ks
      do j=js,je
      do i=is,ie
         write(unitout,'(7(1x,E12.3))') x1b(i),x2b(j),d(i,j,k),v1(i,j,k),v2(i,j,k),v3(i,j,k),p(i,j,k)
      enddo
         write(unitout,*)
      enddo
      close(unitout)

      write(6,*) "output:",nout,time

      nout=nout+1
      tout=time

      return
      end subroutine Output

      subroutine makedirs(outdir)
        implicit none
        character(len=*), intent(in) :: outdir
        character(len=256) command
        write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
        write(*, *) trim(command)
        call system(command)
      end subroutine makedirs
