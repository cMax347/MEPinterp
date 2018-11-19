c**********************************c
c  read in input file and set up   c
c  the Hamiltonian and velocities  c
c                  J.-P. Hanke     c
c**********************************c
      module m_setham_FeMn
      implicit none

      real, public, save :: t,t2
      real, public, save :: lambda,ax,ay,az
      real, public, save :: minef,maxef
      real, public, save :: phi(4),theta(4)
      logical, public, save :: l_path
      integer, public, save :: numsteps,grid(3)

      contains

      subroutine read_inp(irank,l_write)
      implicit none
      integer,intent(in) :: irank
      logical,intent(in) :: l_write
      logical :: l_exist
      integer :: i
      real :: cfac, mag(3,4)

      inquire(file='inp_params_3q',exist=l_exist)
      if(.not.l_exist) stop 'provide inp_params_3q'

      phi=0.0
      open(444,file='inp_params_3q')
      read(444,*)grid(1),grid(2),grid(3)
      read(444,*)t,t2,lambda
      read(444,*)phi(1),theta(1) !A
      read(444,*)phi(2),theta(2) !B
      read(444,*)phi(3),theta(3) !C
      read(444,*)phi(4),theta(4) !D
      read(444,*)minef,maxef,numsteps
      read(444,*)l_path
      if(l_path) then 
       read(444,*)grid(1)
       grid(1) = 4*(grid(1)/4)
       grid(2) = 1
       grid(3) = 1
      endif
      close(444)

      ax = 1.
      ay = 2./sqrt(3.)
      az = sqrt(2./3.)*t/t2

      if(irank.eq.0 .and. l_write) then
       write(*,'(a)')'Parameters:'
       write(*,'(a,3(i5,1x))')'k-grid   =',grid
       write(*,'(a,1x,2(f8.3,1x))')'hoppings =',t,t2
       write(*,'(a,1x,1f8.3)')'exchange =',lambda

       write(*,*)'  S4 ------ S2  '
       write(*,*)'   \        /   '
       write(*,*)'    \  S1  /    '
       write(*,*)'y    \    /     '
       write(*,*)'|     \  /      '
       write(*,*)'|      S3       '
       write(*,*)'O---x           '

       cfac = 3.1415926535897931/180.
       do i=1,4
        mag(1,i) = cos(phi(i)*cfac)*sin(theta(i)*cfac)
        mag(2,i) = sin(phi(i)*cfac)*sin(theta(i)*cfac)
        mag(3,i) =                  cos(theta(i)*cfac)
        write(*,'(a,1x,i1,1x,a,1x,3(f10.5,1x))')'S',i,' = ',mag(1:3,i)
       enddo
      endif

      end subroutine read_inp


      subroutine init_ham(
     >           kpt,num_wann,ham,vx,vy,vz )
      implicit none
      real, intent(in) :: kpt(3)
      integer, intent(in) :: num_wann
      complex, intent(inout) :: ham(num_wann,num_wann)
      complex, optional,intent(inout) :: vx(num_wann,num_wann)
      complex, optional,intent(inout) :: vy(num_wann,num_wann)
      complex, optional,intent(inout) :: vz(num_wann,num_wann)
     
      real :: ht1(3),ht2(3),convfac
      real :: pi,kx,ky,kz,sintheta(4),costheta(4)
      real :: thconv(4),phconv(4),faz
      complex :: phasphi(4)
      integer :: i,j

       pi = 3.1415926535897931
       convfac = pi/180.0
       thconv = theta*convfac
       phconv = phi*convfac
       kx=kpt(1)
       ky=kpt(2)
       kz=kpt(3)

       ! exchange
       do i=1,4
        sintheta(i) = sin(thconv(i))
        costheta(i) = cos(thconv(i))
        phasphi(i) = cmplx(cos(phconv(i)),-sin(phconv(i)))
       enddo

       ! intra-layer hopping
       ht1(1) = -2*t*cos(pi*kx-3./2.*pi*ky)
       ht1(2) = -2*t*cos(2.*pi*kx)
       ht1(3) = -2*t*cos(pi*kx+3./2.*pi*ky)

       ! inter-layer hopping
       ht2(1) =-2*t2*cos(-pi*kx-pi/2.*ky-2.*pi*kz)
       ht2(2) =-2*t2*cos(          pi*ky-2.*pi*kz)
       ht2(3) =-2*t2*cos( pi*kx-pi/2.*ky-2.*pi*kz)

       ! HAMILTONIAN
       ham = cmplx(0.,0.)
       do i=1,4 ! exchange
        ham(i,i) = lambda*costheta(i)
        ham(i+4,i+4) = -lambda*costheta(i)
        ham(i,i+4) = lambda*phasphi(i)*sintheta(i)
       enddo
       do i=0,4,4 ! intra- and inter-layer hopping
        ham(1+i,2+i)=ht1(1)+ht2(1)
        ham(1+i,3+i)=ht1(2)+ht2(2)
        ham(1+i,4+i)=ht1(3)+ht2(3)
        ham(2+i,3+i)=ht1(3)+ht2(3)
        ham(2+i,4+i)=ht1(2)+ht2(2)
        ham(3+i,4+i)=ht1(1)+ht2(1)
       enddo

       ! VELOCITY VX
       if(present(vx)) then
       vx = cmplx(0.,0.)

       ht1(1) = ax*t*sin(pi*kx-3./2.*pi*ky)
       ht1(2) = 2.*ax*t*sin(2.*pi*kx)
       ht1(3) = ax*t*sin(pi*kx+3./2.*pi*ky)

       ht2(1) = -ax*t2*sin(-pi*kx-pi/2.*ky-2.*pi*kz)
       ht2(2) =  0.0
       ht2(3) =  ax*t2*sin( pi*kx-pi/2.*ky-2.*pi*kz)

       do i=0,4,4 ! intra- and inter-layer hopping
        vx(1+i,2+i)=ht1(1)+ht2(1)
        vx(1+i,3+i)=ht1(2)+ht2(2)
        vx(1+i,4+i)=ht1(3)+ht2(3)
        vx(2+i,3+i)=ht1(3)+ht2(3)
        vx(2+i,4+i)=ht1(2)+ht2(2)
        vx(3+i,4+i)=ht1(1)+ht2(1)
       enddo
       endif

       ! VELOCITY VY
       if(present(vy))then
       vy = cmplx(0.,0.)

       ht1(1) = -3./2.*ay*t*sin(pi*kx-3./2.*pi*ky)
       ht1(2) =  0.0
       ht1(3) =  3./2.*ay*t*sin(pi*kx+3./2.*pi*ky)

       ht2(1) = -ay/2.*t2*sin(-pi*kx-pi/2.*ky-2.*pi*kz)
       ht2(2) =  ay*t2*sin(             pi*ky-2.*pi*kz)
       ht2(3) = -ay/2.*t2*sin( pi*kx-pi/2.*ky-2.*pi*kz)

       do i=0,4,4 ! intra- and inter-layer hopping
        vy(1+i,2+i)=ht1(1)+ht2(1)
        vy(1+i,3+i)=ht1(2)+ht2(2)
        vy(1+i,4+i)=ht1(3)+ht2(3)
        vy(2+i,3+i)=ht1(3)+ht2(3)
        vy(2+i,4+i)=ht1(2)+ht2(2)
        vy(3+i,4+i)=ht1(1)+ht2(1)
       enddo
       endif

       ! VELOCITY VZ
       if(present(vz))then
       vz = cmplx(0.,0.)

       ht1(1) = 0.0
       ht1(2) = 0.0
       ht1(3) = 0.0

       ht2(1) =-az*2.*t2*sin(-pi*kx-pi/2.*ky-2.*pi*kz)
       ht2(2) =-az*2.*t2*sin(          pi*ky-2.*pi*kz)
       ht2(3) =-az*2.*t2*sin( pi*kx-pi/2.*ky-2.*pi*kz)

       do i=0,4,4 ! intra- and inter-layer hopping
        vz(1+i,2+i)=ht1(1)+ht2(1)
        vz(1+i,3+i)=ht1(2)+ht2(2)
        vz(1+i,4+i)=ht1(3)+ht2(3)
        vz(2+i,3+i)=ht1(3)+ht2(3)
        vz(2+i,4+i)=ht1(2)+ht2(2)
        vz(3+i,4+i)=ht1(1)+ht2(1)
       enddo
       endif

       ! HERMITIAN OPERATORS
       do i=1,num_wann
        do j=1,num_wann
         if(j.gt.i) then
          !ham(j,i) = conjg(ham(i,j))
          if(present(vx)) vx(j,i) = conjg(vx(i,j))
          if(present(vy)) vy(j,i) = conjg(vy(i,j))
          if(present(vz)) vz(j,i) = conjg(vz(i,j))
         endif
        enddo
       enddo


      end subroutine init_ham


      end module m_setham_FeMn
