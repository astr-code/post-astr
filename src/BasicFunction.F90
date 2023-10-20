!+---------------------------------------------------------------------+
!| This module contains some basic functions.                          |
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2020-04-18                                     |
!+---------------------------------------------------------------------+
module basicfunction
  !
  implicit none
  !
  Interface integration
    module procedure integration_real8
    module procedure integration_2d
  end Interface integration
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| The subroutine is used to do a digital integration.               |
  !+-------------------------------------------------------------------+
  !| Writen by Jian Fang, 2020-04-17                                   |
  !+-------------------------------------------------------------------+
  function integration_real8(vin) result(vsum)
    !
    ! argument
    real(8) :: vsum
    real(8),intent(in) :: vin(:)
    !
    ! local data
    integer :: dim,i
    !
    dim=size(vin)-1
    !
    vsum=0.d0
    if(dim==0) then
      vsum=vin(1)
    else
      do i=1,dim
        vsum=vsum+vin(i)
      enddo
    endif
    !
  end function integration_real8
  !
  function integration_2d(vin) result(vsum)
    !
    ! argument
    real(8),intent(in) :: vin(:,:)
    real(8) :: vsum
    !
    ! local data
    integer :: dim1,dim2,i,j
    !
    dim1=size(vin,1)-1
    dim2=size(vin,2)-1
    !
    vsum=0.d0
    do j=1,dim2
    do i=1,dim1
      vsum=vsum+vin(i,j)
    enddo
    enddo
    !
  end function integration_2d
  !+-------------------------------------------------------------------+
  !| The end of the function integration_real8.                        |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate second invariant.            |
  !| ref: J. Jeong, and F. Hussain, On the identification of a vortex, |
  !|                                    J. Fluid Mech. 285, 69, (1995).|
  !+-------------------------------------------------------------------+
  function qcrition(du1,du2,du3)
    !
    use commvardefine,only: im,jm,km
    !
    real(8) :: qcrition(0:im,0:jm,0:km)
    !
    real(8),intent(in) :: du1(1:3,0:im,0:jm,0:km),                      &
                          du2(1:3,0:im,0:jm,0:km),                      &
                          du3(1:3,0:im,0:jm,0:km)
    !
    integer :: i,j,k
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      qcrition(i,j,k)=du1(1,i,j,k)*du2(2,i,j,k)+du1(1,i,j,k)*du3(3,i,j,k)  &
                     +du2(2,i,j,k)*du3(3,i,j,k)-du1(2,i,j,k)*du2(1,i,j,k)  &
                     -du1(3,i,j,k)*du3(1,i,j,k)-du2(3,i,j,k)*du3(2,i,j,k)+ &
                      du1(1,i,j,k)**2+du2(2,i,j,k)**2+du3(3,i,j,k)**2      &
                    +(du1(2,i,j,k)*du2(1,i,j,k)+du1(3,i,j,k)*du3(1,i,j,k)  &
                     +du2(3,i,j,k)*du3(2,i,j,k))*2.d0
    end do
    end do
    end do
    !
    return
    !
  end function qcrition
  !+-------------------------------------------------------------------+
  !| The end of the function qcrition.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate swirling strengh.            |
  !| ref1: Pirozzoli, S., M. Bernardini, and F. Grasso,                |
  !|       Characterization of coherent vortical structures in a       |
  !|       supersonic turbulent boundary layer. J. Fluid Mech., 2008,  |
  !|        613: p. 205–231.                                           |
  !| ref2:  Wang, L. and X. Lu, Flow topology in compressible turbulent|
  !|        boundary layer. Journal of Fluid Mechanics, 2012: p.255-278|
  !| https://zh.wikipedia.org/wiki/%E4%B8%89%E6%AC%A1%E6%96%B9%E7%A8%8B|
  !+-------------------------------------------------------------------+
  function lambdaci(du,dv,dw)
    !
    use commvardefine,only: im,jm,km
    !
    real(8) :: lambdaci(0:im,0:jm,0:km)
    !
    real(8),intent(in) :: du(1:3,0:im,0:jm,0:km),                      &
                          dv(1:3,0:im,0:jm,0:km),                      &
                          dw(1:3,0:im,0:jm,0:km)
    !
    integer :: i,j,k
    real(8) :: del,a11,a12,a13,a21,a22,a23,a31,a32,a33,Q,R,var1,var2,  &
               var3,delta
    !
    !$omp parallel default(shared) private(i,j,k,del,a11,a12,a13,a21,  &
    !$omp    a22,a23,a31,a32,a33,q,r,var1,var2,var3,delta)
    !
    !$omp do
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      del=-(du(1,i,j,k)+dv(2,i,j,k)+dw(3,i,j,k))/3.d0
      !
      a11=du(1,i,j,k)+del; a12=du(2,i,j,k);     a13=du(3,i,j,k)
      a21=dv(1,i,j,k);     a22=dv(2,i,j,k)+del; a23=dv(3,i,j,k)
      a31=dw(1,i,j,k);     a32=dw(2,i,j,k);     a33=dw(3,i,j,k)+del

      Q=-0.5d0*(a11*a11+a12*a21+a13*a31+ a21*a12+a22*a22+a23*a32+      &
                a31*a13+a32*a23+a33*a33 )
      R=-1.d0/3.d0*(a11*a11*a11+a11*a12*a21+a11*a13*a31+               &
                    a12*a21*a11+a12*a22*a21+a12*a23*a31+               &
                    a13*a31*a11+a13*a32*a21+a13*a33*a31+               &
                    a21*a11*a12+a21*a12*a22+a21*a13*a32+               &
                    a22*a21*a12+a22*a22*a22+a22*a23*a32+               &
                    a23*a31*a12+a23*a32*a22+a23*a33*a32+               &
                    a31*a11*a13+a31*a12*a23+a31*a13*a33+               &
                    a32*a21*a13+a32*a22*a23+a32*a23*a33+               &
                    a33*a31*a13+a33*a32*a23+a33*a33*a33                )
      !
      delta=Q**3/27.d0+R**2/4.d0
      !
      ! To solve the eqation: lambda**3+Q*lambda+R=0
      if(delta>=1.d-10 ) then
        var1=cube_root(-0.5d0*R+sqrt(delta))
        var2=cube_root(-0.5d0*R-sqrt(delta))
        var3=0.5d0*sqrt(3.d0)*(var1-var2)
        lambdaci(i,j,k)=var3*var3
        !
      else
        lambdaci(i,j,k)=0.d0
      end if
      !
    end do
    end do
    end do
    !$omp end do
    !$omp end parallel
    !
    print*,' ** Swirling strength λ calculated'
    !
  end function lambdaci
  !+-------------------------------------------------------------------+
  !| The end of the function lambdaci.                                 |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate vorticity.                   |
  !+-------------------------------------------------------------------+
  function omega(du,dv,dw,dir)
    !
    use commvardefine,only: im,jm,km
    !
    real(8) :: omega(0:im,0:jm,0:km)
    !
    real(8),intent(in) :: du(1:3,0:im,0:jm,0:km),                      &
                          dv(1:3,0:im,0:jm,0:km),                      &
                          dw(1:3,0:im,0:jm,0:km)
    character(len=1),intent(in) :: dir
    !
    if(dir=='x') then
      omega=dw(2,:,:,:)-dv(3,:,:,:)
    elseif(dir=='y') then
      omega=du(3,:,:,:)-dw(1,:,:,:)
    elseif(dir=='z') then
      omega=dv(1,:,:,:)-du(2,:,:,:)
    elseif(dir=='e') then
      omega=0.5d0*( (dw(2,:,:,:)-dv(3,:,:,:))**2+                      &
                    (du(3,:,:,:)-dw(1,:,:,:))**2+                      &
                    (dv(1,:,:,:)-du(2,:,:,:))**2 )
    elseif(dir=='m') then
      omega=sqrt( (dw(2,:,:,:)-dv(3,:,:,:))**2+                      &
                  (du(3,:,:,:)-dw(1,:,:,:))**2+                      &
                  (dv(1,:,:,:)-du(2,:,:,:))**2 )
    else
      stop ' !! error in defining direction in function omega'
    endif
    !
    print*,' ** vorticity ω_',dir,' calculated'
    !
  end function omega
  !+-------------------------------------------------------------------+
  !| The end of the function omega.                                    |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| The subroutine is to calculate cube root of a real8 variable.     |
  !+-------------------------------------------------------------------+
  real(8) function cube_root (original_value)
    !
     real(8),intent(in) :: original_value
     real(8),parameter :: cube_root_power = 0.333333333333333d0
     if(original_value>=0.d0) then
       cube_root=original_value**cube_root_power
     else
       cube_root=-((-original_value)**cube_root_power)
     endif
     !
  end function cube_root
  !+-------------------------------------------------------------------+
  !| The end of the function cube_root.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate density schlieren.           |
  !+-------------------------------------------------------------------+
  function rnsxy(ro,x,y)
    !
    use commvardefine, only: im,jm
    use gradsolver, only: grad_xy
    !
    ! arguments
    real(8) :: rnsxy(0:im,0:jm)
    real(8),intent(in) :: ro(0:im,0:jm)
    real(8),intent(in),optional :: x(0:im,0:jm),y(0:im,0:jm)
    !
    ! local data
    integer :: i,j
    real(8),allocatable :: dro(:,:,:)
    real(8) ::var1
    real(8) ::rnsmin=0.d0,rnsmax=0.d0,c1=1.d0,c2=10.d0
    save rnsmin,rnsmax
    !
    allocate(dro(1:2,0:im,0:jm))
    !
    if(present(x) .and. present(y)) then
      dro=grad_xy(ro,x,y)
    else
      dro=grad_xy(ro)
    endif
    !
    if(rnsmax<=1.d-10 .and. rnsmin<1.d-10) then
      !
      rnsmax=0.d0
      rnsmin=1.d10
      do j=0,jm
      do i=0,im
        !
        var1=dro(1,i,j)**2+dro(2,i,j)**2
        rnsmax=max(rnsmax,var1)
        rnsmin=min(rnsmin,var1)
        !
      end do
      end do
      rnsmin=sqrt(rnsmin)
      rnsmax=sqrt(rnsmax)
      !
      print*,'rnsmin',rnsmin,'rnsmax',rnsmax
    endif
    !
    rnsxy=c1*exp(-c2*(sqrt(dro(1,:,:)**2+dro(2,:,:)**2)-rnsmin)/       &
                                                       (rnsmax-rnsmin))
    !
    deallocate(dro)
    !
    return
    !
  end function rnsxy
  !+-------------------------------------------------------------------+
  !| The end of the function rnsxy.                                    |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate density schlieren.           |
  !+-------------------------------------------------------------------+
  function rns3d(ro,droin)
    !
    use commvardefine, only: im,jm,km
    use gradsolver, only: grad_3d
    !
    ! arguments
    real(8) :: rns3d(0:im,0:jm,0:km)
    real(8),intent(in) :: ro(0:im,0:jm,0:km)
    real(8),intent(in),optional :: droin(1:3,0:im,0:jm,0:km)
    !
    ! local data
    integer :: i,j,k
    real(8),allocatable :: dro(:,:,:,:)
    real(8) ::c1,c2,var1
    real(8),save :: rnsmin=-1.d0,rnsmax=-1.d0
    !
    allocate(dro(1:3,0:im,0:jm,0:km))
    !
    if(present(droin)) then
      dro=droin
    else
      print*,' ** calculating density gradient'
      dro=grad_3d(ro)
    endif
    !
    c1=1.d0
    c2=1.d0
    ! !
    ! if(rnsmax<=1.d-10 .and. rnsmin<1.d-10) then
    !   rnsmax=0.d0
    !   rnsmin=1.d10
    !   do k=0,km
    !   do j=0,jm
    !   do i=0,im
    !     !
    !     var1=sqrt(dro(1,i,j,k)**2+dro(2,i,j,k)**2+dro(3,i,j,k)**2)
    !     !
    !     rnsmax=max(rnsmax,var1)
    !     rnsmin=min(rnsmin,var1)
    !     !
    !     rns3d=var1
    !     !
    !   end do
    !   end do
    !   end do
    !   print*,'rnsmin',rnsmin,'rnsmax',rnsmax
    ! endif
    !
    print*,' ** calculating rns ... '
    rns3d=sqrt(dro(1,:,:,:)**2+dro(2,:,:,:)**2+dro(3,:,:,:)**2)
    !
    if(rnsmax<0.d0) rnsmax=maxval(rns3d)
    if(rnsmin<0.d0) rnsmin=minval(rns3d)
    !
    rns3d=c1*exp(-c2*(rns3d-rnsmin)/(rnsmax-rnsmin))
    !
    print*,' ** calculating rns ... done '
    !
    deallocate(dro)
    !
    return
    !
  end function rns3d
  !+-------------------------------------------------------------------+
  !| The end of the function rnsxy.                                    |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !|calculate dynamic viscosity coefficient at different temperature.  |
  !+-------------------------------------------------------------------+
  real(8) function miucal(temper)
    !
    use commvardefine,only: ref_t,tempconst,tempconst1
    !
    real(8),intent(in) :: temper
    !
    ! temper represent temperature, dimensionless
    ! below calculate miucal using sutherland's law
    ! tempconst=110.4d0/ref_t
    ! tempconst1=1.d0+tempconst
    !
    tempconst=110.4d0/ref_t
    tempconst1=1.d0+tempconst
    !
    miucal=temper*sqrt(temper)*tempconst1/(temper+tempconst)
    !
    return
  end function MiuCal
  !+-------------------------------------------------------------------+
  ! End of function MiuCal.                                            |
  !+-------------------------------------------------------------------+
  !!
  subroutine uplus_yplus(uplus,yplus,u,y,ro,twall,utaw,Re)
    !
    use gradsolver,only: grad_y
    !
    real(8),intent(out),allocatable :: uplus(:),yplus(:)
    real(8),intent(out) :: utaw
    real(8),intent(in) :: u(:),y(:)
    real(8),intent(in) :: Re
    real(8),optional,intent(in) :: ro(:),twall
    !
    integer :: dim
    real(8) :: miu
    real(8),allocatable :: dudy(:)
    !
    dim=size(u)
    !
    allocate(dudy(1:dim))
    dudy=grad_y(u,y)
    !
    ! print*,y
    ! print*,u
    ! print*,dudy
    !
    if(present(ro) .and. present(twall)) then
      ! compressible profile
      miu=miucal(twall)/Re
    else
      ! incompressible profile
    endif
    !
  end subroutine uplus_yplus
  !!
  !+-------------------------------------------------------------------+
  !| this subroutine is to calculate uplus-yplse                       |
  !+-------------------------------------------------------------------+
  subroutine upluscal(uplus,yplus,u,y,ro,tw,taw,utaw)
    !
    use commvardefine,only: Reynolds
    use gradsolver,only: grad_y
    !
    ! arguments
    real(8),intent(out),allocatable :: uplus(:),yplus(:)
    real(8),intent(out),optional :: utaw
    real(8),intent(in) :: u(0:),y(0:),ro(0:)
    real(8),intent(in) :: tw
    real(8),intent(in),optional :: taw
    !
    ! local data
    integer :: dim,j,j1
    real(8),allocatable :: dudy(:),uvd(:),yin(:)
    real(8) :: miu,utaw_local,var1,var2,tawx
    !
    dim=size(u)-1
    !
    miu=miucal(tw)/Reynolds
    !
    if(present(taw)) then
      tawx=taw
    else
      allocate(dudy(0:dim))
      !
      allocate(yin(0:dim))
      yin=y
      !
      dudy=grad_y(u,yin)
      !
      tawx=dabs(miu*dudy(0))
      ! 
    endif
    !
    utaw_local=dsqrt(tawx/ro(0))
    if(present(utaw)) utaw=utaw_local
    !
    allocate(uplus(0:dim),yplus(0:dim),uvd(0:dim))
    !
    uvd(0)=0.d0
    do j=1,dim
      uvd(j)=0.d0
      do j1=1,j
        var1=dsqrt(0.5d0*(ro(j1)+ro(j1-1))/ro(0))
        var2=u(j1)-u(j1-1)
        uvd(j)=uvd(j)+var1*var2
      end do
    end do
    !
    do j=0,dim
      yplus(j)=ro(0)*utaw_local*yin(j)/miu
      uplus(j)=uvd(j)/utaw_local
    end do
    !
    if(allocated(dudy)) deallocate(dudy)
    
    deallocate(uvd)
    !
    write(*,'(A)')'  ------ first point from wall ------'
    write(*,'(A)')'                y1+               u1+'
    write(*,'(1x,F18.7,F18.7)')yplus(1),uplus(1)
    write(*,'(A)')'  -----------------------------------'
    print*,' ** tawx=',tawx
    print*,' ** utaw=',utaw_local
    print*,' ** uplus-yplus calculated.'
    !
  end subroutine upluscal
  !+-------------------------------------------------------------------+
  ! End of subroutine upluscal.                                        |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is to calculate uplus-yplse                       |
  !+-------------------------------------------------------------------+
  subroutine upluscal_imc(uplus,yplus,u,y,taw,utaw,Re)
    !
    use gradsolver,only: grad_y
    !
    ! arguments
    real(8),intent(out),allocatable :: uplus(:),yplus(:)
    real(8),intent(out),optional :: utaw
    real(8),intent(in) :: u(0:),y(0:)
    real(8),intent(in) :: Re
    real(8),intent(in),optional :: taw
    !
    ! local data
    integer :: dim,j,j1
    real(8),allocatable :: dudy(:)
    real(8) :: miu,utaw_local,var1,var2,tawx
    !
    dim=size(u)-1
    !
    miu=1.d0/Re
    !
    if(present(taw)) then
      tawx=taw
    else
      allocate(dudy(0:dim))
      !
      dudy=grad_y(u,y)
      !
      tawx=abs(miu*dudy(0))
      !
    endif
    !
    utaw_local=sqrt(tawx)
    if(present(utaw)) utaw=utaw_local
    !
    allocate(uplus(0:dim),yplus(0:dim))
    !
    do j=0,dim
      yplus(j)=utaw_local*y(j)/miu
      uplus(j)=u(j)/utaw_local
    end do
    !
    if(allocated(dudy)) deallocate(dudy)
    
    write(*,'(A)')'  ------ first point from wall ------'
    write(*,'(A)')'                y1+               u1+'
    write(*,'(1x,F18.7,F18.7)')yplus(1),uplus(1)
    write(*,'(A)')'  -----------------------------------'
    print*,' ** tawx=',tawx
    print*,' ** utaw=',utaw_local
    print*,' ** uplus-yplus calculated.'
    !
  end subroutine upluscal_imc
  !+-------------------------------------------------------------------+
  ! End of subroutine upluscal.                                        |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This function is the the inverse function of the function tanh(x).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function argtanh(var)
    !
    real(8) :: var,argtanh
    !
    argtanh=0.5d0*dlog((1.d0+var)/(1.d0-var))
    !
    return
    !
  end function
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function argtanh.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  real(8) function expfun(xmin,ymin,xmax,ymax,x)
    !
    real(8),intent(in) :: xmin,ymin,xmax,ymax,x
    !
    real(8) :: a,b,var1
    !
    b=ymin
    var1=log(abs(ymax/ymin))
    var1=var1/(xmax-xmin)
    a=exp(var1)
    !
    expfun=b*a**(x-xmin)
    !
  end function expfun
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function expfun.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  real(8) function sqrtn(var,n)
    !
    real(8) :: var
    integer :: n
    !
    integer :: m,i
    !
    if(mod(n,2)==0) then
      sqrtn=sqrt(abs(var))
      !
      m=n/2
      do i=2,m
        sqrtn=sqrt(sqrtn)
      enddo
      !
    else
      stop ' cant do odd n now'
    endif
    !
  end function sqrtn
  !
  real(8) function avec(a)
    !
    use commvardefine,only:pi
    !
    real(8) :: a(2)
    !
    real(8) :: var1,var2
    !
    var1=a(1)
    var2=sqrt(a(1)**2+a(2)**2)
    !
    !
    if(a(2)>=0.d0) then
      avec=acos(var1/var2)
    else
      avec=2.d0*pi-acos(var1/var2)
    endif
    !
  end function avec
  !
  !+---------------------------------------------------------+
  !| This fundtion is a general f(x)=ftarget solver using    |
  !| Newton's method                                         |
  !+---------------------------------------------------------+
  !| writen by J. Fang, 2021-12-26.                          |
  !+---------------------------------------------------------+
  function fxsolver(xmin,xmax,ftarget,fx) result(x)
    !
    use interpolation, only : linear1d
    !
    real(8),intent(in) :: xmin,xmax,ftarget
    real(8) :: x
    !
    real(8) :: fx1,fx2,xans,fxans,x1,x2
    real(8) :: error,ethres
    integer :: counter
    !
    interface
      !
      function fx(xx)
        real(8) :: xx,fx
      end function fx
      !
    end interface
    !
    ethres=1.d-10
    !
    x1=xmin
    x2=xmax
    fx1=fx(x1)
    fx2=fx(x2)
    !
    error=1.d0
    counter=0
    !
    write(*,'(A4,1X,A15)')'n','error'
    do while(abs(error)>ethres .and. counter<1000)
      !
      xans=linear1d(fx1,fx2,x1,x2,ftarget)
      fxans=fx(xans)
      !
      ! to ensure the xmin,xmax is within the effective/monotonicity range
      if(ftarget>=fx1 .and. ftarget<=fx2) then
        !
        ! print*,' ** xans',xans,fxans
        !
        if(fxans>=ftarget) then
          x2 =xans
          fx2=fxans
        elseif(fxans<=ftarget) then
          x1 =xans
          fx1=fxans
        endif
        !
      elseif(ftarget<=fx1 .and. ftarget>=fx2) then
        ! print*,' ** xans',xans,fxans
        !
        if(fxans>=ftarget) then
          x1 =xans
          fx1=fxans
        elseif(fxans<=ftarget) then
          x2 =xans
          fx2=fxans
        endif
        !
      else
        print*,' !! xmin,xmax not in the effective range '
        print*,' ** fx1=',fx1
        print*,' ** fx2=',fx2
        print*,' ** ftarget=',ftarget
        !
        stop
        !
      endif
      !
      counter=counter+1
      error=fxans-ftarget
      !
      write(*,'(I4,1X,E15.7E3)')counter,error
      !
    enddo
    !
    print*,' ** answer found x=',xans,'fx=',fxans
    !
    x=xans
    !
    return
    !
  end function fxsolver
  !
  subroutine quadfit(x1,y1,x2,y2,x3,y3,xh,yh)
    !
    real(8),intent(in) :: x1,y1,x2,y2,x3,y3
    real(8),intent(out) :: xh,yh
    !
    real(8) :: a(3),mat(3,3),y(3)
    !
    mat(1,1)=x1**2
    mat(1,2)=x1
    mat(1,3)=1.d0
    !
    mat(2,1)=x2**2
    mat(2,2)=x2
    mat(2,3)=1.d0
    !
    mat(3,1)=x3**2
    mat(3,2)=x3
    mat(3,3)=1.d0
    !
    y(1)=y1
    y(2)=y2
    y(3)=y3
    !
    a=matmul(matinv3(mat),y)
    !
    xh=-0.5d0*a(2)/a(1)
    yh=-0.25d0*(a(2)**2-4.d0*a(1)*a(3))/a(1)
    !
    return
    !
  end subroutine quadfit
  !
  pure function matinv3(A) result(B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(8), intent(in) :: A(3,3)   !! Matrix
    real(8)             :: B(3,3)   !! Inverse matrix
    real(8)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1.d0/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                 - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                 + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    !
    return
    !
  end function matinv3
  !
  ! ref: https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter#Moving_average
  function SavitzkyGolayFilter(f) result(ff)
    !
    real(8),intent(in) :: f(:)
    real(8) :: ff(1:size(f))
    !
    integer :: i,m
    !
    m=size(ff)
    !
    ff=f
    do i=3,m-2
      ff(i)=-3.d0*f(i-2)+12.d0*f(i-1)+17.d0*f(i)+12.d0*f(i+1)-3.d0*f(i+2)
      ff(i)=ff(i)/35.d0
    enddo
    !
    return
    !
  end function SavitzkyGolayFilter
  !
  function spafilter10_basic(f) result(ff)
    !
    use singleton
    !
    ! arguments
    real(8),intent(in) :: f(0:)
    real(8) :: ff(0:size(f)-1)
    !
    ! local data
    real(8) :: b(0:size(f)-1)
    integer :: k,dim
    real(8) :: aff(3)
    real(8),allocatable :: fc(:,:)
    !
    aff(1)=0.49
    aff(2)=0.49
    aff(3)=0.49
    !
    dim=size(f)-1
    !
    call ptdsfilter_ini(fc,dim)
    !
    b =PFilterRHS2(f)
    ff=ptds_cal(b,aff,fc)
    !
    return
    !
  end function spafilter10_basic
  !
  subroutine ptdsfilter_ini(cc,dim)
    !
    integer,intent(in) :: dim
    real(8),allocatable,intent(out) :: cc(:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af: input dat
    ! cc: output array
    ! dim: input dat
    ! l: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! local data
    integer :: l
    real(8) :: af,af_1,af_2,af_3
    !
    af=0.49d0
    !
    allocate(cc(1:2,0:dim))
    !
    if(dim==0) then
      cc=0.d0
      return
    endif
    !
    af_1=af
    af_2=af
    af_3=af
    !
    cc(1,0)=1.d0
    cc(2,0)=af_1/cc(1,0)
    !
    cc(1,1)=1.d0-af_2*cc(2,0)
    cc(2,1)=af_2/cc(1,1)
    !
    cc(1,2)=1.d0-af_3*cc(2,1)
    !
    do l=2,dim-3
      cc(2,l)=af_3/cc(1,l)
      cc(1,l+1)=1.d0-af_3*cc(2,l)
    end do
    cc(2,dim-2)=af_3/cc(1,dim-2)
    !
    cc(1,dim-1)=1.d0-af_2*cc(2,dim-2)
    cc(2,dim-1)=af_2/cc(1,dim-1)
    !
    cc(1,dim)=1.d0-af_1*cc(2,dim-1)
    !
    cc(1,:)=1.d0/cc(1,:)
    !
    return
    !
  end subroutine ptdsfilter_ini
  !
  subroutine compactfilter(v,vf)
    !
    real(8),intent(in) :: v(0:)
    real(8),intent(out) :: vf(0:size(v)-1)
    !
  end subroutine compactfilter
  !
  function ptds_cal(bd,af,cc) result(xd)
    !
    real(8),intent(in) :: af(3)
    real(8),intent(in) :: bd(0:),cc(1:,0:)
    real(8) :: xd(0:size(bd)-1)
    !
    ! local data
    integer :: l,dim
    real(8),allocatable,dimension(:) :: yd
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af(3): input dat
    ! bd: input array
    ! xd: output array
    ! cc: input array
    ! dim: input dat
    ! l, yd: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    dim=size(bd)-1
    !
    allocate ( yd(0:dim) )
    !
    ! the block with boundary at i=0 and i=im
    !
    yd(0)=bd(0)*cc(1,0)
    yd(1)=(bd(1)-af(2)*yd(0))*cc(1,1)
    do l=2,dim-2
      yd(l)=(bd(l)-af(3)*yd(l-1))*cc(1,l)
    end do
    yd(dim-1)=(bd(dim-1)-af(2)*yd(dim-2))*cc(1,dim-1)
    yd(dim)=(bd(dim)-af(1)*yd(dim-1))*cc(1,dim)
    !
    xd(dim)=yd(dim)
    do l=dim-1,0,-1
      xd(l)=yd(l)-cc(2,l)*xd(l+1)
    end do
    !
    deallocate( yd )
    !
    return
    !
  end function ptds_cal
  !
  function pfilterrhs2(var) result(b)
    !
    real(8),intent(in) :: var(0:)
    real(8) :: b(0:size(var)-1)
    !
    ! local data
    integer :: l,dim
    real(8) :: var0,var1,var2,var3,var4,var5,alfa
    real(8) :: coef6i(0:3),coef8i(0:4),coef10i(0:5),coefb(0:4,0:8)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! dim: dimensions in the direction.
    ! var: the variable of the filetering.
    ! var* : temporary variable.
    ! b: return variable.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    dim=size(var)-1
    alfa=0.49d0
    !
    ! 6th-order
    coef6i(0)=(11.d0 +10.d0*alfa)   /32.d0
    coef6i(1)=(15.d0 +34.d0*alfa)   /64.d0
    coef6i(2)=(-3.d0 + 6.d0*alfa)   /32.d0
    coef6i(3)=( 1.d0 - 2.d0*alfa)   /64.d0
    !
    ! 8th-order
    coef8i(0)=(93.d0 +70.d0*alfa)   /256.d0
    coef8i(1)=( 7.d0 +18.d0*alfa)   /32.d0
    coef8i(2)=(-7.d0 +14.d0*alfa)   /64.d0
    coef8i(3)=( 1.d0 - 2.d0*alfa)   /32.d0
    coef8i(4)=(-1.d0 + 2.d0*alfa)   /256.d0
    !
    ! 10th-order
    coef10i(0)=(193.d0 +126.d0*alfa)/512.d0
    coef10i(1)=(105.d0 +302.d0*alfa)/512.d0
    coef10i(2)=(-15.d0 + 30.d0*alfa)/128.d0
    coef10i(3)=( 45.d0 - 90.d0*alfa)/1024.d0
    coef10i(4)=( -5.d0 + 10.d0*alfa)/512.d0
    coef10i(5)=(  1.d0 -  2.d0*alfa)/1024.d0
    !
    ! boundary coef10ifficent
    !
    ! 8th-order for the point3
    coefb(3,0)=(  1.d0 - 2.d0*alfa)  /256.d0
    coefb(3,1)=( -1.d0 + 2.d0*alfa)  /32.d0
    coefb(3,2)=(  7.d0 +50.d0*alfa)  /64.d0
    coefb(3,3)=( 25.d0 +14.d0*alfa)  /32.d0
    coefb(3,4)=( 35.d0 +58.d0*alfa)  /128.d0
    coefb(3,5)=( -7.d0 +14.d0*alfa)  /32.d0
    coefb(3,6)=(  7.d0 -14.d0*alfa)  /64.d0
    coefb(3,7)=( -1.d0 + 2.d0*alfa)  /32.d0
    coefb(3,8)=(  1.d0 - 2.d0*alfa)  /256.d0
    !
    ! 6-order asymmetry scheme for point 2
    coefb(2,0)=( -1.d0 + 2.d0*alfa)  /64.d0
    coefb(2,1)=(  3.d0 +26.d0*alfa)  /32.d0
    coefb(2,2)=( 49.d0 +30.d0*alfa)  /64.d0
    coefb(2,3)=(  5.d0 + 6.d0*alfa)  /16.d0
    coefb(2,4)=(-15.d0 +30.d0*alfa)  /64.d0
    coefb(2,5)=(  3.d0 - 6.d0*alfa)  /32.d0
    coefb(2,6)=( -1.d0 + 2.d0*alfa)  /64.d0
    coefb(2,7)=0.d0
    coefb(2,8)=0.d0
    !
    ! 6-order asymmetry scheme for point 1
    coefb(1,0)=( 1.d0 +62.d0*alfa)  /64.d0
    coefb(1,1)=(29.d0 + 6.d0*alfa)  /32.d0
    coefb(1,2)=(15.d0 +34.d0*alfa)  /64.d0
    coefb(1,3)=(-5.d0 +10.d0*alfa)  /16.d0
    coefb(1,4)=(15.d0 -30.d0*alfa)  /64.d0
    coefb(1,5)=(-3.d0 + 6.d0*alfa)  /32.d0
    coefb(1,6)=( 1.d0 - 2.d0*alfa)  /64.d0
    coefb(1,7)=0.d0
    coefb(1,8)=0.d0
    !
    ! 6-order asymmetry scheme for point 0
    coefb(0,0)=( 63.d0 + 1.d0*alfa)  /64.d0
    coefb(0,1)=(  3.d0 +29.d0*alfa)  /32.d0
    coefb(0,2)=(-15.d0 +15.d0*alfa)  /64.d0
    coefb(0,3)=(  5.d0 - 5.d0*alfa)  /16.d0
    coefb(0,4)=(-15.d0 +15.d0*alfa)  /64.d0
    coefb(0,5)=(  3.d0 - 3.d0*alfa)  /32.d0
    coefb(0,6)=( -1.d0 + 1.d0*alfa)  /64.d0
    coefb(0,7)=0.d0
    coefb(0,8)=0.d0
    !
    ! the block with boundary at i=0 and i=im
    !
    b(0)=coefb(0,0)*var(0)+coefb(0,1)*var(1)+coefb(0,2)*var(2)+      &
         coefb(0,3)*var(3)+coefb(0,4)*var(4)+coefb(0,5)*var(5)+      &
         coefb(0,6)*var(6)
    !
    b(1)=coefb(1,0)*var(0)+coefb(1,1)*var(1)+coefb(1,2)*var(2)+      &
         coefb(1,3)*var(3)+coefb(1,4)*var(4)+coefb(1,5)*var(5)+      &
         coefb(1,6)*var(6)
    !
    b(2)=coefb(2,0)*var(0)+coefb(2,1)*var(1)+coefb(2,2)*var(2)+      &
         coefb(2,3)*var(3)+coefb(2,4)*var(4)+coefb(2,5)*var(5)+      &
         coefb(2,6)*var(6)
    !
    var0=var(3)+var(3)
    var1=var(4)+var(2)
    var2=var(5)+var(1)
    var3=var(6)+var(0)
    b(3)=coef6i(0)*var0+coef6i(1)*var1+coef6i(2)*var2+coef6i(3)*var3
    !
    var0=var(4)+var(4)
    var1=var(5)+var(3)
    var2=var(6)+var(2)
    var3=var(7)+var(1)
    var4=var(8)+var(0)
    b(4)=coef8i(0)*var0+coef8i(1)*var1+coef8i(2)*var2 +              &
         coef8i(3)*var3+coef8i(4)*var4
    !
    ! inner block
    do l=5,dim-5
      !
      var0=var(l)+var(l)
      var1=var(l+1)+var(l-1)
      var2=var(l+2)+var(l-2)
      var3=var(l+3)+var(l-3)
      var4=var(l+4)+var(l-4)
      var5=var(l+5)+var(l-5)
      !
      b(l)=coef10i(0)*var0+coef10i(1)*var1+coef10i(2)*var2+          &
           coef10i(3)*var3+coef10i(4)*var4+coef10i(5)*var5
      !
    end do
    !
    var0=var(dim-4)+var(dim-4)
    var1=var(dim-3)+var(dim-5)
    var2=var(dim-2)+var(dim-6)
    var3=var(dim-1)+var(dim-7)
    var4=var(dim)  +var(dim-8)
    b(dim-4)=coef8i(0)*var0+coef8i(1)*var1+coef8i(2)*var2 +          &
             coef8i(3)*var3+coef8i(4)*var4
    !
    var0=var(dim-3)+var(dim-3)
    var1=var(dim-2)+var(dim-4)
    var2=var(dim-1)+var(dim-5)
    var3=var(dim)  +var(dim-6)
    b(dim-3)=coef6i(0)*var0+coef6i(1)*var1+                          &
             coef6i(2)*var2+coef6i(3)*var3
    !
    b(dim-2)=coefb(2,0)*var(dim)  +coefb(2,1)*var(dim-1)+            &
             coefb(2,2)*var(dim-2)+coefb(2,3)*var(dim-3)+            &
             coefb(2,4)*var(dim-4)+coefb(2,5)*var(dim-5)+            &
             coefb(2,6)*var(dim-6)
    !
    b(dim-1)=coefb(1,0)*var(dim)  +coefb(1,1)*var(dim-1)+            &
             coefb(1,2)*var(dim-2)+coefb(1,3)*var(dim-3)+            &
             coefb(1,4)*var(dim-4)+coefb(1,5)*var(dim-5)+            &
             coefb(1,6)*var(dim-6)
    !
    b(dim)=coefb(0,0)*var(dim)  +coefb(0,1)*var(dim-1)+            &
           coefb(0,2)*var(dim-2)+coefb(0,3)*var(dim-3)+            &
           coefb(0,4)*var(dim-4)+coefb(0,5)*var(dim-5)+            &
           coefb(0,6)*var(dim-6)
    !
    !
  end function pfilterrhs2
  !
end module basicfunction
!+---------------------------------------------------------------------+
!| The end of the module basicfunction.                                |
!+---------------------------------------------------------------------+