!+---------------------------------------------------------------------+
!| This module contains some common functions                          |
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2020-04-17                                     |
!+---------------------------------------------------------------------+
module gradsolver
  !
  implicit none
  !
  Interface grad
    !
    module procedure grad_3d
    module procedure grad_xy
    module procedure grad_y
    !
  end Interface grad
  !
  real(8),allocatable :: cci(:,:),ccj(:,:),cck(:,:),afi(:),afj(:),afk(:)
  !+----------------+--------------------------------------------------+
  !|     cci,ccj,cck| parameters for gradient calculator.              |
  !|     afi,afj,afk|                                                  |
  !+----------------+--------------------------------------------------+
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !|this function is used to calculate field gradient.                 |
  !+-------------------------------------------------------------------+
  !|Writen by Jian Fang , 2020-04-18.                                  |
  !+-------------------------------------------------------------------+
  function grad_3d(f,x,y,z) result(df)
    !
    use commvardefine, only: lihomo,ljhomo,lkhomo,im,jm,km
    !  
    ! arguments
    real(8) :: df(1:3,0:im,0:jm,0:km)
    real(8),intent(in) :: f(0:im,0:jm,0:km)
    real(8),intent(in),optional :: x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),&
                                   z(0:im,0:jm,0:km)
    !

    ! local data
    integer :: i,j,k
    real(8),allocatable :: vtemp(:)
    real(8),allocatable :: ddi(:,:,:,:),ddj(:,:,:,:),ddk(:,:,:,:)
    save ddi,ddj,ddk
    !$ SAVE vtemp
    !$OMP THREADPRIVATE(vtemp)
    !
    if(.not.allocated(afi)) then
      allocate(afi(0:im),cci(1:2,0:im))
      call tds_ini(afi,cci,im,lihomo)
    endif
    if(.not.allocated(afj)) then
      allocate(afj(0:jm),ccj(1:2,0:jm))
      call tds_ini(afj,ccj,jm,ljhomo)
    endif
    if(.not.allocated(afk)) then
      allocate(afk(0:km),cck(1:2,0:km))
      call tds_ini(afk,cck,km,lkhomo)
    endif
    !
    if(.not. allocated(ddi)) then
      allocate(ddi(1:3,0:im,0:jm,0:km),ddj(1:3,0:im,0:jm,0:km),        &
               ddk(1:3,0:im,0:jm,0:km) )
      !
      if(present(x) .and. present(y) .and. present(z)) then
        call gridjacobian_3d(x,y,z,ddi,ddj,ddk)
      else
        stop ' !! x,y,z needed in function grad_3d !!'
      endif
      !
    else
      ! print*,' ** Use the saved grid jacobian matrix'
    endif
    !
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
    !
    allocate(vtemp(0:im))
    !$OMP DO
    do k=0,km
    do j=0,jm
      vtemp=dfdi(f(:,j,k))
      df(1,:,j,k)=vtemp(:)*ddi(1,:,j,k)
      df(2,:,j,k)=vtemp(:)*ddi(2,:,j,k)
      df(3,:,j,k)=vtemp(:)*ddi(3,:,j,k)
    end do
    end do
    !$OMP END DO
    deallocate(vtemp)
    !
    allocate(vtemp(0:jm))
    !$OMP DO
    do k=0,km
    do i=0,im
      vtemp=dfdj(f(i,:,k))
      df(1,i,:,k)=df(1,i,:,k)+vtemp(:)*ddj(1,i,:,k)
      df(2,i,:,k)=df(2,i,:,k)+vtemp(:)*ddj(2,i,:,k)
      df(3,i,:,k)=df(3,i,:,k)+vtemp(:)*ddj(3,i,:,k)
    end do
    end do
    !$OMP END DO
    deallocate(vtemp)
    !
    allocate(vtemp(0:km))
    !$OMP DO
    do j=0,jm
    do i=0,im
      vtemp=dfdk(f(i,j,:))
      df(1,i,j,:)=df(1,i,j,:)+vtemp(:)*ddk(1,i,j,:)
      df(2,i,j,:)=df(2,i,j,:)+vtemp(:)*ddk(2,i,j,:)
      df(3,i,j,:)=df(3,i,j,:)+vtemp(:)*ddk(3,i,j,:)
    end do
    end do
    !$OMP END DO
    deallocate(vtemp)
    !
    !$OMP END PARALLEL
    !
    print*,' ** 3D gradient field calculated'
    !
  end function grad_3d
  !
  function grad_xy(f,x,y) result(df)
    !
    use commvardefine, only: lihomo,ljhomo,im,jm
    !  
    ! arguments
    real(8) :: df(1:2,0:im,0:jm)
    real(8),intent(in) :: f(0:im,0:jm)
    real(8),intent(in),optional :: x(0:im,0:jm),y(0:im,0:jm)
    !
    ! local data
    integer :: i,j
    real(8),allocatable :: vtemp(:)
    real(8),allocatable :: ddi(:,:,:),ddj(:,:,:)
    save ddi,ddj
    !
    if(.not.allocated(afi)) then
      allocate(afi(0:im),cci(1:2,0:im))
      call tds_ini(afi,cci,im,lihomo)
    endif
    if(.not.allocated(afj)) then
      allocate(afj(0:jm),ccj(1:2,0:jm))
      call tds_ini(afj,ccj,jm,ljhomo)
    endif
    !
    if(.not. allocated(ddi)) then
      allocate(ddi(1:2,0:im,0:jm),ddj(1:2,0:im,0:jm) )
      !
      if(present(x) .and. present(y)) then
        call gridjacobian_xy(x,y,ddi,ddj)
      else
        stop ' !! x,y needed in function grad_xy !!'
      endif
      !
    else
      ! print*,' ** Use the saved grid jacobian matrix'
    endif
    !
    allocate(vtemp(0:im))
    do j=0,jm
      vtemp=dfdi(f(:,j))
      df(1,:,j)=vtemp(:)*ddi(1,:,j)
      df(2,:,j)=vtemp(:)*ddi(2,:,j)
    end do
    deallocate(vtemp)
    !
    allocate(vtemp(0:jm))
    do i=0,im
      vtemp=dfdj(f(i,:))
      df(1,i,:)=df(1,i,:)+vtemp(:)*ddj(1,i,:)
      df(2,i,:)=df(2,i,:)+vtemp(:)*ddj(2,i,:)
    end do
    deallocate(vtemp)
    !
    ! print*,' ** Spatial gradient at x-y plane calculated'
    !
  end function grad_xy
  !+-------------------------------------------------------------------+
  ! end of the function grad_xy.                                       |
  !+-------------------------------------------------------------------+
  !
  function grad_x(f,y) result(df)
    !
    use commvardefine, only: lihomo,im
    !  
    ! arguments
    real(8) :: df(0:im)
    real(8),intent(in) :: f(0:im)
    real(8),intent(in),optional :: y(0:im)
    !
    ! local data
    integer :: i
    real(8),allocatable :: ddi(:)
    save ddi
    !
    if(.not.allocated(afi)) then
      allocate(afi(0:im),cci(1:2,0:im))
      call tds_ini(afi,cci,im,lihomo)
    endif
    !
    if(.not. allocated(ddi)) allocate(ddi(0:im) )
    if( present(y)) ddi=1.d0/dfdi(y)
    ! !
    ! if(.not. allocated(ddi)) then
    !   allocate(ddi(0:im) )
    !   if( present(y)) then
    !     ddi=1.d0/dfdi(y)
    !   else
    !     stop ' !! y needed in function grad_y !!'
    !   endif
    !   !
    ! else
    !   print*,' ** Use the saved grid iacobian matrix'
    ! endif
    !
    df=dfdi(f)*ddi
    !
    print*,' ** Spatial gradient along x calculated'
    !
  end function grad_x
  !+-------------------------------------------------------------------+
  ! end of the function grad_y.                                        |
  !+-------------------------------------------------------------------+
  !
  function grad_y(f,y) result(df)
    !
    use commvardefine, only: ljhomo
    !  
    ! arguments
    real(8),intent(in) :: f(0:)
    real(8),intent(in),optional :: y(0:)
    real(8) :: df(0:size(f)-1)
    !
    ! local data
    integer :: jm,j
    real(8),allocatable :: ddj(:)
    save ddj
    !
    jm=size(f)-1
    !
    if(allocated(afj) .and. size(afj).ne.jm+1) then
      deallocate(afj,ccj)
    endif
    !
    if(.not. allocated(afj)) allocate(afj(0:jm),ccj(1:2,0:jm))
    !
    call tds_ini(afj,ccj,jm,ljhomo)
    !
    if(.not. allocated(ddj)) allocate(ddj(0:jm) )
    if( present(y)) ddj=1.d0/dfdj(y)
    !
    ! if(.not. allocated(ddj)) then
    !   allocate(ddj(0:jm) )
    !   if( present(y)) then
    !     ddj=1.d0/dfdj(y)
    !   else
    !     stop ' !! y needed in function grad_y !!'
    !   endif
    !   !
    ! else
    !   print*,' ** Use the saved grid jacobian matrix'
    ! endif
    !
    df=dfdj(f)*ddj
    !
    deallocate(afj,ccj)
    !
    print*,' ** Spatial gradient along y calculated'
    !
  end function grad_y
  !+-------------------------------------------------------------------+
  ! end of the function grad_y.                                        |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !|this subroutine calculates grid jacobian matrix.                   |
  !+-------------------------------------------------------------------+
  !|Writen by Jian Fang , 2020-04-18.                                  |
  !+-------------------------------------------------------------------+
  subroutine gridjacobian_xy(x,y,ddi,ddj,wallnormal)
    !
    use commvardefine, only: im,jm,lihomo,ljhomo
    !
    ! arguments
    real(8),intent(in) :: x(0:im,0:jm),y(0:im,0:jm)
    real(8),intent(out),optional :: ddi(2,0:im,0:jm),ddj(2,0:im,0:jm)
    real(8),intent(out),optional :: wallnormal(2,0:im)
    !
    ! local data
    integer :: i,j
    real(8) :: var1
    real(8),allocatable :: dx(:,:,:),dy(:,:,:)
    real(8),allocatable :: di(:,:,:),dj(:,:,:)
    !
    if(.not.allocated(afi)) then
      allocate(afi(0:im),cci(1:2,0:im))
      call tds_ini(afi,cci,im,lihomo)
    endif
    if(.not.allocated(afj)) then
      allocate(afj(0:jm),ccj(1:2,0:jm))
      call tds_ini(afj,ccj,jm,ljhomo)
    endif
    !
    allocate(dx(2,0:im,0:jm),dy(2,0:im,0:jm),                          &
             di(2,0:im,0:jm),dj(2,0:im,0:jm))
    !
    do j=0,jm
      dx(1,:,j)=dfdi(x(:,j))
      dy(1,:,j)=dfdi(y(:,j))
    end do
    !
    do i=0,im
      dx(2,i,:)=dfdj(x(i,:))
      dy(2,i,:)=dfdj(y(i,:))
    end do
    !
    do j=0,jm
    do i=0,im
      !
      var1=dx(1,i,j)*dy(2,i,j)-dx(2,i,j)*dy(1,i,j)
      !
      var1=1.d0/var1
      !
      di(1,i,j)= var1*dy(2,i,j)
      di(2,i,j)=-var1*dx(2,i,j)
      !
      dj(1,i,j)=-var1*dy(1,i,j)
      dj(2,i,j)= var1*dx(1,i,j)
      !
    end do
    end do
    !
    if(present(ddi) .and. present(ddj)) then
      ddi=di
      ddj=dj
    endif
    !
    if(present(wallnormal)) then
      do i=0,im
        !
        var1=dsqrt(dj(1,i,0)**2+dj(2,i,0)**2)
        !
        wallnormal(1,i)=dj(1,i,0)/var1
        wallnormal(2,i)=dj(2,i,0)/var1
        !
      end do
    endif
    !
    deallocate(dx,dy)
    !
    print*, ' ** Grid Jacobian matrix in a x-y plane is calculated'
    !
  end subroutine gridjacobian_xy
  !
  subroutine gridjacobian_3d(x,y,z,ddi,ddj,ddk)
    !
    use commvardefine, only: im,jm,km
    !
    ! arguments
    real(8),intent(in) :: x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),         &
                          z(0:im,0:jm,0:km)
    real(8),intent(out) :: ddi(1:3,0:im,0:jm,0:km),                    &
                           ddj(1:3,0:im,0:jm,0:km),                    &
                           ddk(1:3,0:im,0:jm,0:km)
    !
    ! local data
    integer :: i,j,k
    real(8) :: var1
    real(8),allocatable :: dx(:,:,:,:),dy(:,:,:,:),dz(:,:,:,:)
    !
    allocate(dx(1:3,0:im,0:jm,0:km),dy(1:3,0:im,0:jm,0:km),            &
             dz(1:3,0:im,0:jm,0:km))
    !
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,var1)
    !
    !$OMP DO
    do k=0,km
    do j=0,jm
      dx(1,:,j,k)=dfdi(x(:,j,k))
      dy(1,:,j,k)=dfdi(y(:,j,k))
      dz(1,:,j,k)=dfdi(z(:,j,k))
    end do
    end do
    !$OMP END DO
    !
    !$OMP DO
    do k=0,km
    do i=0,im
      dx(2,i,:,k)=dfdj(x(i,:,k))
      dy(2,i,:,k)=dfdj(y(i,:,k))
      dz(2,i,:,k)=dfdj(z(i,:,k))
    end do
    end do
    !$OMP END DO
    !
    !$OMP DO
    do j=0,jm
    do i=0,im
      dx(3,i,j,:)=dfdk(x(i,j,:))
      dy(3,i,j,:)=dfdk(y(i,j,:))
      dz(3,i,j,:)=dfdk(z(i,j,:))
    end do
    end do
    !$OMP END DO
    !
    !$OMP DO
    do k=0,km
    do j=0,jm
    do i=0,im
      var1= dx(1,i,j,k)*dy(2,i,j,k)*dz(3,i,j,k)                        &
           +dx(2,i,j,k)*dy(3,i,j,k)*dz(1,i,j,k)                        &
           +dx(3,i,j,k)*dy(1,i,j,k)*dz(2,i,j,k)                        &
           -dx(3,i,j,k)*dy(2,i,j,k)*dz(1,i,j,k)                        &
           -dx(2,i,j,k)*dy(1,i,j,k)*dz(3,i,j,k)                        &
           -dx(1,i,j,k)*dy(3,i,j,k)*dz(2,i,j,k)
      !
      ddi(1,i,j,k)=dy(2,i,j,k)*dz(3,i,j,k)-dy(3,i,j,k)*dz(2,i,j,k)
      ddi(2,i,j,k)=dx(3,i,j,k)*dz(2,i,j,k)-dx(2,i,j,k)*dz(3,i,j,k)
      ddi(3,i,j,k)=dx(2,i,j,k)*dy(3,i,j,k)-dx(3,i,j,k)*dy(2,i,j,k)
      !
      ddj(1,i,j,k)=dy(3,i,j,k)*dz(1,i,j,k)-dy(1,i,j,k)*dz(3,i,j,k)
      ddj(2,i,j,k)=dx(1,i,j,k)*dz(3,i,j,k)-dx(3,i,j,k)*dz(1,i,j,k)
      ddj(3,i,j,k)=dx(3,i,j,k)*dy(1,i,j,k)-dx(1,i,j,k)*dy(3,i,j,k)
      !
      ddk(1,i,j,k)=dy(1,i,j,k)*dz(2,i,j,k)-dy(2,i,j,k)*dz(1,i,j,k)
      ddk(2,i,j,k)=dx(2,i,j,k)*dz(1,i,j,k)-dx(1,i,j,k)*dz(2,i,j,k)
      ddk(3,i,j,k)=dx(1,i,j,k)*dy(2,i,j,k)-dx(2,i,j,k)*dy(1,i,j,k)
      !
      ddi(1,i,j,k)=ddi(1,i,j,k)/var1
      ddi(2,i,j,k)=ddi(2,i,j,k)/var1
      ddi(3,i,j,k)=ddi(3,i,j,k)/var1
      ddj(1,i,j,k)=ddj(1,i,j,k)/var1
      ddj(2,i,j,k)=ddj(2,i,j,k)/var1
      ddj(3,i,j,k)=ddj(3,i,j,k)/var1
      ddk(1,i,j,k)=ddk(1,i,j,k)/var1
      ddk(2,i,j,k)=ddk(2,i,j,k)/var1
      ddk(3,i,j,k)=ddk(3,i,j,k)/var1
    end do
    end do
    end do
    !$OMP END DO
    !
    !$OMP END PARALLEL
    !
    deallocate(dx,dy,dz)
    !
    print*, ' ** Grid Jacobian matrix is calculated'
    !
  end subroutine gridjacobian_3d
  !+-------------------------------------------------------------------+
  ! end of the function gridjacobian.                                  |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !|these functions are used to get spatial derivatives                |
  !+-------------------------------------------------------------------+
  function dfdi(var)
    !
    use commvardefine, only: im,lihomo
    !
    real(8),allocatable :: dfdi(:)
    real(8),intent(in) :: var(0:im)
    !
    real(8),allocatable :: vtemp(:)
    !
    allocate(dfdi(0:im))
    !
    vtemp=tds_rhs(var,afi,642,geom=lihomo)
    dfdi=tds_cal(vtemp,cci,afi)
    !
    deallocate(vtemp)
    !
  end function dfdi
  !
  function dfdj(var)
    !
    use commvardefine, only: jm,ljhomo
    !
    real(8),allocatable :: dfdj(:)
    real(8),intent(in) :: var(0:jm)
    !
    real(8) ,allocatable :: vtemp(:)
    !
    allocate(dfdj(0:jm))
    !
    vtemp=tds_rhs(var,afj,642,geom=ljhomo)
    dfdj =tds_cal(vtemp,ccj,afj)
    !
    deallocate(vtemp)
    !
  end function dfdj
  !
  function dfdk(var)
    !
    use commvardefine, only: km,lkhomo
    !
    real(8) ,allocatable :: dfdk(:)
    real(8),intent(in) :: var(0:km)
    !
    real(8) ,allocatable :: vtemp(:)
    !
    if(.not.allocated(afk)) then
      allocate(afk(0:km),cck(1:2,0:km))
      call tds_ini(afk,cck,km,lkhomo)
    endif
    !
    allocate(dfdk(0:km))
    !
    vtemp=tds_rhs(var,afk,642,geom=lkhomo)
    dfdk =tds_cal(vtemp,cck,afk)
    !
    deallocate(vtemp)
    !
  end function dfdk
  !+-------------------------------------------------------------------+
  !|end of the subroutine dfd*.                                        |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  ! this subroutine is used to initial for solving the tridiagonal 
  ! matrix with two layer of boundary scheme:
  ! a*x=b
  !   |1,af1,.............|
  !   |af2,1,af2,.... ....|
  !   |..af,1,af,.........|
  !   |...................|
  ! a=|...,af,1,af,.......|
  !   |...................|
  !   |.........,af,1,af..|
  !   |.........,af2,1,af2|
  !   |.............,af1,1|
  !
  !+-------------------------------------------------------------------+
  !|writen by fang jian, 2008-11-04.                                   |
  !+-------------------------------------------------------------------+
  subroutine tds_ini(af,cc,dim,homo)
    !
    ! arguments
    integer,intent(in) ::  dim 
    logical,intent(in) :: homo
    !
    real(8),intent(out) :: af(0:dim) 
    real(8),intent(out) :: cc(1:2,0:dim) 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af1,af2,af: input dat
    ! cc: output array
    ! dim: input dat
    ! l: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: l
    !
    if(homo) then
      af(0)=0.d0
      af(1)=1.d0/3.d0
      af(dim-1)=1.d0/3.d0
      af(dim)=0.d0
    else
      af(0)=1.d0
      af(1)=0.25d0
      af(dim-1)=0.25d0
      af(dim)=1.d0
    endif
    !
    af(2:dim-2)=1.d0/3.d0
    !
    l=0
    cc(1,l)=af(l)
    do l=1,dim-1
      cc(1,l)=af(l)/(1.d0-af(l)*cc(1,l-1))
    enddo
    !
  end subroutine tds_ini
  !+-------------------------------------------------------------------+
  ! end of the subroutine tds_ini.                                     |
  !+-------------------------------------------------------------------+
  !+-------------------------------------------------------------------+
  ! this subroutine is used to calculate the rhs of tridiagonal matrix |
  ! for a*x=b.                                                         |
  !+-------------------------------------------------------------------+
  function tds_rhs(vin,af,ns,geom) result(vout)
    !
    real(8),allocatable ::  vout(:)
    !
    real(8),intent(in) ::  vin(0:),af(0:)
    integer,intent(in) :: ns
    logical,intent(in),optional :: geom
    !
    integer :: dim,l,lp1,lp2,lp3,lm1,lm2,lm3
    logical :: ge
    real(8),allocatable ::  vos(:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ns=601: 6-o compact central scheme with
    ! boundary scheme: 2-4-6.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if (present(geom)) then
       ge = geom
    else
       ge = .false.
    end if
    !
    dim=size(af)-1
    !
    allocate(vout(0:dim),vos(-3:dim+3))
    !
    vos(0:dim)=vin(0:dim)
    if(ge) then
      do l=1,3
        vos(0-l)=vin(dim-l)-vin(dim)+vin(0)
        vos(dim+l)=vin(l)-vin(0)+vin(dim)
      enddo
    else
      do l=1,3
        vos(0-l)=vin(dim-l)
        vos(dim+l)=vin(l)
      enddo
    endif
    !
    do l=0,dim
      !
      lp1=l+1
      lp2=l+2
      lp3=l+3
      lm1=l-1
      lm2=l-2
      lm3=l-3
      !
      if(abs(af(l))<1.d-16) then
        ! interface point
        vout(l)=0.75d0*(vos(lp1)-vos(lm1))-                            &
                0.15d0*(vos(lp2)-vos(lm2))+                            &
                1.666666666666667d-2*(vos(lp3)-vos(lm3))
      elseif(abs(af(l)-1.d0/3.d0)<1.d-16) then

        vout(l)=0.777777777777778d0*(vos(lp1)-vos(lm1))+               &
               2.777777777777778d-2*(vos(lp2)-vos(lm2))
      elseif(abs(af(l)-0.25d0)<1.d-16) then
        vout(l)=0.75d0*(vos(lp1)-vos(lm1))
      elseif(abs(af(l)-1.d0)<1.d-16) then
        if(l==0) then
          vout(0)=  -2.d0*    vos(0)+  2.d0*  vos(1)
        elseif(l==dim) then
          vout(dim)=  -2.d0*vos(dim-1)+  2.d0*vos(dim)
        else
          stop 'error 1 in tds_rhs'
        endif
        !
      else
        print*,'af=',abs(af(l))
        stop 'non-defined scheme in tds_rhs'
      endif
      !
    enddo
    !
    deallocate(vos)
    !
  end function tds_rhs
  !+-------------------------------------------------------------------+
  ! end of the  function tds_rhs                                       |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  ! this subroutine is used to solve the tridiagonal matrix in the 
  ! programme with two layer of boundary scheme:
  ! a*x=b
  !   |1,af1,.............|
  !   |af2,1,af2,.... ....|
  !   |..af,1,af,.........|
  !   |...................|
  ! a=|...,af,1,af,.......|
  !   |...................|
  !   |.........,af,1,af..|
  !   |.........,af2,1,af2|
  !   |.............,af1,1|
  !
  !+-------------------------------------------------------------------+
  ! writen by fang jian, 2008-11-04.                                   |
  !+-------------------------------------------------------------------+
  function tds_cal(bd,cc,af) result(xd)
    !
    real(8),allocatable :: xd(:)
    real(8),intent(in) :: af(0:),cc(:,0:),bd(0:)
    !
    integer :: l,dim
    real(8),allocatable,dimension(:) :: yd
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af: input dat
    ! bd: input array
    ! xd: output array
    ! cc: input array
    ! dim: input dat
    ! l, yd: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !af1=1.d0
    !af2=0.25d0
    !af=0.333333333333333d0
    !
    dim=size(af)-1
    allocate( xd(0:dim),yd(0:dim) )
    !
    yd(0)=bd(0)
    do l=1,dim
      yd(l)=(bd(l)-af(l)*yd(l-1))/(1.d0-af(l)*cc(1,l-1))
    enddo
    !
    xd(dim)=yd(dim)
    do l=dim-1,0,-1
      xd(l)=yd(l)-cc(1,l)*xd(l+1)
    end do
    !
    deallocate( yd )
    !
  end function tds_cal
  !+-------------------------------------------------------------------+
  !|end of the subroutine tds_cal.                                     |
  !+-------------------------------------------------------------------+
  !!
  !
end module gradsolver
!+---------------------------------------------------------------------+
!| The end of the module gradsolver.                                   |
!+---------------------------------------------------------------------+