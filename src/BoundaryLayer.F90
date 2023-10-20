!+---------------------------------------------------------------------+
!| This module contains subroutines to process BL data.                |
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2020-04-19                                     |
!+---------------------------------------------------------------------+
module boundarylayer
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to get statistics for channel flow.            |
  !+-------------------------------------------------------------------+
  subroutine chanmean
    !
    use commvardefine,only: im,jm,km,Reynolds,Mach,pinf,gridfile
    use h5readwrite
    use basicfunction
    use interpolation 
    use writetec
    !
    ! local data
    real(8),allocatable,dimension(:) :: x,y,z,ro,u1,p,t,dvdx,dudy,marms,malo
    real(8),allocatable,dimension(:) :: yplus,uplus
    real(8),allocatable,dimension(:) :: u11,u22,u33,u12,tke,miut
    real(8),allocatable,dimension(:) :: convect,pre_ace,pre_dil,     &
                                        pre_tra,product,tur_tra,     &
                                        vis_ace,vis_dif,dissipa,balance
    real(8),allocatable,dimension(:,:) :: convect_xy,pre_ace_xy,pre_dil_xy,     &
                                          pre_tra_xy,product_xy,tur_tra_xy,     &
                                          vis_ace_xy,vis_dif_xy,dissipa_xy,     &
                                          balance_xy,x2d,y2d
    !
    real(8) :: utaw,taw,var1,lvis,miu
    integer :: j,jm2
    !
    allocate( y(0:jm),ro(0:jm),u1(0:jm),p(0:jm),t(0:jm),dvdx(0:jm),    &
              dudy(0:jm),miut(0:jm),x(0:im),z(0:km) )
    !
    call H5ReadSubset(x,im,jm,km,'x',trim(gridfile),jslice=1,kslice=1)
    call H5ReadSubset(z,im,jm,km,'z',trim(gridfile),islice=1,jslice=1)
    call H5ReadSubset(y,im,jm,km,'y',trim(gridfile),islice=1,kslice=1)
    !
    call H5ReadArray(ro,  jm,   'ro','Results/mean.fav.xzm.h5')
    call H5ReadArray(u1,  jm,   'u1','Results/mean.fav.xzm.h5')
    call H5ReadArray(p,   jm,    'p','Results/mean.fav.xzm.h5')
    call H5ReadArray(t,   jm,    't','Results/mean.fav.xzm.h5')
    call H5ReadArray(dudy,jm,'du1dy','Results/mean.fav.xzm.h5')
    call H5ReadArray(dvdx,jm,'du2dx','Results/mean.fav.xzm.h5')
    !
    allocate(u11(0:jm),u22(0:jm),u33(0:jm),u12(0:jm),tke(0:jm))
    !
    call H5ReadArray(u11,jm,'u11','Results/2order.fav.xzm.h5')
    call H5ReadArray(u22,jm,'u22','Results/2order.fav.xzm.h5')
    call H5ReadArray(u33,jm,'u33','Results/2order.fav.xzm.h5')
    call H5ReadArray(u12,jm,'u12','Results/2order.fav.xzm.h5')
    !
    jm2=jm/2
    !
    print*,Reynolds
    !
    ! taw=MiuCal(t(0))/Reynolds*dudy(0)
    call upluscal(uplus,yplus,u1,y,ro,t(0),utaw=utaw)
    !
    print*,' ** utau=',utaw
    ! 
    lvis=ro(0)*utaw/miucal(T(0))*Reynolds
    open(18,file='Results/mesh_report.dat')
    write(18,"(A)")'mesh resolution'
    write(18,"(A)")'----------------------------------------------------'
    write(18,"(A)")'          Δx+         Δy1+         Δym+          Δz+'
    write(18,"(A)")'----------------------------------------------------'
    write(18,"(4(1X,F12.7))")(x(1)-x(0))*lvis,(y(1)-y(0))*lvis,           &
                             (y(jm2)-y(jm2-1))*lvis,(z(1)-z(0))*lvis
    write(18,"(A)")'----------------------------------------------------'
    close(18)
    print*,' << mesh_report.dat ... done !'
    !
    open(18,file='Results/uplus.dat')
    write(18,"(2(1X,A15))")'yplus','uplus'
    write(18,"(2(1X,E15.7E3))")(yplus(j),uplus(j),j=0,jm2)
    close(18)
    print*,' << Results/uplus.dat ... done. '
    !
    tke=0.5d0*(u11+u22+u33)
    !
    open(18,file='Results/Restress.dat')
    write(18,"(6(1X,A15))")'yplus','uu','vv','ww','tke','u12'
    do j=0,jm2
      var1=ro(j)/ro(0)/utaw/utaw
      write(18,"(6(1X,E15.7E3))")yplus(j),u11(j)*var1,u22(j)*var1,     &
                                     u33(j)*var1,tke(j)*var1,u12(j)*var1
    end do
    close(18)
    print*,' << Results/Restress.dat'
    !
    allocate(marms(0:jm),malo(0:jm))
    do j=0,jm
      var1=sqrt(abs(u11(j)+u22(j)+u33(j)))
      marms(j)=var1/sqrt(t(j))*Mach
      malo(j)=u1(j)/sqrt(t(j))*Mach
      ! var1=max(abs(dudy(j)),0.01d0)
      ! var1=sign(1.d0,dudy(j))*var1
      var1=dudy(j)
      miut(j)=-ro(j)*u12(j)/var1*Reynolds
    enddo
    !
    open(18,file='Results/profile.dat')
    write(18,"(7(1X,A15))")'y','u','T','p','rho','Mrms','M'
    do j=0,jm
      write(18,"(7(1X,E15.7E3))")y(j),u1(j),t(j),p(j),ro(j),           &
                                                        marms(j),malo(j)
    end do
    close(18)
    print*,' << Results/profile.dat'
    !
    open(18,file='Results/miut.dat')
    write(18,"(4(1X,A15))")'y','miut','uv','dudy'
    do j=0,jm
      write(18,"(4(1X,E15.7E3))")y(j),miut(j),u12(j),dudy(j)
    end do
    close(18)
    print*,' << Results/miut.dat'
    !
    ! budget terms
    allocate(convect(0:jm),pre_ace(0:jm),pre_dil(0:jm),         &
             pre_tra(0:jm),product(0:jm),tur_tra(0:jm),         &
             vis_ace(0:jm),vis_dif(0:jm),dissipa(0:jm),         &
             balance(0:jm) )
    !
    call H5ReadArray(convect,jm,'convection','Results/tke_budget.xzm.h5')
    call H5ReadArray(dissipa,jm,'dissipation','Results/tke_budget.xzm.h5')
    call H5ReadArray(pre_ace,jm,'pressure_accelaration','Results/tke_budget.xzm.h5')
    call H5ReadArray(pre_dil,jm,'pressure_dilatation','Results/tke_budget.xzm.h5')
    call H5ReadArray(pre_tra,jm,'pressure_transport','Results/tke_budget.xzm.h5')
    call H5ReadArray(product,jm,'production','Results/tke_budget.xzm.h5')
    call H5ReadArray(tur_tra,jm,'turbulent_transport','Results/tke_budget.xzm.h5')
    call H5ReadArray(vis_ace,jm,'viscous_accelaration','Results/tke_budget.xzm.h5')
    call H5ReadArray(vis_dif,jm,'viscous_diffusion','Results/tke_budget.xzm.h5')
    !
    balance=convect+pre_ace+pre_dil+pre_tra+product+tur_tra+vis_ace+vis_dif+dissipa
    !
    open(18,file='Results/budget_profile.dat')
    write(18,"(12(1X,A15))")'y','yplus','convection',                  &
            'pressure_accelaration','pressure_dilatation',             &
            'pressure_transport','production','turbulent_transport',   &
            'viscous_accelaration','viscous_diffusion','dissipation',  &
            'balance'
    do j=0,jm
      miu=miucal(t(0))/Reynolds
      var1=ro(0)**2*utaw**4/miu
      write(18,"(12(1X,E15.7E3))")y(j),yplus(j),                &
                convect(j)/var1,pre_ace(j)/var1,pre_dil(j)/var1, &
                pre_tra(j)/var1,product(j)/var1,tur_tra(j)/var1, &
                vis_ace(j)/var1,vis_dif(j)/var1,dissipa(j)/var1, &
                balance(j)/var1
    end do
    close(18)
    print*,' << budget_profile.dat ... done. '
    !
    call writebudgetlay('budget_profile.dat','Results/plot_budget.lay')
    !
    allocate(convect_xy(0:im,0:jm),pre_ace_xy(0:im,0:jm),pre_dil_xy(0:im,0:jm),     &
             pre_tra_xy(0:im,0:jm),product_xy(0:im,0:jm),tur_tra_xy(0:im,0:jm),     &
             vis_ace_xy(0:im,0:jm),vis_dif_xy(0:im,0:jm),dissipa_xy(0:im,0:jm),     &
             balance_xy(0:im,0:jm),x2d(0:im,0:jm),y2d(0:im,0:jm))
    !
    call H5ReadArray(convect_xy,im,jm,'convection','Results/tke_budget.zm.h5')
    call H5ReadArray(dissipa_xy,im,jm,'dissipation','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_ace_xy,im,jm,'pressure_accelaration','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_dil_xy,im,jm,'pressure_dilatation','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_tra_xy,im,jm,'pressure_transport','Results/tke_budget.zm.h5')
    call H5ReadArray(product_xy,im,jm,'production','Results/tke_budget.zm.h5')
    call H5ReadArray(tur_tra_xy,im,jm,'turbulent_transport','Results/tke_budget.zm.h5')
    call H5ReadArray(vis_ace_xy,im,jm,'viscous_accelaration','Results/tke_budget.zm.h5')
    call H5ReadArray(vis_dif_xy,im,jm,'viscous_diffusion','Results/tke_budget.zm.h5')
    !
    miu=miucal(t(0))/Reynolds
    var1=ro(0)**2*utaw**4/miu
    !
    convect_xy=convect_xy/var1
    dissipa_xy=dissipa_xy/var1
    pre_ace_xy=pre_ace_xy/var1
    pre_dil_xy=pre_dil_xy/var1
    pre_tra_xy=pre_tra_xy/var1
    product_xy=product_xy/var1
    tur_tra_xy=tur_tra_xy/var1
    vis_ace_xy=vis_ace_xy/var1
    vis_dif_xy=vis_dif_xy/var1
    !
    balance_xy=convect_xy+pre_ace_xy+pre_dil_xy+pre_tra_xy + &
               product_xy+tur_tra_xy+vis_ace_xy+vis_dif_xy+dissipa_xy
    !
    call H5ReadSubset(x2d,im,jm,km,'x',trim(gridfile),kslice=0)
    call H5ReadSubset(y2d,im,jm,km,'y',trim(gridfile),kslice=0)
    !
    call writetecbin('Results/tecbudget.plt', x2d,'x',y2d,'y', &
                                      convect_xy,'convection', &
                                      pre_ace_xy,'pressure_accelaration', &
                                      pre_dil_xy,'pressure_dilatation',   &
                                      pre_tra_xy,'pressure_transport',    &
                                      product_xy,'production',            &
                                      tur_tra_xy,'turbulent_transport',   &
                                      vis_ace_xy,'viscous_accelaration',  &
                                      vis_dif_xy,'viscous_diffusion',     &
                                      dissipa_xy,'dissipation',           &
                                      balance_xy,'balance',im,jm)
    !
  end subroutine chanmean
  !+-------------------------------------------------------------------+
  !| The end of the subroutine chanmean.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to get statistics for boundary layer.          |
  !+-------------------------------------------------------------------+
  subroutine blmean(iref)
    !
    use commvardefine,only: im,jm,km,Reynolds,Prandtl,Mach,gamma,pinf,gridfile
    use h5readwrite
    use basicfunction
    use interpolation 
    use writetec
    !
    ! argument
    integer,intent(in) :: iref
    !
    ! local data
    integer :: i,j,je99,jeir
    ! jeir: edge of the rotational part of the boundary layer (δe) is 
    ! defined as the point where the mean spanwise vorticity becomes 
    ! less than a suitable threshold value (here set to 0.005 u∞/δin).
    !
    real(8),allocatable,dimension(:,:) :: x,y,ro,u1,u2,p,t,dudy,dvdx,omegaz,dtdy
    real(8),allocatable,dimension(:,:) :: u11,u22,u33,u12,tke,tt,pp
    real(8),allocatable,dimension(:) :: z,yplus,uplus,bl_delta99,      &
                                        bl_dstar,bl_theta,cf,shapefac, &
                                        Retau,Re_thetai,cfi,ch
    real(8),allocatable,dimension(:,:) :: convect,pre_ace,pre_dil,     &
                                          pre_tra,product,tur_tra,     &
                                          vis_ace,vis_dif,dissipa,balance
    real(8) :: utaw,lvis,omegaz_ref,miu,var1,var2,lref
    character(len=5) :: irefname
    !
    print*,' ** get statistics at the location i=',iref
    if(iref<=0) stop
    write(irefname,'(1hi,I4.4)')iref
    !
    allocate(x(0:im,0:jm),y(0:im,0:jm),z(0:km))
    call H5ReadSubset(x,im,jm,km,'x',trim(gridfile),kslice=0)
    call H5ReadSubset(y,im,jm,km,'y',trim(gridfile),kslice=0)
    call H5ReadSubset(z,im,jm,km,'z',trim(gridfile),islice=iref,jslice=0)
    !
    ! open(18,file='yn.dat')
    ! do j=0,jm
    ! write(18,*)j,y(i,j)
    ! enddo
    ! close(18)
    ! print*,' << yn.dat'
    !
    allocate(ro(0:im,0:jm),u1(0:im,0:jm),p(0:im,0:jm),                 &
             t(0:im,0:jm),dudy(0:im,0:jm),dvdx(0:im,0:jm),             &
             omegaz(0:im,0:jm),u2(0:im,0:jm),dtdy(0:im,0:jm))
    !
    call H5ReadArray(ro,im,jm,'ro','Results/mean.fav.zm.h5')
    call H5ReadArray(u1,im,jm,'u1','Results/mean.fav.zm.h5')
    call H5ReadArray(u2,im,jm,'u2','Results/mean.fav.zm.h5')
    call H5ReadArray(p, im,jm, 'p','Results/mean.fav.zm.h5')
    call H5ReadArray(t, im,jm, 't','Results/mean.fav.zm.h5')
    call H5ReadArray(dudy,im,jm,'du1dy','Results/mean.fav.zm.h5')
    call H5ReadArray(dvdx,im,jm,'du2dx','Results/mean.fav.zm.h5')
    call H5ReadArray(dtdy,im,jm,'du1dy','Results/mean.fav.zm.h5')
    !
    allocate(cf(0:im),Retau(0:im),cfi(0:im),ch(0:im))
    do i=0,im
      miu=miucal(t(i,0))/Reynolds
      cf(i)=miu*dudy(i,0)*2.d0
      !
      utaw=sqrt(miu*dudy(i,0)/ro(i,0))
      !
      Retau(i)=ro(i,0)*utaw/miu
      !
      var1=(t(i,0)-1.d0)/sqrt(t(i,0)*(t(i,0)-1.d0))
      var2=(t(i,0)-1.d0)/asin(var1)/asin(var1)
      !
      cfi(i)=var2*cf(i)
      !
      ch(i)=2.d0*miu*dtdy(i,0)/(Mach**2*Prandtl*(gamma-1))
    enddo
    !
    allocate(bl_delta99(0:im),bl_dstar(0:im),bl_theta(0:im),shapefac(0:im),Re_thetai(0:im))
    !
    do i=0,im
      bl_delta99(i)=0.d0
      do j=1,jm
        if(u1(i,j-1)<=0.99d0 .and. u1(i,j)>=0.99d0) then
          bl_delta99(i)=linear1d(u1(i,j-1),u1(i,j),y(i,j-1),y(i,j),0.99d0)
          exit
        endif
      enddo
      if(i==iref) je99=j
    enddo
    omegaz_ref=0.005d0/bl_delta99(iref)
    !
    omegaz=dudy-dvdx
    !
    do i=0,im
      !
      do j=1,jm
        if(omegaz(i,j-1)>=omegaz_ref .and. omegaz(i,j)<=omegaz_ref) then
          jeir=j
          exit
        endif
      enddo
      !
      bl_dstar(i)=0.d0
      bl_theta(i)=0.d0
      do j=1,jeir
        var1=0.5d0*(u1(i,j)*ro(i,j)+u1(i,j-1)*ro(i,j-1))
        bl_dstar(i)=bl_dstar(i)+(1.d0-var1)*(y(i,j)-y(i,j-1))
        !
        var1=0.5d0*(u1(i,j)*ro(i,j)+u1(i,j-1)*ro(i,j-1))
        var2=1.d0-0.5d0*(u1(i,j)+u1(i,j-1))
        bl_theta(i)=bl_theta(i)+var1*var2*(y(i,j)-y(i,j-1))
      enddo
      !
      shapefac(i)=bl_dstar(i)/bl_theta(i)
      !
      Re_thetai(i)=Reynolds*bl_theta(i)/MiuCal(t(i,0))
      !
    enddo
    !
    write(*,"(A,F7.3,A)")'  ---------- BL thickness at x=',x(iref,0),' ----------'
    write(*,"(A)")'            δ          δ*           θ           H'
    write(*,"(1X,4(F12.7))")bl_delta99(iref),bl_dstar(iref),           &
                                           bl_theta(iref),shapefac(iref)
    write(*,"(A)")'  -----------------------------------------------'
    !
    open(18,file='Results/BLParameter.dat')
    write(18,"(13(1X,A15))")'x','Cf','delta','delta*','theta','H','Re_delta', &
                                            'Re_delta*','Re_theta','Re_tau','cfi','Re_thetai','Ch'
    write(18,"(13(1X,E15.7E3))")(x(i,0),cf(i),bl_delta99(i),            &
                               bl_dstar(i),bl_theta(i),shapefac(i),     &
                        Reynolds*bl_delta99(i),Reynolds*bl_dstar(i),    &
                        Reynolds*bl_theta(i),Retau(i)*bl_delta99(i),    &
                        cfi(i),Re_thetai(i),ch(i),i=0,im)
    close(18)
    print*,' << BLParameter.dat ... done !'
    !
    !+-----------------------------------+
    !| calculate profile at iref station.|
    !+-----------------------------------+
    lref=bl_delta99(iref)
    i=iref
    !
    write(*,'(A,I4)')'  ** Now processing profile at i=',i
    !
    allocate(yplus(0:jm),uplus(0:jm))
    !
    call upluscal(uplus,yplus,u1(i,:),y(i,:),ro(i,:),t(i,0),utaw=utaw)
    !
    print*,' ** utau=',utaw
    !
    open(18,file='Results/profile.'//irefname//'.dat')
    write(18,"(7(1X,A15))")'y','y/delta','yplus','ro','u','T','du/dy'
    write(18,"(7(1X,E15.7E3))")(y(i,j),y(i,j)/lref,yplus(j),ro(i,j),   &
                                        u1(i,j),t(i,j),dudy(i,j),j=0,jm)
    close(18)
    print*,' << profile.',irefname,'.dat ... done. '
    !
    open(18,file='inlet.dat')
    write(18,"(5(1X,A15))")'y','ro','u','v','T'
    write(18,"(5(1X,E15.7E3))")(y(i,j),ro(i,j),u1(i,j),u2(i,j),t(i,j),j=0,jm)
    close(18)
    print*,' << inlet.dat ... done. '
    !
    open(18,file='Results/uplus.'//irefname//'.dat')
    write(18,"(2(1X,A15))")'yplus','uplus'
    write(18,"(2(1X,E15.7E3))")(yplus(j),uplus(j),j=0,jm)
    close(18)
    print*,' << uplus.',irefname,'.dat ... done. '
    !
    call writeupluslay('uplus.'//irefname//'.dat',                    &
                       'Results/uplus-yplus.lay')
    !
    lvis=ro(iref,0)*utaw/miucal(T(iref,0))*Reynolds
    open(18,file='Results/mesh_report.'//irefname//'.dat')
    write(18,"(A,F12.7,A,I4)")'mesh resolution at x=',x(i,0),', i=',i
    write(18,"(A)")'-----------------------------------------------------------------'
    write(18,"(A)")'          Δx+         Δy1+         Δye+          Δz+          Lz+'
    write(18,"(A)")'-----------------------------------------------------------------'
    write(18,"(5(1X,F12.7))")0.5d0*(x(i+1,0)-x(i-1,0))*lvis,           &
                                       (y(i,1)-y(i,0))*lvis,           &
                               (y(i,je99)-y(i,je99-1))*lvis,           &
                                           (z(1)-z(0))*lvis,           &
                                          (z(km)-z(0))*lvis
    write(18,"(A)")'-----------------------------------------------------------------'
    close(18)
    print*,' << mesh_report.',irefname,'.dat ... done !'
    !
    allocate(u11(0:im,0:jm),u22(0:im,0:jm),u33(0:im,0:jm),             &
             u12(0:im,0:jm),tke(0:im,0:jm),tt(0:im,0:jm),pp(0:im,0:jm))
    !
    call H5ReadArray(u11,im,jm,'u11','Results/2order.fav.zm.h5')
    call H5ReadArray(u22,im,jm,'u22','Results/2order.fav.zm.h5')
    call H5ReadArray(u33,im,jm,'u33','Results/2order.fav.zm.h5')
    call H5ReadArray(u12,im,jm,'u12','Results/2order.fav.zm.h5')
    call H5ReadArray( tt,im,jm, 'tt','Results/2order.fav.zm.h5')
    call H5ReadArray( pp,im,jm, 'pp','Results/2order.fav.zm.h5')
    !
    tke=0.5d0*(u11+u22+u33)
    !
    open(18,file='Results/Restress_profile.'//irefname//'.dat')
    write(18,"(9(1X,A15))")'y/delta','yplus','uu','vv','ww','tke','u12','tt','pp'
    do j=0,jm
      var1=ro(i,j)/ro(i,0)/utaw/utaw
      write(18,"(9(1X,E15.7E3))")y(i,j)/lref,yplus(j),                 &
             u11(i,j)*var1,u22(i,j)*var1,u33(i,j)*var1,tke(i,j)*var1,  &
             u12(i,j)*var1,tt(i,j),pp(i,j)/(pinf**2)
    end do
    close(18)
    print*,' << Restress_profile.',irefname,'.dat ... done. '
    !
    call writerestresslay('Restress_profile.'//irefname//'.dat',       &
                          'Results/Restress.lay')
    !
    ! budget terms
    allocate(convect(0:im,0:jm),pre_ace(0:im,0:jm),pre_dil(0:im,0:jm), &
             pre_tra(0:im,0:jm),product(0:im,0:jm),tur_tra(0:im,0:jm), &
             vis_ace(0:im,0:jm),vis_dif(0:im,0:jm),dissipa(0:im,0:jm), &
             balance(0:im,0:jm) )
    !
    call H5ReadArray(convect,im,jm,'convection','Results/tke_budget.zm.h5')
    call H5ReadArray(dissipa,im,jm,'dissipation','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_ace,im,jm,'pressure_accelaration','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_dil,im,jm,'pressure_dilatation','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_tra,im,jm,'pressure_transport','Results/tke_budget.zm.h5')
    call H5ReadArray(product,im,jm,'production','Results/tke_budget.zm.h5')
    call H5ReadArray(tur_tra,im,jm,'turbulent_transport','Results/tke_budget.zm.h5')
    call H5ReadArray(vis_ace,im,jm,'viscous_accelaration','Results/tke_budget.zm.h5')
    call H5ReadArray(vis_dif,im,jm,'viscous_diffusion','Results/tke_budget.zm.h5')
    !
    balance=convect+pre_ace+pre_dil+pre_tra+product+tur_tra+vis_ace+vis_dif+dissipa
    !
    i=iref
    open(18,file='Results/budget_profile.'//irefname//'.dat')
    write(18,"(12(1X,A15))")'y/delta','yplus','convection',            &
            'pressure_accelaration','pressure_dilatation',             &
            'pressure_transport','production','turbulent_transport',   &
            'viscous_accelaration','viscous_diffusion','dissipation',  &
            'balance'
    do j=0,jm
      miu=miucal(t(i,0))/Reynolds
      var1=ro(i,0)**2*utaw**4/miu
      write(18,"(12(1X,E15.7E3))")y(i,j)/lref,yplus(j),                &
                convect(i,j)/var1,pre_ace(i,j)/var1,pre_dil(i,j)/var1, &
                pre_tra(i,j)/var1,product(i,j)/var1,tur_tra(i,j)/var1, &
                vis_ace(i,j)/var1,vis_dif(i,j)/var1,dissipa(i,j)/var1, &
                balance(i,j)/var1
    end do
    close(18)
    print*,' << budget_profile.',irefname,'.dat ... done. '
    !
    call writebudgetlay('budget_profile.'//irefname//'.dat',           &
                        'Results/plot_budget.lay')
    !
    open(18,file='Results/mesh_report.'//irefname//'.dat',position='append')
    !
    write(18,"(A)")' ratio of mesh to kolmogorov scale at wall          '
    write(18,"(A)")'----------------------------------------------------'
    write(18,"(A)")'         Δx/η         Δy/η         Δz/η         Δe/η'
    write(18,"(A)")'----------------------------------------------------'
    !
    miu=miucal(t(i,0))/Reynolds
    var1=sqrt(sqrt(miu**3/ro(i,0)**2/(-dissipa(i,0))))
    write(18,"(4(1X,F12.7))")0.5d0*(x(i+1,0)-x(i-1,0))/var1,           &
                                       (y(i,1)-y(i,0))/var1,           &
                                           (z(1)-z(0))/var1,           &
                         cube_root(0.5d0*(x(i+1,0)-x(i-1,0))*          &
                                (y(i,1)-y(i,0))*(z(1)-z(0)))/var1
    write(18,"(A)")'----------------------------------------------------'
    !
    write(18,"(A)")' ratio of mesh to kolmogorov scale at BL edge       '
    write(18,"(A)")'----------------------------------------------------'
    write(18,"(A)")'         Δx/η         Δy/η         Δz/η         Δe/η'
    write(18,"(A)")'----------------------------------------------------'
    !
    miu=miucal(t(i,je99))/Reynolds
    var1=sqrt(sqrt(miu**3/ro(i,je99)**2/(-dissipa(i,je99))))
    write(18,"(4(1X,F12.7))")0.5d0*(x(i+1,je99)-x(i-1,je99))/var1,     &
                             0.5d0*(y(i,je99+1)-y(i,je99-1))/var1,     &
                                                 (z(1)-z(0))/var1,     &
                       cube_root(0.25d0*(x(i+1,je99)-x(i-1,je99))*     &
                           (y(i,je99+1)-y(i,je99-1))*(z(1)-z(0)))/var1
    write(18,"(A)")'----------------------------------------------------'
    close(18)
    print*,' << mesh_report.',irefname,'.dat ... done !'
    !
    !
  end subroutine blmean
  !+-------------------------------------------------------------------+
  !| The end of the subroutine blmean.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is used to calculate the thickness of a BL profile. |
  !+-------------------------------------------------------------------+
  function blthickness(y,u,uinf,ro,omegaz,jedge,btype)
    !
    use basicfunction
    use interpolation
    !
    real(8) :: blthickness
    real(8),intent(in) :: y(0:),u(0:)
    real(8),intent(in),optional :: ro(0:),omegaz(0:)
    integer,intent(out),optional :: jedge
    character(len=*),intent(in) :: btype
    !
    real(8),intent(in),optional :: uinf
    !
    ! local data
    integer :: dim,j,jeir
    real(8) :: vinf,roinf,normthick,omegaz_ref
    real(8) :: var1,var2
    !
    if(present(uinf)) then
      vinf=uinf
    else
      vinf=1.d0
    endif
    roinf=1.d0
    !
    dim=size(y)-1
    !
    do j=0,dim
      if(u(j-1)<=0.99d0*vinf .and. u(j)>=0.99d0*vinf) then
        normthick=linear1d(u(j-1),u(j),y(j-1),y(j),0.99d0*vinf)
        if(present(jedge)) jedge=j
        exit
      endif
    enddo
    !
    if(btype=='nominal') then
      !
      blthickness=normthick
      !
    elseif(btype=='displacement') then
      !
      if(present(omegaz)) then
        omegaz_ref=0.005d0*vinf/normthick
        do j=1,dim
          if(omegaz(j-1)>=omegaz_ref .and. omegaz(j)<=omegaz_ref) then
            jeir=j
            exit
          endif
        enddo
        print*,' .. edge of rotational boundar layer is at j=',jeir,'y=',y(jeir)
      else
        jeir=dim
      endif
      !
      blthickness=0.d0
      do j=1,jeir
        var1=0.5d0*(u(j)*ro(j)+u(j-1)*ro(j-1))/(roinf*vinf)
        blthickness=blthickness+(1.d0-var1)*(y(j)-y(j-1))
      enddo
      !
    elseif(btype=='momentum') then
      !
      if(present(omegaz)) then
        omegaz_ref=0.005d0*vinf/normthick
        do j=1,dim
          if(omegaz(j-1)>=omegaz_ref .and. omegaz(j)<=omegaz_ref) then
            jeir=j
            exit
          endif
        enddo
        print*,' .. edge of rotational boundar layer is at j=',jeir,'y=',y(jeir)
      else
        jeir=dim
      endif
      !
      blthickness=0.d0
      do j=1,jeir
        var1=0.5d0*(u(j)*ro(j)+u(j-1)*ro(j-1))/(roinf*vinf)
        var2=1.d0-0.5d0*(u(j)+u(j-1))/vinf
        blthickness=blthickness+var1*var2*(y(j)-y(j-1))
      enddo
      !
    else
      stop ' !! BL thickness not defined !!'
    endif
    !
  end function blthickness
  !+-------------------------------------------------------------------+
  !| The end of the function blthickness.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this function is used to give a turbulent boundary layer profiles |
  !+-------------------------------------------------------------------+
  !| ref: Rona, Aldo and Monti, Manuele and Airiau, Christophe,
  !| On the generation of the mean velocity profile for turbulent
  !| boundary layers with pressure gradient under equilibrium conditions
  !| (2012) The Aeronautical Journal, vol. 116 (n° 1180). ISSN 0001-9240
  !| eq.(6). 
  !+-------------------------------------------------------------------+
  function tblvelocity(y,utau,delta,rowall,twall,re) result(u)
    !
    use basicfunction, only: miucal
    !
    ! arguments
    real(8),allocatable :: u(:)
    real(8),intent(in) :: y(0:)
    real(8),intent(in) :: utau,delta,rowall,twall,re
    !
    ! local data
    integer :: j,m,je
    real(8) :: miu
    real(8),allocatable :: yplus(:),ydelta(:),ui(:),uplus(:)
    real(8) :: kamma,a,b,alpham,betam,eta,ppi,retau
    real(8) :: var1,var2,ueplus,uinfplus
    parameter(kamma=0.45d0,a=10.75d0,b=3.17d0)
    !
    m=size(y)-1
    allocate(u(0:m),yplus(0:m),ydelta(0:m),ui(0:m),uplus(0:m))
    u=0.d0
    !
    !
    print*,' ** boundary layer thickness=',delta
    !
    miu=miucal(twall)
    do j=0,m
      yplus(j)=rowall*utau*y(j)/miu*re
      ydelta(j)=y(j)/delta
    enddo
    retau=rowall*utau*delta/miu*re
    !
    write(*,"(2(A,F12.7))")'  ** y+_max= ',yplus(m),'; y_max=',ydelta(m)
    !
    ! Eq. (5) The profile by Musker(30) is a profile of velocity 
    ! gradient with an integral solution by partial fractions given in 
    ! Monkewitz et al(33) as,
    alpham=0.5d0*(a-1.d0/kamma)
    betam=sqrt(2.d0*a*alpham-alpham**2)
    !
    ! write(*,"(2(A,F10.7))")' ** αm=',alpham,'; βm=',betam
    !
    do j=0,m
      ui(j)=1.d0/kamma*log((yplus(j)+a)/a)
      !
      var1=a*((yplus(j)-alpham)**2+betam**2)/(2.d0*alpham*(yplus(j)+a)**2)
      var1=(a-4.d0*alpham)*log(var1)
      !
      var2=atan((yplus(j)-alpham)/betam)+atan(alpham/betam)
      var2=var2*2.d0*alpham*(5.d0*a-4.d0*alpham)/betam
      !
      ui(j)=ui(j)+(var1+var2)*alpham/(a+4.d0*alpham)
    enddo
    !
    do j=0,m
      ! eta=yplus(j)/retau
      eta=y(j)/delta*0.86d0
      !
      ueplus=0.99d0/utau
      uinfplus=1.d0/utau
      ppi=0.5d0*kamma*(ueplus-1.d0/kamma*log(retau)-b)
      !
      var1=2.d0*ppi/kamma*eta*eta*(3.d0-2.d0*eta)
      !
      if(eta<=1.d0) then
        uplus(j)=ui(j)+1.d0/kamma*eta*eta*(1.d0-eta)+var1
        je=j
      else
        var1=uplus(j-1)-uplus(j-2)
        var2=exp(-1.d0*(yplus(j)-yplus(je)))
        uplus(j)=uplus(j-1) !+var1*var2
      endif
      !
    enddo
    !
    do j=0,m
      u(j)=uplus(j)*utau
    enddo
    ! var1=u(m) !/0.99d0
    ! u=u/var1
    !
    ! print*,' ** dime=',m
    open(18,file='velocity_profile.dat')
    write(18,"(5(1X,A15))")'y','y/delta','yplus','uplus','u'
    write(18,"(5(1X,E15.7E3))")(y(j),ydelta(j),yplus(j),uplus(j),u(j),j=0,m)
    close(18)
    print*,' << velocity_profile.dat'
    !
  end function tblvelocity
  !+-------------------------------------------------------------------+
  !| The end of the function tblvelocity.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is used to generate Reynolds stress profiels.     |
  !+-------------------------------------------------------------------+
  subroutine tblrestress(y,utau,delta,rowall,twall,re,uu,vv,ww,uv)
    !
    use basicfunction, only: miucal
    !
    ! arguments
    real(8),intent(in) :: y(0:)
    real(8),intent(in) :: utau,delta,rowall,twall,re
    real(8),allocatable,intent(out) :: uu(:),vv(:),ww(:),uv(:)
    !
    ! local data
    integer :: j,m,je,i,h
    real(8) :: miu,alfa,beter,var1
    real(8),allocatable :: yplus(:),ydelta(:)
    real(8) :: ares(1:4,0:9),bres(1:4,0:9)
    real(8),allocatable :: restress(:,:),restressin(:,:),restressou(:,:)
    real(8),allocatable :: weig1(:),weig2(:),weig(:)
    !
    m=size(y)-1
    allocate(yplus(0:m),ydelta(0:m))
    allocate(uu(0:m),vv(0:m),ww(0:m),uv(0:m))
    !
    miu=miucal(twall)
    do j=0,m
      yplus(j)=rowall*utau*y(j)/miu*re
      ydelta(j)=y(j)/delta
    enddo
    !
    write(*,"(2(A,F12.7))")'  ** y+_max= ',yplus(m),'; y_max=',ydelta(m)
    !
    ares(1,0)=0.0d0           ;  bres(1,0)=3.29569d0   
    ares(1,1)=0.54433d0       ;  bres(1,1)=-18.42051d0 
    ares(1,2)=-0.03905d0      ;  bres(1,2)=91.29933d0  
    ares(1,3)=0.00135d0       ;  bres(1,3)=-267.43365d0
    ares(1,4)=-2.66693d-5     ;  bres(1,4)=484.83535d0 
    ares(1,5)=3.20392d-7      ;  bres(1,5)=-578.84056d0
    ares(1,6)=-2.37906d-9     ;  bres(1,6)=459.33499d0 
    ares(1,7)=1.06538d-11     ;  bres(1,7)=-233.5512d0 
    ares(1,8)=-2.63502d-14    ;  bres(1,8)=68.86463d0  
    ares(1,9)=2.76277d-17     ;  bres(1,9)=-8.9459d0
    !
    ares(2,0)=-0.0d0          ;  bres(2,0)=0.30039d0   
    ares(2,1)=0.03934d0       ;  bres(2,1)=10.64176d0  
    ares(2,2)=1.14714d-4      ;  bres(2,2)=-66.37375d0 
    ares(2,3)=-3.85848d-5     ;  bres(2,3)=225.86454d0 
    ares(2,4)=1.07969d-6      ;  bres(2,4)=-480.39234d0
    ares(2,5)=-1.51441d-8     ;  bres(2,5)=650.63211d0 
    ares(2,6)=1.22937d-10     ;  bres(2,6)=-560.54251d0
    ares(2,7)=-5.83882d-13    ;  bres(2,7)=297.57773d0 
    ares(2,8)=1.50637d-15     ;  bres(2,8)=-88.75901d0 
    ares(2,9)=-1.63084d-18    ;  bres(2,9)=11.37986d0 
    !
    ares(3,0)=0.0d0           ;  bres(3,0)=1.19685d0   
    ares(3,1)=0.16499d0       ;  bres(3,1)=1.4586d0    
    ares(3,2)=-0.00998d0      ;  bres(3,2)=-14.13002d0 
    ares(3,3)=3.28209d-4      ;  bres(3,3)=58.62495d0  
    ares(3,4)=-6.39718d-6     ;  bres(3,4)=-153.12029d0
    ares(3,5)=7.68168d-8      ;  bres(3,5)=237.76412d0 
    ares(3,6)=-5.72803d-10    ;  bres(3,6)=-221.67287d0
    ares(3,7)=2.57995d-12     ;  bres(3,7)=122.62448d0 
    ares(3,8)=-6.41997d-15    ;  bres(3,8)=-37.20488d0 
    ares(3,9)=6.77108d-18     ;  bres(3,9)=4.77732d0   
    !
    ares(4,0)=0.0d0           ;  bres(4,0)=-0.77067d0    
    ares(4,1)=-0.06311d0      ;  bres(4,1)=-3.15123d0    
    ares(4,2)=2.71427d-4      ;  bres(4,2)=27.29017d0    
    ares(4,3)=7.09745d-5      ;  bres(4,3)=-108.34305d0  
    ares(4,4)=-2.42809d-6     ;  bres(4,4)=270.30229d0   
    ares(4,5)=3.78387d-8      ;  bres(4,5)=-414.51648d0  
    ares(4,6)=-3.28545d-10    ;  bres(4,6)=387.55797d0   
    ares(4,7)=1.63456d-12     ;  bres(4,7)=-215.80852d0  
    ares(4,8)=-4.36111d-15    ;  bres(4,8)=65.89566d0    
    ares(4,9)=4.84169d-18     ;  bres(4,9)=-8.50571d0
    !
    allocate(restress(1:4,0:m),restressin(1:4,0:m),restressou(1:4,0:m))
    !
    restressin=0.d0
    !
    do j=0,m
      !
      if(yplus(j)<=150.d0) then
        !
        do i=1,4
          !
          do h=0,9
            !
            restressin(i,j)=restressin(i,j)+ares(i,h)*yplus(j)**h
            !
          end do
          !
        end do
        !
      endif
      !
    end do
    !
    restressou=0.d0
    do j=0,m
      !
      do i=1,4
        !
        do h=0,9
          !
          restressou(i,j)=restressou(i,j)+bres(i,h)*ydelta(j)**h
          !
        end do
        !
      end do
      !
    end do
    !
    allocate(weig1(0:m),weig2(0:m),weig(0:m))
    !
    alfa=2.d0
    beter=0.01d0
    do j=0,m
      var1=alfa*(ydelta(j)-beter)/((1.d0-2.d0*beter)*ydelta(j)+beter)
      weig1(j)=0.5d0+0.5d0*tanh(var1)/tanh(alfa)
    end do
    !
    alfa=15.d0
    beter=1.d0
    do j=0,m
      var1=(ydelta(j)-beter)*alfa
      weig2(j)=1.d0-(1.d0/(1.d0+dexp(-var1)))
    end do
    !
    do j=0,m
      weig(j)=(weig1(j)*weig2(j))
    end do
    !
    open(18,file='weight.dat')
    write(18,"(4(1X,A20))")'y','weight1','weight2','weight'
    write(18,"(4(1X,E14.7E2))")(yplus(j),weig1(j),weig2(j),weig(j),j=0,m)
    close(18)
    print*,' <<< weight.dat ... done.'
    !
    do j=0,m
      do i=1,4
        restress(i,j)=(1.d0-weig(j))*restressin(i,j)+weig(j)*restressou(i,j)
      end do
    end do
    !
    open(18,file='restress_profile.dat')
    write(18,"(6(1X,A15))")'yplus','ydelta','uu','vv','ww','uv'
    write(18,"(6(1X,E15.7E3))")(yplus(j),ydelta(j),(restress(i,j)**2,i=1,3),restress(4,j),j=0,m)
    close(18)
    print*,' << restress_profile_in.dat'
    !
    do j=0,m
      uu(j)=(restress(1,j)*utau)**2
      vv(j)=(restress(2,j)*utau)**2
      ww(j)=(restress(3,j)*utau)**2
      uv(j)=restress(4,j)*utau**2
    end do
    !
  end subroutine tblrestress
  !+-------------------------------------------------------------------+
  !| The end of the subroutine tblrestress.                            |
  !+-------------------------------------------------------------------+
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to analysis compressible laminar boundary 
  ! layer.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2010-07-07.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine lamianrBL(mode)
    !
    use commvardefine,only: im,jm,km,Reynolds,pinf,Mach,Prandtl,const2,const5
    use basicfunction
    !
    character(len=*),intent(in) :: mode
    !
    real(8),allocatable,dimension(:) :: ye,y2,f,f1,f2,f3,t,t1,t2,u,v,ro,uplus,yplus
                                        
    real(8) :: dye,pre,den,xin,thick1,thick2,thick3,miu
    !
    integer :: ns1,ns2,j1,j,i,ios
    real(8) :: error1,error2,p0,p1,p2,p3,p4,q0,q1,q2,var0,var1,var2,   &
               var3,var4,var5,var6,df2,dt,f210,f211,t00,t01,rom,rom2,utaw,lref
    !
    print*,' >> xin' 
    read(*,*)xin
    !
    if(trim(mode)=='gen') then
      print*,' will solve the blasisu equation'
      !
      jm=260
      !
      allocate( ye(0:jm),y2(0:jm),f(0:jm),f1(0:jm),f2(0:jm),f3(0:jm),    &
                t(0:jm),t1(0:jm),t2(0:jm),u(0:jm),v(0:jm),ro(0:jm)       )
      !
      do j=0,jm
        ye(j)=20.d0/jm*j
      end do
      dye=ye(1)-ye(0)
      !
      ! for boundary thickness=1.d0
      ! xin=20.d0
      !
      pre=pinf
      !
      ! Initial Boundary value
      f(0)=0.d0
      f1(0)=0.d0
      ! shooting
      f2(0)=0.33206d0 
      !
      f3(0)=-1.d0*f(0)*f2(0)
      !
      df2=0.01d0
      !
      ! shooting
      ! t(0)=1.d0+0.5d0*const5*dsqrt(Prandtl)
      t(0)=1.d0
      ! t(0)=0.801d0
      ! t(0)=4.39d0
      !t(0)=1522.44/ref_t
      t(jm)=1.d0
      print*,' ** Twall=',t(0)
      !
      den=pinf/t(0)*const2
      miu=MiuCal(t(0))
      rom=den*miu
      !
      ! adiabatic wall
      ! t1(0)=0.d0
      t1(0)= 1.d0
      ! t2(0)=-const5*rom*f2(0)**2-f(0)*t1(0)
      t2(0)=0.08d0
      !
      dt=0.1d0
      !
      ! starting calculationg
      !
      open(13,file='stat_laminar.dat')
      write(13,"(A7,2(1X,A20))")'n','twall','error'
      error1=1.d0
      ns1=0
      do while(error1>1.d-5)
        ns1=ns1+1
        error1=0.d0
        f210=f1(jm)-1.d0
        t00=t(jm)-1.d0
        !
        error2=1.d0
        ns2=0
        do while(error2>1d-6)
          !
          ns2=ns2+1
          error2=0.d0
          !
          do j=1,jm
            !
            ! Predictor
            miu=MiuCal(t(j-1))
            den=pinf/t(j-1)*const2
            rom=miu*den
            rom2=miu*den/prandtl
            !
            p0= f(j-1)+f1(j-1)*dye
            p1=f1(j-1)+f2(j-1)*dye
            p2=f2(j-1)+f3(j-1)/rom*dye
            p3=-1d0*p0*p2
            !
            q0=t(j-1)+t1(j-1)*dye
            miu=MiuCal(q0)
            den=pinf/q0*const2
            rom=miu*den
            rom2=miu*den/prandtl
            q1=t1(j-1)+t2(j-1)/rom2*dye
            q2=-p0*q1-const5*rom*p2**2
            !
            ! Corrector
            miu=MiuCal(q0)
            den=pinf/q0*const2
            rom=miu*den
            rom2=miu*den/prandtl
            !
            var0= f(j-1)+0.5d0*(f1(j-1)+p1)*dye
            var1=f1(j-1)+0.5d0*(f2(j-1)+p2)*dye
            var2=f2(j-1)+0.5d0*(f3(j-1)+p3)/rom*dye
            var3=-1d0*f(j)*f2(j)
            !
            var4= t(j-1)+0.5d0*(t1(j-1)+q1)*dye
            miu=MiuCal(var4)
            den=pinf/var4*const2
            rom=miu*den
            rom2=miu*den/prandtl
            var5=t1(j-1)+0.5d0*(t2(j-1)+q2)/rom2*dye
            var6=-var0*var5-const5*rom*var2**2
            !
            ! error calculation
            error2=max( error2,dabs(var0- f(j)),dabs(var1-f1(j)),        &
                               dabs(var2-f2(j)),dabs(var3-f3(j)),        &
                               dabs(var4- t(j)),dabs(var5-t1(j)),        &
                               dabs(var6-t2(j))                          )
            !
             f(j)=var0
            f1(j)=var1
            f2(j)=var2
            f3(j)=var3
            !
             t(j)=var4
            t1(j)=var5
            t2(j)=var6
            !
          end do
          !
          !
          !print*,'ns2',ns2,'error2',error2
          !
        end do
        !
        ! shooting f2(0) to satisify f1(jm)=1.d0
        error1=max(error1,dabs(f1(jm)-1.d0))
        !
        f211=f1(jm)-1.d0
        !
        if(f210*f211>0.d0) then
        else
          df2=df2*0.25d0
        end if
        !
        if(f1(jm)<1.d0) then
          f2(0)=f2(0)+df2
        elseif(f1(jm)>1.d0) then
          f2(0)=f2(0)-df2
        end if
        !
        !
        ! shooting t1(0) to satisify t(jm)=1.d0
        error1=max(error1,dabs(t(jm)-1.d0))
        t01=t(jm)-1.d0
        if(t00*t01>0.d0) then
        else
          dt=dt*0.5d0
        end if
        !
        if(t(jm)<1.d0) then
          t1(0)=t1(0)+dt
        elseif(t(jm)>1.d0) then
          t1(0)=t1(0)-dt
        end if
        !
        ! print*,f1(jm),t(jm),df2
        !! shooting t(0) to satisify t(jm)=1.d0
        !error1=max(error1,dabs(t(jm)-1.d0))
        !t01=t(jm)-1.d0
        !!
        !if(t00*t01<0.d0) then
        !  dt=dt*0.5d0
        !end if
        !!
        !if(t(jm)<1.d0) then
        !  t(0)=t(0)+dt
        !elseif(t(jm)>1.d0) then
        !  t(0)=t(0)-dt
        !end if
        !
        write(13,"(I7,2(1X,E20.13E2))")ns1,t(0),error1
        !
        print*,ns1,error1,t1(0),t(jm),dt
        !
      end do
      close(13)
      !
    else
      print*,' will process a data from the file:',trim(mode)
      open(12,file=trim(mode))
      read(12,*)
      jm=-1
      ios=0
      do while(ios==0)
        read(12,*,iostat=ios)
        if(ios==0) jm=jm+1
      enddo
      close(12)
      !
      print*,'jm=',jm
      !
      allocate( ye(0:jm),y2(0:jm),f(0:jm),f1(0:jm),f2(0:jm),f3(0:jm),    &
                t(0:jm),t1(0:jm),t2(0:jm),u(0:jm),v(0:jm),ro(0:jm)       )
      open(12,file=trim(mode))
      read(12,*)
      do j=0,jm
        read(12,*)ye(j),f1(j),t(j)
      enddo
      close(12)
      print*,' >> ',trim(mode)
      !
      dye=ye(1)-ye(0)
      !
    endif
    !
    var0=1.d0/dsqrt(2.d0*xin*Reynolds)
    do j=0,jm
      u(j)=f1(j)
      v(j)=0.d0
      do i=1,j
        var1=0.5d0*(t(i)+t(i-1))*dye
        var2=0.5d0*(u(i)+u(i-1))*dye
        v(j)=v(j)+var0*(u(j)*var1-t(j)*var2)
      end do
    end do
    !
    y2(0)=0.d0
    var1=dsqrt(2.d0*xin/Reynolds)
    do j=1,jm
      !y2(j)=var1*ye(j)
      y2(j)=0.d0
      do i=1,j
        y2(j)=y2(j)+var1*0.5d0*(t(i)+t(i-1))*dye
      end do
    end do
    !
    ! Calculationg of boundary thickness
    thick1=0.d0
    thick2=0.d0
    thick3=0.d0
    !
    var0=u(jm)*0.99d0
    do j=1,jm
      if(u(j-1)<=var0 .and. u(j)>=var0) then
        thick1=(y2(j)-y2(j-1))/(u(j)-u(j-1))*(var0-u(j-1))+y2(j-1)
        j1=j
        exit
      end if
    end do
    !
    ro=1.d0/t
    !
    var0=ro(jm)*u(jm)
    do j=1,jm
      var1=0.5d0*(ro(j)+ro(j-1))
      var2=0.5d0*(u(j)+u(j-1))
      thick2=thick2+(1.d0-var1*var2/var0)*(y2(j)-y2(j-1))
    end do
    !
    var0=ro(jm)*u(jm)**2
    do j=1,jm
      var1=0.5d0*(ro(j)+ro(j-1))
      var2=0.5d0*(u(j)+u(j-1))
      thick3=thick3+(var1*var2*(u(jm)-var2))/var0*(y2(j)-y2(j-1))
    end do
    !
    open(18,file='BoundaryStat.dat')
    write(18,*)'norminal thickness=',thick1
    write(18,*)'displacement thickness=',thick2
    write(18,*)'momentum thickness=',thick3
    write(18,*)'Reynolds Based on δ=',Reynolds*thick1
    write(18,*)'Reynolds Based on δ*=',Reynolds*thick2
    write(18,*)'Reynolds Based on θ=',Reynolds*thick3
    close(18) 
    !
    print*,'norminal thickness=',thick1
    print*,'displacement thickness=',thick2
    print*,'momentum thickness=',thick3
    !
    print*,'Reynolds Based on x=',Reynolds*xin
    print*,'Reynolds Based on δ=',Reynolds*thick1
    print*,'Reynolds Based on δ*=',Reynolds*thick2
    print*,'Reynolds Based on θ=',Reynolds*thick3
    !
    open(18,file='CompBlasius.dat')
    write(18,"(5(1X,A15))")'y','yn','u','v','t'
    do j=0,jm
      write(18,"(5(1X,E15.7E3))")y2(j),y2(j)*dsqrt(Reynolds/xin),u(j),v(j),t(j)
    end do
    close(18)
    print*,' << CompBlasius.dat'
    !
    print*,'Reynolds ',Reynolds
    !
    call upluscal(uplus,yplus,u,y2,ro,t(0),utaw=utaw)
    !
    lref=1.d0
    open(18,file='profile.bl')
    write(18,*)
    write(18,"(4(1X,A15))")'δ','δ*','θ','utaw'
    write(18,"(4(1X,E15.7E3))")thick1/lref,thick2/lref,thick3/lref,utaw/lref
    write(18,"(4(1X,A15))")'y','u','v','t'
    do j=0,jm
      write(18,"(4(1X,E15.7E3))")y2(j)/lref,u(j),v(j),t(j)
    end do
    close(18)
    print*,' >> profile.bl'
    deallocate( ye,f,f1,f2,f3,t,t1,t2,u,v,y2 )
    !
  end subroutine lamianrBL
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine lamianrBL.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
end module boundarylayer
!+---------------------------------------------------------------------+
!| The end of the module boundarylayer.                                |
!+---------------------------------------------------------------------+