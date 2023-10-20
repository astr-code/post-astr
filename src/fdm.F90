module fdm
  !
  use commvardefine
  !
  implicit none
  !
  integer,parameter :: hm=4
  real(8),parameter :: bfacmpld=0.8d0
  !
  contains
  !
  ! this function return the derivative
  function dfcal(vin,dim,scheme) result(vout)
    !
    integer,intent(in) :: dim
    real(8),intent(in) :: vin(-hm:dim+hm)
    character(len=*),intent(in) :: scheme
    real(8) :: vout(0:dim)
    !
    real(8) :: fh(-1:dim)
    !
    integer :: i
    !
    call recons(vin,fh,dim,scheme=scheme,mode='sequential')
    !
    do i=0,dim
      vout(i)=(fh(i)-fh(i-1))
    enddo
    !
    return
    !
  end function dfcal
  !
  ! this function return the reconstructured interface value using compact scheme
  subroutine recons(vin,vout,dim,scheme,mode,cstep)
    !
    use LinearAlgegra
    !
    integer,intent(in) :: dim
    real(8),intent(in) :: vin(-hm:dim+hm)
    real(8),intent(inout) :: vout(-1:dim)
    character(len=*),intent(in) :: scheme,mode
    character(len=*),intent(in),optional :: cstep
    !
    !
    real(8),allocatable :: b(:)
    !
    ! local data
    integer :: i,m
    !
    if(scheme=='eu5') then
      !
      do i=-1,dim
        vout(i)=suw5(vin(i-2:i+2))
      enddo
      !
    elseif(scheme=='mp5') then
      !
      do i=-1,dim
        vout(i)=mp5(vin(i-2:i+2))
      enddo
      !
    elseif(scheme=='mp5ld') then
      !
      do i=-1,dim
        vout(i)=mp5ld(vin(i-2:i+3))
      enddo
      !
    elseif(scheme=='eu5-ld') then
      !
      do i=-1,dim
        vout(i)=suw5ld(vin(i-2:i+3))
      enddo
      !
    elseif(scheme=='cmp5') then
      !
      b=compa_recon_rhs(vin,dim,'cu5',mode=mode)
      vout(1:dim)=qtds(b,'cu5',opt='none')
      !
      vout(0)=vout(dim)
      vout(-1)=vout(dim-1)
      !
      do i=-1,dim
        vout(i)=mp5(vin(i-2:i+2),vout(i))
      enddo
      !
    elseif(scheme=='cmp5ld') then
      !
      b=compa_recon_rhs(vin,dim,'cu5-ld',mode=mode)
      vout(1:dim)=qtds(b,'cu5-ld',bf=bfacmpld,opt='none')
      !
      vout(0)=vout(dim)
      vout(-1)=vout(dim-1)
      !
      do i=-1,dim
        vout(i)=mp5(vin(i-2:i+2),vout(i))
      enddo
      !
    elseif(scheme=='cu5') then
      !
      if(mode=='sequential') then
        !
        b=compa_recon_rhs(vin,dim,scheme,mode=mode)
        vout(1:dim)=qtds(b,scheme,opt='none')
        !
        vout(0)=vout(dim)
        vout(-1)=vout(dim-1)
        !
      elseif(mode=='parallel') then
        !
        b=compa_recon_rhs(vin,dim,scheme,mode=mode)
        !
        if(cstep=='pred') then 
          !
          i=-2
          b(1)=b(1)-0.5d0*suw5(vin(i-2:i+2))
          !
          i=dim+1
          b(size(b))=b(size(b))-1.d0/6.d0*suw5(vin(i-2:i+2))
          !
        elseif(cstep=='corr') then 
          !
          i=-2
          b(1)=b(1)-0.5d0*vout(i)
          !
          i=dim+1
          b(size(b))=b(size(b))-1.d0/6.d0*vout(i)
          !
        endif
        !
        vout(-1:dim)=tds(b,scheme,opt='once')
        !
        ! vout(-1)=suw5(vin(-3:1))
        ! vout(dim)=suw5(vin(dim-2:dim+2))
      endif
      !
    elseif(scheme=='cu5-ld') then
      !
      b=compa_recon_rhs(vin,dim,scheme,mode=mode)
      vout(1:dim)=qtds(b,scheme,opt='once',bf=bfacmpld)
      !
      vout(0)=vout(dim)
      vout(-1)=vout(dim-1)
      !
    else
      stop ' !! scheme not defined @ reconcc'
    endif
    !
    ! print*,' ** interface reconstructured'
    !
    return
    !
  end subroutine recons
  !
  !+-------------------------------------------------------------------+
  !| This function is get the r.h.s. of a compact scheme.              |
  !+-------------------------------------------------------------------+
  function compa_recon_rhs(f,m,scheme,mode) result(b)
    !
    ! arguments
    integer,intent(in) :: m
    real(8),intent(in) :: f(-hm:m+hm)
    character(len=*),intent(in) :: scheme,mode
    !
    real(8),allocatable :: b(:)
    !
    ! local data
    integer :: i
    real(8) :: var1
    !
    if(mode=='sequential') then
      allocate(b(1:m))
    elseif(mode=='parallel') then
      allocate(b(-1:m))
    endif
    !
    if(scheme=='cu5') then
      !
      if(mode=='sequential') then
        do i=1,m
          b(i)=1.d0/18.d0*f(i-1)+19.d0/18.d0*f(i)+5.d0/9.d0*f(i+1)
        enddo
      elseif(mode=='parallel') then
        !
        do i=-1,m
          b(i)=1.d0/18.d0*f(i-1)+19.d0/18.d0*f(i)+5.d0/9.d0*f(i+1)
        enddo
        !
      endif
      !
    elseif(scheme=='cu5-ld') then
      !
      if(mode=='sequential') then
        do i=1,m
          b(i)=(1.d0/18.d0 -1.d0/36.d0*bfacmpld)*f(i-1) + &
               (19.d0/18.d0-9.d0/36.d0*bfacmpld)*f(i)   + &
               (5.d0/9.d0  +9.d0/36.d0*bfacmpld)*f(i+1) + &
                            1.d0/36.d0*bfacmpld *f(i+2)
        enddo
      endif
      !
    else
      stop ' !! scheme not defined @ qtds_rec_rhs'
    endif
    !
    return
    !
  end function compa_recon_rhs
  !+-------------------------------------------------------------------+
  !| The end of the function compa_recon_rhs                           |
  !+-------------------------------------------------------------------+
  !
  function diff5ec(vin,dim) result(vout)
    !
    integer,intent(in) :: dim
    real(8),intent(in) :: vin(-hm:dim+hm)
    real(8) :: vout(0:dim)
    !
    ! local data
    integer :: i
    !
    do i=0,dim
      vout(i)  = -1.d0/30.d0*vin(i-3) + &
                  1.d0/4.d0 *vin(i-2)   &
                           - vin(i-1) + &
                  1.d0/3.d0 *vin(i)   + &
                  1.d0/2.d0 *vin(i+1)   &
                 -1.d0/20.d0*vin(i+2)
    enddo
    !
  end function diff5ec
  !
  function diff6ec(vin,dim) result(vout)
    !
    integer,intent(in) :: dim
    real(8),intent(in) :: vin(-hm:dim+hm)
    real(8) :: vout(0:dim)
    !
    ! local data
    integer :: i
    !
    do i=0,dim
      vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-                         &
                0.15d0 *(vin(i+2)-vin(i-2))+                         &
                1.d0/60.d0*(vin(i+3)-vin(i-3))
    enddo
    !
  end function diff6ec
  !
  pure function suw5(u) result(uh)
    !
    real(8),intent(in) :: u(1:5)
    real(8) :: uh
    !
    uh=(2.d0*u(1)-13.d0*u(2)+47.d0*u(3)+27.d0*u(4)-3.d0*u(5))/60.d0
    !
    return
    !
  end function suw5
  !
  pure function suw5ld(u) result(uh)
    !
    real(8),intent(in) :: u(1:6)
    real(8) :: uh
    real(8) :: vadp,b1,b2,b3,b4,b5,b6
    !
    uh=(2.d0*u(1)-13.d0*u(2)+47.d0*u(3)+27.d0*u(4)-3.d0*u(5))/60.d0
    !
    vadp=1.666666666666667d-2*bfacmpld
    !
    b1=      -vadp+3.333333333333333d-2
    b2=  5.d0*vadp-2.166666666666667d-1
    b3=-10.d0*vadp+7.833333333333333d-1
    b4= 10.d0*vadp+0.45d0
    b5= -5.d0*vadp-5.d-2
    b6=       vadp
    !
    uh=b1*u(1)+b2*u(2)+b3*u(3)+b4*u(4)+b5*u(5)+b6*u(6)
    !
    return
    !
  end function suw5ld
  !
  pure function suw7(u) result(uh)
    !
    real(8),intent(in) :: u(1:7)
    real(8) :: uh
    !
    uh=(-3.d0*u(1)+25.d0*u(2)-101.d0*u(3)+319.d0*u(4)+214.d0*u(5)-     &
        38.d0*u(6)+ 4.d0*u(7))/420.d0
    !
    return
    !
  end function suw7
  !
  pure function MP5(u,uhat)  result(uh)
    !
    real(8),intent(in) :: u(1:5)
    real(8),intent(in),optional :: uhat
    real(8) :: uh
    !
    real(8) :: ulinear,uMP,uUL,uAV,uMD,uLC,uMIN,uMAX
    real(8) :: var1,var2,dm1,d0,d1,dhm1,dh0
    !
    if(present(uhat)) then
      ulinear=uhat
    else
      ulinear=(2.d0*u(1)-13.d0*u(2)+47.d0*u(3)+27.d0*u(4)-3.d0*u(5))/60.d0
    endif
    !
    var1=u(4)-u(3)
    var2=4.d0*(u(3)-u(2))
    uMP=u(3)+minmod2(var1,var2)
    !
    var1=(ulinear-u(3))*(ulinear-uMP)
    if(var1>=1.d-10) then
    !if(switch(i)) then
      dm1=u(1)-2.d0*u(2)+u(3)
      d0= u(2)-2.d0*u(3)+u(4)
      d1= u(3)-2.d0*u(4)+u(5)
      !
      dhm1=minmod4( 4.d0*dm1-d0,4.d0*d0-dm1,dm1,d0 )
      dh0 =minmod4(  4.d0*d0-d1, 4.d0*d1-d0, d0,d1 )
      !
      uUL=u(3)+4.d0*(u(3)-u(2))
      uAV=0.5d0*(u(3)+u(4))
      uMD=uAV-0.5d0*dh0
      uLC=u(3)+0.5d0*(u(3)-u(2))+1.333333333333333d0*dhm1
      !
      var1=min(u(3),u(4),uMD)
      var2=min(u(3),uUL,uLC)
      uMIN=max(var1,var2)
      !
      var1=max(u(3),u(4),uMD)
      var2=max(u(3),uUL,uLC)
      uMAX=min(var1,var2)
      !
      var1=uMIN-ulinear
      var2=uMAX-ulinear
      uh=ulinear+minmod2(var1,var2)
    else
      uh=ulinear
      ! No limiter is needed
    end if
    !
    return
    !
  end function MP5
  !
  function MP5LD(u) result(uh)
    !
    real(8),intent(in) :: u(1:6)
    real(8) :: uh
    !
    ! local data
    real(8) :: weightBW
    real(8) :: ulinear,uMP,uUL,uAV,uMD,uLC,uMIN,uMAX
    real(8) :: var1,var2,dm1,d0,d1,dhm1,dh0
    !
    real(8) :: b1,b2,b3,b4,b5,b6,vadp
    !
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! vadp: optimize parameter, to control the dissipation
    ! of the linear scheme. 
    ! vadp from 0 to 1/60
    ! vadp=0: standard 5 order upwind. Max dissipation
    ! vadp=1/60: standard 6 order center. 0 dissipation
    ! lskt: shock or not
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    weightBW=0.7d0
    ! if(lsod) then
    !   ulinear=-num1d6*u(2)+num5d6*u(3)+num1d3*u(4)
    !   ! ulinear=u(3)
    ! else
      vadp=1.666666666666667d-2*weightBW
      !vadp=0.0d0*1.666666666666667d-2
      !
      b1=      -vadp+3.333333333333333d-2
      b2=  5.d0*vadp-2.166666666666667d-1
      b3=-10.d0*vadp+7.833333333333333d-1
      b4= 10.d0*vadp+0.45d0
      b5= -5.d0*vadp-5.d-2
      b6=       vadp
      !
      ulinear=b1*u(1)+b2*u(2)+b3*u(3)+b4*u(4)+b5*u(5)+b6*u(6)
    ! endif
    !!
    var1=u(4)-u(3)
    var2=4.d0*(u(3)-u(2))
    uMP=u(3)+minmod2(var1,var2)
    !
    var1=(ulinear-u(3))*(ulinear-uMP)
    ! if(.false.) then
    if(var1>=1.d-10) then
      dm1=u(1)-2.d0*u(2)+u(3)
      d0= u(2)-2.d0*u(3)+u(4)
      d1= u(3)-2.d0*u(4)+u(5)
      !
      dhm1=minmod4( 4.d0*dm1-d0,4.d0*d0-dm1,dm1,d0 )
      dh0 =minmod4(  4.d0*d0-d1, 4.d0*d1-d0, d0,d1 )
      !
      uUL=u(3)+4.d0*(u(3)-u(2))
      uAV=0.5d0*(u(3)+u(4))
      uMD=uAV-0.5d0*dh0
      uLC=u(3)+0.5d0*(u(3)-u(2))+1.333333333333333d0*dhm1
      !
      var1=min(u(3),u(4),uMD)
      var2=min(u(3),uUL,uLC)
      uMIN=max(var1,var2)
      !
      var1=max(u(3),u(4),uMD)
      var2=max(u(3),uUL,uLC)
      uMAX=min(var1,var2)
      !
      var1=uMIN-ulinear
      var2=uMAX-ulinear
      uh=ulinear+minmod2(var1,var2)
    else
      uh=ulinear
      !
    end if
    !
    return
    !
  end function MP5LD
  !
  pure function MP7(u) result(uh)
    !
    real(8),intent(in) :: u(0:6)
    real(8) :: uh
    !
    real(8) :: ulinear,uMP,uUL,uAV,uMD,uLC,uMIN,uMAX
    real(8) :: var1,var2,dm1,d0,d1,dhm1,dh0
    !
    ulinear=(-3.d0*u(0)+25.d0*u(1)-101.d0*u(2)+319.d0*u(3)+          &
               214.d0*u(4)-38.d0*u(5)+ 4.d0*u(6))/420.d0
    !
    var1=u(4)-u(3)
    var2=4.d0*(u(3)-u(2))
    uMP=u(3)+minmod2(var1,var2)
    !
    var1=(ulinear-u(3))*(ulinear-uMP)
    if(var1>=1.d-10) then
      dm1=u(1)-2.d0*u(2)+u(3)
      d0= u(2)-2.d0*u(3)+u(4)
      d1= u(3)-2.d0*u(4)+u(5)
      !
      dhm1=minmod4( 4.d0*dm1-d0,4.d0*d0-dm1,dm1,d0 )
      dh0 =minmod4(  4.d0*d0-d1, 4.d0*d1-d0, d0,d1 )
      !
      uUL=u(3)+4.d0*(u(3)-u(2))
      uAV=0.5d0*(u(3)+u(4))
      uMD=uAV-0.5d0*dh0
      uLC=u(3)+0.5d0*(u(3)-u(2))+1.333333333333333d0*dhm1
      !
      var1=min(u(3),u(4),uMD)
      var2=min(u(3),uUL,uLC)
      uMIN=max(var1,var2)
      !
      var1=max(u(3),u(4),uMD)
      var2=max(u(3),uUL,uLC)
      uMAX=min(var1,var2)
      !
      var1=uMIN-ulinear
      var2=uMAX-ulinear
      uh=ulinear+minmod2(var1,var2)
    else
      uh=ulinear
    end if
    !
    return
    !
  end function MP7
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this fuction is the 3-variable median() function.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure function median(var1,var2,var3)
    !
    real(8),intent(in) :: var1,var2,var3
    real(8) :: median
    !
    median= var1+minmod2(var2-var1,var3-var1)
    !
    return
    !
  end function median
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function median.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this fuction is the 2-variable minmod() function.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure function minmod2(var1,var2)
    !
    real(8),intent(in) :: var1,var2
    real(8) :: minmod2
    !
    if(var1>0.d0 .and. var2>0.d0) then
      minmod2=min(dabs(var1),dabs(var2))
    elseif(var1<0.d0 .and. var2<0.d0) then
      minmod2=-1.d0*min(dabs(var1),dabs(var2))
    else
      minmod2=0.d0
    end if
    !
    return
    !
  end function minmod2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function minmod2.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this fuction is the 4-variable minmod() function.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure function minmod4(var1,var2,var3,var4)
    !
    real(8),intent(in) :: var1,var2,var3,var4
    real(8) :: minmod4
    !
    if(var1>0.d0 .and. var2>0.d0 .and. var3>0.d0 .and. var4>0.d0) then
      minmod4=min(dabs(var1),dabs(var2),dabs(var3),dabs(var4))
    elseif(var1<0.d0 .and. var2<0.d0 .and. var3<0.d0 .and. var4<0.d0)&
                                                                  then
      minmod4=-1.d0*min(dabs(var1),dabs(var2),dabs(var3),dabs(var4))
    else
      minmod4=0.d0
    end if
    !
    return
    !
  end function minmod4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function minmod4.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine rk4(f,rhs,deltat,rkstep)
    !
    real(8),intent(inout) :: f(:)
    real(8),intent(in) :: rhs(:),deltat
    integer,intent(in) :: rkstep
    !
    ! local data
    integer :: dim
    logical,save :: firstcall = .true.
    real(8),save :: rkcoe(2,4)
    real(8),allocatable,save :: fsave(:),rhsav(:)
    !
    if(firstcall) then
      !
      rkcoe(1,1)=0.5d0
      rkcoe(1,2)=0.5d0
      rkcoe(1,3)=1.d0
      rkcoe(1,4)=1.d0/6.d0
      !
      rkcoe(2,1)=1.d0
      rkcoe(2,2)=2.d0
      rkcoe(2,3)=2.d0
      rkcoe(2,4)=1.d0
      !
      dim=size(f)
      !
      allocate(fsave(1:dim),rhsav(1:dim))
      !
      firstcall=.false.
      !
    endif
    !
    if(rkstep==1) then
      fsave=f
      rhsav=0.d0
      !
    endif
    !
    if(rkstep<=3) then
      f(:)=fsave(:)+rkcoe(1,rkstep)*deltat*rhs(:)
      !
      rhsav(:)=rhsav(:)+rkcoe(2,rkstep)*rhs(:)
    elseif(rkstep==4) then
      f(:)=fsave(:)+rkcoe(1,rkstep)*deltat*(rhs(:)+rhsav(:))
    else
      print*,rkstep
      stop ' !! error rkstep @ rk4'
    endif
    !
  end subroutine rk4
  !
end module fdm