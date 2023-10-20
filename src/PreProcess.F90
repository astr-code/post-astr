!+---------------------------------------------------------------------+
!| This module contains subroutines for preprocess.                    |
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2020-04-16                                     |
!+---------------------------------------------------------------------+
module preprocess
  !
  use h5readwrite
  use writetec
  !
  implicit none
  !
  contains
  !
  subroutine onedflame_process(infilename)
    !
    use h5readwrite
    use basicfunction
    !
    character(len=*),intent(in) :: infilename
    !
    integer :: im,jm,km,dims(3),num_species
    real(8),allocatable,dimension(:,:) :: var
    real(8),allocatable,dimension(:) :: x,ro,u1,p,t
    real(8),allocatable,dimension(:,:) :: spc
    !
    integer :: i,j,k,jsp,nlen
    character(len=10) :: strmp
    character(len=3) :: spname
    character(len=32) :: databasename
    character(len=4096) :: chem
    !
    dims=h5_getdimensio('x','datin/grid.h5')
    im=dims(1)-1
    jm=dims(2)-1
    km=dims(3)-1
    print*,' ** dimension of grid: ',im,jm,km
    !
    num_species=11
    !
    allocate(x(0:im))
    !
    allocate(var(0:im,0:jm))
    call h5_read2dfrom3d(var,im,jm,km,'x','datin/grid.h5',kslice=0)
    x(:)=var(:,0)
    !
    allocate(ro(0:im),u1(0:im),p(0:im),t(0:im))
    !
    call H5ReadArray(ro,im,'ro',infilename)
    call H5ReadArray(u1,im,'u1',infilename)
    call H5ReadArray( t,im, 't',infilename)
    call H5ReadArray( p,im, 'p',infilename)
    !
    allocate(spc(0:im,1:num_species))
    do jsp=1,num_species
      write(spname,'(i3.3)') jsp
      call H5ReadArray(spc(:,jsp),im,'sp'//spname,infilename)
    enddo
    !
    open(12,file='chemrep.dat')
    read(12,*)
    read(12,*)
    read(12,*)
    !
    nlen=1
    !
    do jsp=1,num_species
      !
      write(spname,'(i3.3)') jsp
      read(12,*)i,strmp
      !
      strmp=trim(strmp)//', '
      !
      chem(nlen:nlen+len(trim(strmp)))=trim(strmp)
      !
      nlen=nlen+len(trim(strmp))+1
      !
    enddo
    close(12)
    print*,' >> chemrep.dat'
    !
    print*,trim(chem(1:nlen))
    !
    databasename='h2flame_ratio0.6_1atm.h5'
    call H5WriteArray(num_species,'num_species',trim(databasename))
    call H5WriteArray( x,im, 'x',trim(databasename))
    call H5WriteArray(ro,im,'ro',trim(databasename))
    call H5WriteArray(u1,im,'u1',trim(databasename))
    call H5WriteArray( t,im, 't',trim(databasename))
    call H5WriteArray( p,im, 'p',trim(databasename))
    !
    call H5WriteArray(chem,'chem_list',trim(databasename))
    do jsp=1,num_species
      write(spname,'(i3.3)') jsp
      call H5WriteArray(spc(:,jsp),im,'sp'//spname,trim(databasename))
    enddo
    !
    ! u1=spafilter10_basic(u1)
    !
    open(18,file='profile.dat')
    write(18,"(5(1X,A20))")'x','ro','u1','T','p'
    write(18,"(5(1X,E20.12E3))")(x(i),ro(i),u1(i),t(i),p(i),i=0,im)
    close(18)
    print*,' << profile.dat ... done.'
    !
  end subroutine onedflame_process
  !
  subroutine onedflame_int(infilename)
    !
    use h5readwrite
    use basicfunction
    use interpolation
    !
    character(len=*),intent(in) :: infilename
    !
    integer :: im,im2,jm,km,dims(3),num_species
    real(8),allocatable,dimension(:) :: x,ro,u1,p,t
    real(8),allocatable,dimension(:,:) :: spc
    !
    integer :: i,j,k,jsp
    real(8) :: x0,var1,var2,var3
    character(len=3) :: spname
    character(len=32) :: databasename
    !
    real(8),allocatable,dimension(:) :: x2,ro2,u12,p2,t2
    real(8),allocatable,dimension(:,:) :: spc2
    !
    dims=h5_getdimensio('x',infilename)
    im=dims(1)-1
    !
    allocate(x(0:im),ro(0:im),u1(0:im),p(0:im),t(0:im))
    !
    call H5ReadArray( x,im, 'x',infilename)
    call H5ReadArray(ro,im,'ro',infilename)
    call H5ReadArray(u1,im,'u1',infilename)
    call H5ReadArray( t,im, 't',infilename)
    call H5ReadArray( p,im, 'p',infilename)
    !
    call H5ReadArray(num_species,'num_species',infilename)
    !
    allocate(spc(0:im,1:num_species))
    do jsp=1,num_species
      write(spname,'(i3.3)') jsp
      call H5ReadArray(spc(:,jsp),im,'sp'//spname,infilename)
    enddo
    !
    im2=2*im
    allocate(x2(0:im2),ro2(0:im2),u12(0:im2),p2(0:im2),t2(0:im2))
    !
    x0=0.01d0
    !
    x=x-x0
    !
    open(18,file='profile_ref.dat')
    write(18,"(5(1X,A20))")'x','ro','u1','T','p'
    write(18,"(5(1X,E20.12E3))")(x(i),ro(i),u1(i),t(i),p(i),i=0,im)
    close(18)
    print*,' << profile_ref.dat ... done.'
    !
    do i=0,im2
      x2(i)=0.025d0/dble(im2)*dble(i)
    enddo
    !
    call regularlinearinterp(x,ro,im, x2,ro2,im2, mode='|->')
    call regularlinearinterp(x,u1,im, x2,u12,im2, mode='|->')
    call regularlinearinterp(x, p,im, x2, p2,im2, mode='|->')
    call regularlinearinterp(x, t,im, x2, t2,im2, mode='|->')
    ! u1=spafilter10_basic(u1)
    !
    open(18,file='profile.dat')
    write(18,"(5(1X,A20))")'x','ro','u1','T','p'
    write(18,"(5(1X,E20.12E3))")(x2(i),ro2(i),u12(i),t2(i),p2(i),i=0,im2)
    close(18)
    print*,' << profile.dat ... done.'
    !
  end subroutine onedflame_int
  !
  subroutine onedflame_turb_int(infilename)
    !
    use h5readwrite
    use basicfunction
    use interpolation
    !
    character(len=*),intent(in) :: infilename
    !
    integer :: im,im2,jm,jm2,km,km2,dims(3),num_species
    real(8),allocatable,dimension(:,:,:) :: x,y,z,x2,y2,z2
    !
    real(8),allocatable,dimension(:) :: rop,u1p,pp,tp
    real(8),allocatable,dimension(:,:) :: spcp
    !
    integer :: i,j,k,jsp
    real(8) :: xflame,dflame,x0,y1max,y2max,scale
    real(8) :: uref,pref,roref,tref,var1,var2,var3,xshift,xs
    character(len=3) :: spname
    character(len=32) :: databasename
    !
    real(8),allocatable,dimension(:,:,:) :: ro,u,v,w,p,t
    real(8),allocatable,dimension(:,:,:) :: u2,v2,w2,ro2,p2,t2
    real(8),allocatable,dimension(:,:,:,:) :: spc
    !
    xflame=0.017d0
    dflame=0.0003763d0
    !
    dims=h5_getdimensio('x','datin/grid.h5')
    im=dims(1)-1
    jm=dims(2)-1
    km=dims(3)-1
    !
    allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
    !
    call H5ReadArray(x,im,jm,km,'x','datin/grid.h5')
    call H5ReadArray(y,im,jm,km,'y','datin/grid.h5')
    call H5ReadArray(z,im,jm,km,'z','datin/grid.h5')
    !
    y1max=y(0,jm,0)
    !
    allocate(rop(0:im),u1p(0:im),pp(0:im),tp(0:im))
    !
    call H5ReadArray(rop,im,'ro','datin/flowini1d.h5')
    call H5ReadArray(u1p,im,'u1','datin/flowini1d.h5')
    call H5ReadArray( tp,im, 't','datin/flowini1d.h5')
    call H5ReadArray( pp,im, 'p','datin/flowini1d.h5')
    !
    num_species=11
    !
    allocate(spcp(0:im,1:num_species))
    do jsp=1,num_species
      write(spname,'(i3.3)') jsp
      call H5ReadArray(spcp(:,jsp),im,'sp'//spname,'datin/flowini1d.h5')
    enddo
    !
    dims=h5_getdimensio('x',infilename//'/datin/grid.h5')
    im2=dims(1)-1
    jm2=dims(2)-1
    km2=dims(3)-1
    !
    allocate(x2(0:im2,0:jm2,0:km2),y2(0:im2,0:jm2,0:km2),z2(0:im2,0:jm2,0:km2))
    call H5ReadArray(x2,im2,jm2,km2,'x',infilename//'datin/grid.h5')
    call H5ReadArray(y2,im2,jm2,km2,'y',infilename//'datin/grid.h5')
    call H5ReadArray(z2,im2,jm2,km2,'z',infilename//'datin/grid.h5')
    !
    y2max=y2(0,jm2,0)
    scale=y1max/y2max
    !
    x2=x2*scale
    y2=y2*scale
    z2=z2*scale
    !
    allocate(u2(0:im2,0:jm2,0:km2),v2(0:im2,0:jm2,0:km2),w2(0:im2,0:jm2,0:km2))
    allocate(p2(0:im2,0:jm2,0:km2),t2(0:im2,0:jm2,0:km2),ro2(0:im2,0:jm2,0:km2))
    call H5ReadArray(u2,im2,jm2,km2,'u1',infilename//'outdat/flowfield.h5')
    call H5ReadArray(v2,im2,jm2,km2,'u2',infilename//'outdat/flowfield.h5')
    call H5ReadArray(w2,im2,jm2,km2,'u3',infilename//'outdat/flowfield.h5')
    call H5ReadArray(p2,im2,jm2,km2,'p',infilename//'outdat/flowfield.h5')
    call H5ReadArray(t2,im2,jm2,km2,'t',infilename//'outdat/flowfield.h5')
    call H5ReadArray(ro2,im2,jm2,km2,'ro',infilename//'outdat/flowfield.h5')
    !
    print*,' x2 max=',x2(im2,0,0)
    !
    ! x0=xflame-1.d0*dflame
    x0=0.019d0
    xshift=x2(600,0,0)-x0
    !
    x2=x2-xshift
    print*, ' ** the turbulence box shift in x:',xshift
    !
    xs=5.d0*scale-xshift
    print*, ' ** the turbulence force applied at x:',xs
    !
    uref =u1p(0)
    roref=rop(0)
    tref =tp(0)
    pref=roref*uref**2
    print*,' **     reference velocity:',uref,' m/s'
    print*,' **      reference density:',roref,' kg/m^3'
    print*,' **  reference temperature:',tref,' K'
    print*,' **     reference pressure:',pref,' N/m^2'
    !
    do i=0,im2
      !
      var1=0.d0
      var2=0.d0
      var3=0.d0
      do k=1,km2
      do j=1,jm2
        var1=var1+p2(i,j,k)
        var2=var2+t2(i,j,k)
        var3=var3+ro2(i,j,k)
      enddo
      enddo
      !
      p2(i,:,:) =p2(i,:,:) -var1/dble(jm2*km2)
      t2(i,:,:) =t2(i,:,:) -var2/dble(jm2*km2)
      ro2(i,:,:)=ro2(i,:,:)-var3/dble(jm2*km2)
      !
    enddo
    !
    u2 =(u2-1.d0)*uref
    v2 =v2*uref
    w2 =w2*uref
    p2 =p2*pref
    t2 =t2*tref
    ro2=ro2*roref
    !
    allocate(u(0:im,0:jm,0:km),v(0:im,0:jm,0:km),w(0:im,0:jm,0:km))
    allocate(ro(0:im,0:jm,0:km),t(0:im,0:jm,0:km),p(0:im,0:jm,0:km),spc(0:im,0:jm,0:km,1:num_species))
    !
    do k=0,km
    do j=0,jm
      call regularlinearinterp(x2(:,0,0),u2(:,j,k),im2, x(:,0,0),u(:,j,k),im,   mode='|-|')
      call regularlinearinterp(x2(:,0,0),v2(:,j,k),im2, x(:,0,0),v(:,j,k),im,   mode='|-|')
      call regularlinearinterp(x2(:,0,0),w2(:,j,k),im2, x(:,0,0),w(:,j,k),im,   mode='|-|')
      call regularlinearinterp(x2(:,0,0),ro2(:,j,k),im2, x(:,0,0),ro(:,j,k),im, mode='|-|')
      call regularlinearinterp(x2(:,0,0),p2(:,j,k),im2, x(:,0,0),p(:,j,k),im,   mode='|-|')
      call regularlinearinterp(x2(:,0,0),t2(:,j,k),im2, x(:,0,0),t(:,j,k),im,   mode='|-|')
    enddo
    enddo
    !
    do i=0,im
      ro(i,:,:)=rop(i) + ro(i,:,:)
      u(i,:,:) =u1p(i) +  u(i,:,:)
      p(i,:,:) =pp(i)  +  p(i,:,:)
      t(i,:,:) =tp(i)  +  t(i,:,:)
      !
      do jsp=1,num_species
        spc(i,:,:,jsp)=spcp(i,jsp)
      enddo
    enddo
    !
    call writetecbin('tecini.plt',x,'x',y,'y',z,'z',ro,'ro',u,'u', &
                                    v,'v',w,'w',p,'p',t,'T',im,jm,km)
    call writetecbin('tecini2d.plt',x(:,:,km/2),'x',y(:,:,km/2),'y', &
                                    z(:,:,km/2),'z',u(:,:,km/2),'u', &
                                    v(:,:,km/2),'v',w(:,:,km/2),'w', &
                                    ro(:,:,km/2),'ro',p(:,:,km/2),'p', &
                                    t(:,:,km/2),'T', im,jm)
    !
    call H5WriteArray(ro,im,jm,km,'ro','flowini3d.h5')
    call H5WriteArray( u,im,jm,km,'u1','flowini3d.h5')
    call H5WriteArray( v,im,jm,km,'u2','flowini3d.h5')
    call H5WriteArray( w,im,jm,km,'u3','flowini3d.h5')
    call H5WriteArray( p,im,jm,km,'p','flowini3d.h5')
    call H5WriteArray( t,im,jm,km,'t','flowini3d.h5')
    do jsp=1,num_species
      write(spname,'(i3.3)') jsp
      call H5WriteArray(spc(:,:,:,jsp),im,jm,km,'sp'//spname,'flowini3d.h5')
    enddo
    ! call regularlinearinterp(x,u1,im, x2,u12,im2, mode='|->')
    ! call regularlinearinterp(x, p,im, x2, p2,im2, mode='|->')
    ! call regularlinearinterp(x, t,im, x2, t2,im2, mode='|->')
    ! ! u1=spafilter10_basic(u1)
    ! !
    ! open(18,file='profile.dat')
    ! write(18,"(5(1X,A20))")'x','ro','u1','T','p'
    ! write(18,"(5(1X,E20.12E3))")(x2(i),ro2(i),u12(i),t2(i),p2(i),i=0,im2)
    ! close(18)
    ! print*,' << profile.dat ... done.'
    !
  end subroutine onedflame_turb_int
  !
  ! ref: Li, Z., Jaberi, F. 2010. Numerical Investigations of        
  !      Shock-Turbulence Interaction in a Planar Mixing Layer. 
  ! Kourta, A., Sauvage, R. 2002 Computation of supersonic mixing 
  ! layers. PHYS FLUIDS. 14,3790-7.
  subroutine inletmixing
    !
    use commvardefine, only: jm,Reynolds,mach,gamma,const2,const5,pi,Prandtl
    !
    real(8),allocatable :: y(:),u(:),v(:),t(:),ro(:),p(:),uu(:),vv(:),ww(:),uv(:)
    !
    real(8) :: u1,u2,delta0,c,mc,ustar,uinf,var1
    integer :: j,j1
    !
    allocate(y(0:jm),u(0:jm),v(0:jm),t(0:jm),p(0:jm),ro(0:jm))
    !
    open(12,file='gridy.dat')
    do j=0,jm
      read(12,*)j1,y(j)
    end do
    close(12)
    print*,' >> gridy.dat ... done.'
    !
    u1=3.5d0
    u2=1.5d0
    delta0=1.d0
    mc=(u1-u2)/2.d0*mach
    ! delta0=0.01d0
    !
    uinf=0.5d0*(u1-u2)
    !
    print*,'mc=',mc
    !
    do j=0,jm
      !
      u(j)=0.5d0*(u1+u2)+0.5d0*(u1-u2)*tanh(2.d0*y(j)/delta0)
      v(j)=0.d0
      ustar=u(j)-0.5d0*(u1+u2)
      !
      t(j)=1.d0*(1.d0+mc**2*0.5d0*(gamma-1.d0)*(1.d0-(ustar/uinf)**2))
      ro(j)=1.d0/t(j)
      !
    enddo
    !
    open(16,file='inlet.prof')    
    write(16,"(A26)")'# parameters of inlet flow'
    write(16,"(4(1X,A14))")'delta','delta*','theter','utaw'
    write(16,"(4(1X,E14.7E2))")0.d0,0.d0,0.d0,0.d0
    write(16,"(4(1X,A14))")'ro','u1','u2','t'
    write(16,"(4(1X,E14.7E2))")(ro(j),u(j),v(j),t(j),j=0,jm)
    close(16)
    print*,' <<< inlet.prof ... done!'
    !
    allocate(uu(0:jm),vv(0:jm),ww(0:jm),uv(0:jm))
    !
    c=5.d0
    do j=0,jm
      var1=exp(-0.5d0*(y(j)/c)**10)
      uu(j)=0.15d0**2*var1
      vv(j)=0.15d0**2*var1
      ww(j)=0.15d0**2*var1
      uv(j)=0.15d0**2*0.2d0*var1
    enddo
    !
    open(18,file='inlet.res')
    write(18,"(4(1X,A23))")'#                   R11','R22','R33','R12'
    write(18,"(4(1X,E23.15E3))")(uu(j),vv(j),ww(j),uv(j),j=0,jm)
    close(18)
    print*,' <<< inlet.res ... done.'
    !
    open(18,file='inlet.dat')
      write(18,"(6(1X,A15))")'y','ro','u','v','T','uu'
    do j=0,jm
      write(18,"(6(1X,E15.7E3))")y(j),ro(j),u(j),v(j),t(j),uu(j)
    end do
    close(18)
    print*,' << inlet.dat ... done.'
    !
  end subroutine inletmixing
  !
  subroutine inlet_supchan(database)
    !
    use commvardefine, only: gridfile,jm,Reynolds,mach,gamma,const2,   &
                             const5,pi,Prandtl
    use interpolation
    use flowanalyse, only: PostShockCal
    use basicfunction
    !
    character(len=*),intent(in) :: database
    !
    integer :: jn,j1,j2,j,dims(3)
    real(8) :: var1,ro_ps,u_ps,v_ps,p_ps,t_ps,yshock,xshock,alfa,delta,hstep
    real(8) :: thick1,thick2,thick3,utaw,roinf,uinf,vinf,tinf,pinf
    real(8) :: u_ref,l_ref,ro_ref,t_ref
    real(8),allocatable :: yin(:),uin(:),vin(:),tin(:),roin(:),pin(:)
    real(8),allocatable :: y(:),u(:),v(:),t(:),ro(:),p(:),yplus(:),uplus(:)
    real(8),allocatable :: y2d(:,:)
    !
    ! dims=h5_getdimensio('y',database//'datin/grid.h5')
    ! jn=dims(1)-1
    ! print*,' ** dimension of the input data: ',jn
    ! !
    ! print*,' ** reading mean profile data... '
    ! allocate(yin(0:jn),uin(0:jn),vin(0:jn),tin(0:jn),roin(0:jn),pin(0:jn))
    ! call H5ReadArray(roin,jn,'ro',database//'Results/meanprofile.h5')
    ! call H5ReadArray( uin,jn,'u1',database//'Results/meanprofile.h5')
    ! call H5ReadArray( vin,jn,'u2',database//'Results/meanprofile.h5')
    ! call H5ReadArray( tin,jn, 't',database//'Results/meanprofile.h5')
    ! pin=roin*tin/const2
    ! !
    ! open(12,file=database//'Results/meanprofile.dat')
    ! read(12,*)
    ! read(12,*)thick1,thick2,thick3,utaw
    ! read(12,*)
    ! do j=0,jn
    !   read(12,*)yin(j)
    ! enddo
    ! close(12)
    ! print*,' >> ',database,'meanprofile.dat'
    !
    jn=260
    allocate(yin(0:jn),uin(0:jn),vin(0:jn),tin(0:jn),roin(0:jn),pin(0:jn))
    open(12,file='profile.bl')
    read(12,*)
    read(12,*)
    read(12,*)thick1,thick2,thick3,utaw
    read(12,*)
    do j=0,jn
      read(12,*)yin(j),uin(j),vin(j),tin(j)
    end do
    close(12)
    print*,' >> profile.bl'
    roin=1.d0/tin
    !
    delta=thick1
    !
    l_ref=delta/thick1
    u_ref=952.2174d0
    ro_ref=0.1736504d0
    t_ref=998.7451d0
    !
    ! delta=1.466d0
    ! !
    yin=yin*l_ref
    !
    hstep=0.00305d0
    ! yin=yin+3.05d0
    yin=yin+hstep
    !
    ! jm=320
    !
    dims=h5_getdimensio('y','datin/grid.h5')
    dims=dims-1
    jm=dims(2)
    print*,' ** dimension of the input data: ',jm
    !
    allocate(y(0:jm),u(0:jm),v(0:jm),t(0:jm),p(0:jm),ro(0:jm))
    !
    allocate(y2d(0:dims(1),0:dims(2)))
    ! call H5ReadArray(y2d,dims(1),dims(2),'y','datin/grid.h5')
    call h5_read2dfrom3d(y2d,dims(1),dims(2),dims(3),'y','datin/grid.h5',kslice=0)
    !
    y(:)=y2d(0,:)
    !
    ! open(12,file='gridy.dat')
    ! y(0)=0.d0
    ! do j=1,jm
    !   read(12,*)j1,y(j)
    ! end do
    ! close(12)
    ! print*,' >> gridy.dat ... done.'
    !
    call regularlinearinterp(yin,roin,uin,vin,tin,jn,                  &
                             y,ro,u,v,t,jm)
    !
    uinf=1.d0
    vinf=0.d0
    tinf=1.d0
    roinf=1.d0
    pinf=roinf*tinf/const2
    !
    ! do j=0,jm
    !   if(y(j)<hstep) then
    !     ro(j)=roin(0)
    !      u(j)= uin(0)
    !      v(j)= vin(0)
    !      t(j)= tin(0)
    !   endif
    ! enddo
    ! !
    ! do j=0,jm
    !   if(y(j)>hstep+2.d0*delta) then
    !     !
    !     j2=j
    !     !
    !     exit
    !   endif
    ! enddo
    !
    do j=0,jm
      if((y(j)-hstep)>2.d0*delta) then
        !
        j2=j
        !
        exit
      endif
    enddo
    !
    print*,' edge of 2*delta=',j2,y(j2)
    !
    ! app a damping function
    do j=j2,jm
      var1=(y(j)-y(j2))**2
      var1=exp(-1.d0*var1)
      !
      ro(j)=(ro(j)-roinf)*var1+roinf
      u(j)=(u(j)-uinf)*var1+uinf
      v(j)=(v(j)-vinf)*var1+vinf
      t(j)=(t(j)-tinf)*var1+tinf
      !
    enddo
    ! !
    jn=120
    j2=(jm-jn)/2+jn
    !
    print*,'j2=',j2
    !
    do j=jm,j2,-1
      ro(j)=ro(jm-j+jn)
      u(j)=u(jm-j+jn)
      v(j)=-v(jm-j+jn)
      t(j)=t(jm-j+jn)
    enddo
    ! do j=jm,jm/2,-1
    !   ro(j)=ro(jm-j)
    !   u(j)=u(jm-j)
    !   v(j)=-v(jm-j)
    !   t(j)=t(jm-j)
    !   print*,j,ro(j),t(j)
    ! enddo
    !
    ro=ro*ro_ref
    u=u*u_ref
    v=v*u_ref
    t=t*t_ref
    !
    open(18,file='inlet.dat')
    write(18,"(5(1X,A15))")'y','ro','u','v','T'
    do j=0,jm
      write(18,"(5(1X,E15.7E3))")y(j),ro(j),u(j),v(j),t(j)
    end do
    close(18)
    print*,' << inlet.dat ... done.'
    !
    open(16,file='inlet.prof')    
    write(16,"(A26)")'# parameters of inlet flow'
    write(16,"(4(1X,A14))")'delta','delta*','theter','utaw'
    write(16,"(4(1X,E14.7E2))")thick1*l_ref,thick2*l_ref,thick3*l_ref,utaw*u_ref
    write(16,"(4(1X,A14))")'ro','u1','u2','t'
    write(16,"(4(1X,E14.7E2))")(ro(j),u(j),v(j),t(j),j=0,jm)
    close(16)
    print*,' <<< inlet.prof ... done!'
    !
  end subroutine inlet_supchan
  !
  subroutine inletprofile(database)
    !
    use commvardefine, only: gridfile,jm,Reynolds,mach,gamma,const2,   &
                             const5,pi,Prandtl
    use interpolation
    use flowanalyse, only: PostShockCal
    use basicfunction
    !
    character(len=*),intent(in) :: database
    !
    integer :: jn,j1,j2,j,dims(3)
    real(8) :: var1,ro_ps,u_ps,v_ps,p_ps,t_ps,yshock,xshock,alfa
    real(8) :: ufree,tfree,twall,twall0,thick1,thick2,thick3,utaw
    real(8) :: roinf,uinf,vinf,tinf,pinf
    real(8),allocatable :: yin(:),uin(:),vin(:),tin(:),roin(:),pin(:)
    real(8),allocatable :: y(:),u(:),v(:),t(:),ro(:),p(:),yplus(:),uplus(:)
    real(8),allocatable :: y2d(:,:)
    ! !
    ! dims=h5_getdimensio('x',database//'./datin/grid.2d')-1
    ! !
    ! jn=dims(2)
    
    ! print*,dims,'-',jn
    
    ! allocate(yin(0:jn),uin(0:jn),vin(0:jn),tin(0:jn),roin(0:jn),pin(0:jn))
    ! !
    ! call h5_read1dfrom3d(yin,dims(1),dims(2),dims(3),'y','datin/grid.h5',islice=0,kslice=-1)
    ! !
    ! open(12,file='datin/inlet.prof')
    ! read(12,*)
    ! read(12,*)
    ! read(12,*)
    ! read(12,*)
    ! read(12,*)(roin(j),uin(j),vin(j),tin(j),j=0,jn)
    ! close(12)
    ! print*,' >> inlet.prof ... done!'
    !
    ! allocate(yin(0:jn),uin(0:jn),vin(0:jn),tin(0:jn),roin(0:jn),pin(0:jn))
    ! open(12,file=database//'./inlet.dat')
    ! read(12,*)
    ! read(12,*)(yin(j),roin(j),uin(j),vin(j),tin(j),j=0,jn)
    ! close(12)
    ! print*,' << inlet.dat ... done.'
    !
    jn=600
    allocate(yin(0:jn),uin(0:jn),vin(0:jn),tin(0:jn),roin(0:jn),pin(0:jn))
    open(12,file='profile.bl')
    read(12,*)
    read(12,*)
    read(12,*)thick1,thick2,thick3,utaw
    read(12,*)
    do j=0,jn
      read(12,*)yin(j),uin(j),vin(j),tin(j)
    end do
    close(12)
    print*,' >> profile.bl'
    roin=1.d0/tin
    !
    ! dims=h5_getdimensio('y',database//'datin/grid.h5')
    ! jn=dims(2)-1
    ! print*,' ** dimension of the input data: ',jn
    ! !
    ! print*,' ** reading mean profile data... '
    ! allocate(yin(0:jn),uin(0:jn),vin(0:jn),tin(0:jn),roin(0:jn),pin(0:jn))
    ! call H5ReadArray(roin,jn,'ro',database//'Results/meanprofile.h5')
    ! call H5ReadArray( uin,jn,'u1',database//'Results/meanprofile.h5')
    ! call H5ReadArray( vin,jn,'u2',database//'Results/meanprofile.h5')
    ! call H5ReadArray( tin,jn, 't',database//'Results/meanprofile.h5')
    ! pin=roin*tin/const2
    ! !
    ! open(12,file=database//'Results/meanprofile.dat')
    ! read(12,*)
    ! read(12,*)thick1,thick2,thick3,utaw
    ! read(12,*)
    ! do j=0,jn
    !   read(12,*)yin(j)
    ! enddo
    ! close(12)
    ! print*,' << ',database,'meanprofile.dat'
    
    gridfile='./datin/grid.2d'
    dims=h5_getdimensio('x',trim(gridfile))-1
    
    jm=dims(2)
    !
    print*,dims
    allocate(y2d(0:dims(1),0:dims(2)),y(0:jm),u(0:jm),v(0:jm),t(0:jm),p(0:jm),ro(0:jm))
    !   
    call H5ReadSubset(y2d,dims(1),dims(2),dims(3),'y',                &
                                               trim(gridfile),kslice=0)
    !
    y=y2d(0,:)
    ! do j=0,jm
    !   print*,j,y(j)
    ! enddo
    !
    ! open(12,file='gridy.dat')
    ! read(12,*)jm
    ! print*,dims
    ! allocate(y(0:jm),u(0:jm),v(0:jm),t(0:jm),p(0:jm),ro(0:jm))
    ! y(0)=0.d0
    ! do j=0,jm
    !   read(12,*)j1,y(j)
    ! end do
    ! close(12)
    ! print*,' >> gridy.dat ... done.'
    !
    ! !
    ! twall=1.d0+0.5d0*const5*dsqrt(Prandtl)
    ! twall0=t(0)
    ! print*,' ** twall=',twall
    ! !
    ! ufree=0.d0
    ! do j=0,jm
    !   !
    !   if(u(j)>ufree) then
    !     j2=j
    !     ufree=u(j)
    !   endif
    !   !
    ! enddo
    ! !
    ! tfree=t(j2)
    ! !
    ! print*,ufree,tfree,j2
    ! !
    ! do j=j2,jm
    !   !
    !   u(j)=ufree
    !   v(j)=v(j2)
    !   t(j)=tfree
    !   !
    ! enddo
    ! u=u/ufree
    ! v=v/ufree
    ! t=t/tfree
    ! !
    ! do j=0,jm
    !   t(j)=(t(j)-1.d0)/(twall0-1.d0)*(twall-1.d0)+1.d0
    !   ro(j)=1.d0/t(j)
    ! enddo
    !
    ! ro=1.d0/t
    ! p =ro*t/const2
    !
    ! yshock=3.d0
    ! xshock=yshock/tan(alfa/180.d0*pi)
    ! xshock=9.5d0
    ! yshock=xshock*tan(-alfa/180.d0*pi)
    ! print*,' ** distance from inlet to Interaction point=',xshock,yshock
    ! !
    ! yin=yin+3.05d0
    !
    do j=0,jn
        print*,j,yin(j),roin(j)
    enddo
    !
    call regularlinearinterp(yin,roin,uin,vin,tin,jn,                  &
                             y,ro,u,v,t,jm)
    !
    do j=0,jm
      if(y(j)<0.d0) then
        ro(j)=roin(0)
         u(j)= uin(0)
         v(j)= vin(0)
         t(j)= tin(0)
      endif
    enddo
    !
    ! ro=1.d0
    ! u=1.d0
    ! v=0.d0
    ! t=1.d0
    p =ro*t/const2
    !
    ! alfa=-22.76d0
    ! ! alfa=-11.4862741998204d0
    ! !
    ! roinf=1.d0
    ! uinf=1.d0
    ! vinf=0.d0
    ! tinf=1.d0

    ! pinf=roinf*tinf/const2
    ! call PostShockCal(ro(jm),u(jm),v(jm),p(jm),t(jm),        &
    !                    ro_ps, u_ps, v_ps, p_ps,t_ps,alfa,mach)
    ! print*,' ** angel of shock-generator: ',atan(v_ps/u_ps)/pi*180.d0
    ! !
    ! yshock=40.d0*tan(-alfa/180.d0*pi)
    ! print*,' ** shock jump at y=',yshock
    ! !
    ! write(*,'(2X,A)')'-------------------- before shock -----------'
    ! write(*,'(2X,A)')'       ro        u        v        p        t'
    ! write(*,'(2X,5(F9.3))')ro(jm),u(jm),v(jm),p(jm),t(jm)
    ! write(*,'(2X,A)')'--------------------  after shock -----------'
    ! write(*,'(2X,A)')'       ro        u        v        p        t'
    ! write(*,'(2X,5(F9.3))')ro_ps, u_ps, v_ps, p_ps,t_ps
    ! write(*,'(2X,A)')'---------------------------------------------'
    ! !
    call upluscal(uplus,yplus,u,y,ro,t(0),utaw=utaw)
    print*,' ** lvis=',yplus(1)/y(1)
    ! !
    ! do j=0,jm
    !   if(y(j)>=yshock) then
    !     ro(j)=ro_ps
    !     u(j)=u_ps
    !     v(j)=v_ps
    !     p(j)=p_ps
    !     t(j)=t_ps
    !   endif
    ! enddo
    !
    open(18,file='uplus.dat')
    write(18,"(2(1X,A15))")'yplus','uplus'
    write(18,"(2(1X,E15.7E3))")(yplus(j),uplus(j),j=0,jm)
    close(18)
    print*,' << uplus.dat ... done.'
    !
    open(18,file='inlet.dat')
      write(18,"(5(1X,A15))")'y','ro','u','v','T'
    do j=0,jm
      write(18,"(5(1X,E15.7E3))")y(j),ro(j),u(j),v(j),t(j)
      ! print*,j,y(j),ro(j)
    end do
    close(18)
    print*,' << inlet.dat ... done.'
    !
    open(16,file='inlet.prof')    
    write(16,"(A26)")'# parameters of inlet flow'
    write(16,"(4(1X,A14))")'delta','delta*','theter','utaw'
    write(16,"(4(1X,E14.7E2))")thick1,thick2,thick3,utaw
    write(16,"(4(1X,A14))")'ro','u1','u2','t'
    write(16,"(4(1X,E14.7E2))")(ro(j),u(j),v(j),t(j),j=0,jm)
    close(16)
    print*,' <<< inlet.prof ... done!'
    !
  end subroutine inletprofile
  !
  subroutine inletgen
    !
    use commvardefine, only: jm,Reynolds,mach,gamma
    use boundarylayer, only: tblvelocity,tblrestress
    use basicfunction, only: upluscal
    !
    integer :: j,j1,n
    real(8) :: delta,utau,rowall,twall,recovery,var1,var2,var3,error1,error2
    real(8) :: yn(0:jm),u(0:jm),v(0:jm),t(0:jm),ro(0:jm),uvd(0:jm)
    real(8),allocatable :: uplus(:),yplus(:)
    real(8),allocatable,dimension(:) :: uu,vv,ww,uv
    !
    delta=0.8298795d0
    utau=0.05115d0
    !
    open(12,file='yn.dat')
    do j=0,jm
      read(12,*)n,yn(j)
    enddo
    close(12)
    print*,' >> yn.dat'
    !
    rowall=0.4222878d0
    twall =2.427224d0
    !
    uvd=tblvelocity(yn,utau,delta,rowall,twall,Reynolds)
    call tblrestress(yn,utau,delta,rowall,twall,Reynolds,uu,vv,ww,uv)
    !
    recovery=sqrt(0.7d0)
    !
    n=0
    do while(n<15)
      !
      error1=0.d0
      error2=0.d0
      n=n+1
      !
      do j=0,jm
        var1=0.5d0*(gamma-1.d0)*mach**2*recovery
        var2=twall+((1.d0+var1)-twall)*u(j)-var1*u(j)**2
        var3=1.d0/var2
        !
        error1=max(error1,abs(t(j)-var2))
        error2=max(error2,abs(ro(j)-var3))
        t(j)=var2
        ro(j)=var3
        !
      enddo
      !
      v=0.d0
      !
      u=0.d0
      do j=0,jm
        do j1=1,j
          var1=1.d0/sqrt(0.5d0*(ro(j1)+ro(j1-1))/ro(0))
          var2=uvd(j1)-uvd(j1-1)
          u(j)=u(j)+var1*var2
        enddo
      enddo
      !
      var1=u(jm)
      u=u/var1
      !
      print*,' ** n=',n,'error temperature=',error1,'error density=',error2
      !
    enddo
    !
    call upluscal(uplus,yplus,u,yn,ro,twall,utau)
    !
    open(18,file='uplus.dat')
    write(18,"(2(1X,A15))")'yplus','uplus'
    write(18,"(2(1X,E15.7E3))")(yplus(j),uplus(j),j=0,jm)
    close(18)
    print*,' << uplus.dat ... done. '
    !
    open(18,file='generated_profile.dat')
    write(18,"(8(1X,A15))")'y','u','uvd','T','uu','vv','ww','uv'
    write(18,"(8(1X,E15.7E3))")(yn(j),u(j),uvd(j),t(j),uu(j),vv(j),ww(j),uv(j),j=0,jm)
    close(18)
    print*,' << generated_profile.dat'
    !
    open(18,file='generated_profile.dat')
    write(18,"(8(1X,A15))")'y','u','uvd','T','uu','vv','ww','uv'
    write(18,"(8(1X,E15.7E3))")(yn(j),u(j),uvd(j),t(j),uu(j),vv(j),ww(j),uv(j),j=0,jm)
    close(18)
    print*,' << generated_profile.dat'
    !
    open(16,file='inlet.prof')    
    write(16,"(A26)")'# parameters of inlet flow'
    write(16,"(4(1X,A14))")'delta','delta*','theter','utaw'
    write(16,"(4(1X,E14.7E2))")delta,0.25d0,0.d0,utau
    write(16,"(4(1X,A14))")'ro','u1','u2','t'
    write(16,"(4(1X,E14.7E2))")(ro(j),u(j),v(j),t(j),j=0,jm)
    close(16)
    print*,' <<< inlet.prof ... done!'
    !
    open(18,file='inlet.res')
    write(18,"(4(1X,A23))")'#                   R11','R22','R33','R12'
    write(18,"(4(1X,E23.15E3))")(uu(j),vv(j),ww(j),uv(j),j=0,jm)
    close(18)
    print*,' <<< inlet.res ... done.'
    !
  end subroutine inletgen
  !
  subroutine mapislice
    !
    use commvardefine, only: Reynolds,mach,gamma,const2,pi,im,jm,km
    use interpolation
    !
    real(8),allocatable,dimension(:,:,:) :: x,y,z,ro,u,v,w,p,t
    real(8) :: zmax,time,time1,time2,dx,var1,deltat,lenx,tvar,         &
               xintp,yintp,alpha
    integer :: im1,jm1,km1,im2,jm2,km2,in1,ns1,ns2,jn
    integer :: i,j,k,i1,dim(3)
    character(len=128) :: filename
    character(len=5) :: fname
    real(8),allocatable,dimension(:) :: rom,u1m,u2m,pm,tm 
    real(8),allocatable,dimension(:,:) :: rof_1,u1f_1,u2f_1,u3f_1,pf_1,tf_1,  &
                                          rof_2,u1f_2,u2f_2,u3f_2,pf_2,tf_2,  &
                                          roxy,u1xy,u2xy,u3xy,pxy,txy
    real(8),allocatable,dimension(:) :: alpha_up,ro_up,u1_up,u2_up,u3_up,p_up,tup
    !
    allocate( x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km) )
    allocate( ro(0:im,0:jm,0:km),u(0:im,0:jm,0:km),v(0:im,0:jm,0:km),  &
              w(0:im,0:jm,0:km),p(0:im,0:jm,0:km),t(0:im,0:jm,0:km) )
    !
    call H5ReadArray(x,im,jm,km,'x','datin/grid.h5')
    call H5ReadArray(y,im,jm,km,'y','datin/grid.h5')
    call H5ReadArray(z,im,jm,km,'z','datin/grid.h5')
    !
    allocate(rom(0:jm),u1m(0:jm),u2m(0:jm),pm(0:jm),tm(0:jm))
    open(12,file='datin/inlet.prof',action='read')    
    read(12,*)
    read(12,*)
    read(12,*)
    read(12,*)
    read(12,"(4(1X,E14.7E2))")(rom(j),u1m(j),u2m(j),tm(j),j=0,jm)
    close(12)
    print*,' >> datin/inlet.prof'
    !
    allocate(roxy(0:im,0:jm),u1xy(0:im,0:jm),u2xy(0:im,0:jm),          &
             u3xy(0:im,0:jm), pxy(0:im,0:jm), txy(0:im,0:jm) )
    !
    call H5ReadArray(roxy,im,jm,'ro','outdat/flowfield.h5')
    call H5ReadArray(u1xy,im,jm,'u1','outdat/flowfield.h5')
    call H5ReadArray(u2xy,im,jm,'u2','outdat/flowfield.h5')
    call H5ReadArray( pxy,im,jm, 'p','outdat/flowfield.h5')
    call H5ReadArray( txy,im,jm, 't','outdat/flowfield.h5')
    !
    allocate(rof_1(0:jm,0:km),u1f_1(0:jm,0:km),u2f_1(0:jm,0:km), &
             u3f_1(0:jm,0:km),pf_1(0:jm,0:km),tf_1(0:jm,0:km),   &
             rof_2(0:jm,0:km),u1f_2(0:jm,0:km),u2f_2(0:jm,0:km), &
             u3f_2(0:jm,0:km),pf_2(0:jm,0:km),tf_2(0:jm,0:km))
    !
    call H5ReadArray(time,'time','inflow/islice00000.h5')
    call H5ReadArray(deltat,'time','inflow/islice00001.h5')
    deltat=deltat-time
    print*,' ** time step:',deltat
    !
    do i=0,im
      !
      lenx=x(im,0,0)-x(i,0,0)
      tvar=lenx/1.d0
      !
      ns1=tvar/deltat
      ns2=ns1+1
      !
      write(fname,'(i5.5)')ns1
      call H5ReadArray(time1,'time','inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(rof_1,jm,km,'ro','inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(u1f_1,jm,km,'u1','inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(u2f_1,jm,km,'u2','inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(u3f_1,jm,km,'u3','inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray( tf_1,jm,km, 't','inflow/islice'//fname//'.h5',explicit=.false.)
      write(fname,'(i5.5)')ns2
      call H5ReadArray(time2,'time','inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(rof_2,jm,km,'ro','inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(u1f_2,jm,km,'u1','inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(u2f_2,jm,km,'u2','inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(u3f_2,jm,km,'u3','inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray( tf_2,jm,km, 't','inflow/islice'//fname//'.h5',explicit=.false.)
      !
      if(tvar>=time1 .and. tvar<=time2) then
        !
        ro(i,:,:)=linear1d(time1,time2,rof_1,rof_2,tvar)
         u(i,:,:)=linear1d(time1,time2,u1f_1,u1f_2,tvar)
         v(i,:,:)=linear1d(time1,time2,u2f_1,u2f_2,tvar)
         w(i,:,:)=linear1d(time1,time2,u3f_1,u3f_2,tvar)
         t(i,:,:)=linear1d(time1,time2, tf_1, tf_2,tvar)
         !
         if(x(i,0,0)<20.5d0 .or. x(i,0,0)>=42.d0) then
           do j=0,jm
           do k=0,km
             ro(i,j,k)=ro(i,j,k)+roxy(i,j)
              u(i,j,k)= u(i,j,k)+u1xy(i,j)
              v(i,j,k)= v(i,j,k)+u2xy(i,j)
              t(i,j,k)= t(i,j,k)+ txy(i,j)
              !
           enddo
           enddo
         else
           do j=0,jm
             ro(i,j,:)=roxy(i,j)
              u(i,j,:)=u1xy(i,j)
              v(i,j,:)=u2xy(i,j)
              t(i,j,:)= txy(i,j)
              !=
           enddo
         endif
         !
      else
        print*,tvar,'|',time1,time2,ns1,ns2
        stop
      endif
      !
      write(*,'(1A1,2(A,I0),$)')char(13),'  ** ... i:',i,'/',im
      !
    enddo
    write(*,*)''
    p=ro*t/const2
    !
    call H5WriteArray(time,'time','flowini3d.h5')
    print*,'time=',time
    call H5WriteArray(ro,im,jm,km,'ro','flowini3d.h5')
    call H5WriteArray( u,im,jm,km,'u1','flowini3d.h5')
    call H5WriteArray( v,im,jm,km,'u2','flowini3d.h5')
    call H5WriteArray( w,im,jm,km,'u3','flowini3d.h5')
    call H5WriteArray( p,im,jm,km,'p','flowini3d.h5')
    call H5WriteArray( t,im,jm,km,'t','flowini3d.h5')
    !
    call writetecbin('tecinixy.plt',x(:,:,km/2),'x',y(:,:,km/2),'y', &
                      ro(:,:,km/2),'ro',u(:,:,km/2),'u',             &
                      v(:,:,km/2),'v',w(:,:,km/2),'w',               &
                      p(:,:,km/2),'p',t(:,:,km/2),'t',im,jm)
    call writetecbin('tecinixz.plt',x(:,10,:),'x',z(:,10,:),'z',     &
                            ro(:,10,:),'ro',u(:,10,:),'u',           &
                            v(:,10,:),'v',w(:,10,:),'w',             &
                            p(:,10,:),'p',t(:,10,:),'t',im,km)
    call writetecbin('teciniyz.plt',y(im/2,:,:),'y',z(im/2,:,:),'z', &
                        ro(im/2,:,:),'ro',u(im/2,:,:),'u',           &
                        v(im/2,:,:),'v',w(im/2,:,:),'w',             &
                        p(im/2,:,:),'p',t(im/2,:,:),'t',jm,km)
    !
  end subroutine mapislice
  !
  subroutine flowinterp_swtbli
    !
    use commvardefine, only: Reynolds,mach,gamma,const2,pi
    use interpolation
    !
    real(8),allocatable,dimension(:,:,:) :: x1,y1,z1,x2,y2,z2,         &
                                            ro1,u1,v1,w1,p1,t1,        &
                                            ro2,u2,v2,w2,p2,t2
    real(8) :: zmax,time,time1,time2,dx,var1,deltat,lenx,tvar,         &
               xintp,yintp,alpha
    integer :: im1,jm1,km1,im2,jm2,km2,in1,ns1,ns2,jn,ibl
    integer :: i,j,k,i1,dim(3)
    character(len=128) :: filename
    character(len=5) :: fname
    real(8),allocatable,dimension(:) :: rom,u1m,u2m,pm,tm 
    real(8),allocatable,dimension(:,:) :: rof_1,u1f_1,u2f_1,u3f_1,pf_1,tf_1,  &
                                          rof_2,u1f_2,u2f_2,u3f_2,pf_2,tf_2
    real(8),allocatable,dimension(:) :: alpha_up,ro_up,u1_up,u2_up,u3_up,p_up,tup
    !
    dim=h5_getdimensio('x','datin/grid.h5')
    print*,dim-1
    !
    im1=dim(1)-1
    jm1=dim(2)-1
    km1=dim(3)-1
    !
    in1=0
    im1=im1+in1
    !
    allocate( x1(0:im1,0:jm1,0:km1),y1(0:im1,0:jm1,0:km1),             &
              z1(0:im1,0:jm1,0:km1),ro1(0:im1,0:jm1,0:km1),            &
              u1(0:im1,0:jm1,0:km1),v1(0:im1,0:jm1,0:km1),             &
              w1(0:im1,0:jm1,0:km1),p1(0:im1,0:jm1,0:km1),             &
              t1(0:im1,0:jm1,0:km1))
    !
    call H5ReadArray(x1(in1:im1,:,:),im1-in1,jm1,km1,'x','datin/grid.h5')
    call H5ReadArray(y1(in1:im1,:,:),im1-in1,jm1,km1,'y','datin/grid.h5')
    call H5ReadArray(z1(in1:im1,:,:),im1-in1,jm1,km1,'z','datin/grid.h5')
    !
    ! dx=x1(in1+1,0,0)-x1(in1,0,0)
    ! do i=in1-1,0,-1
    !   x1(i,:,:)=x1(i+1,:,:)-dx
    !   y1(i,:,:)=y1(i+1,:,:)
    !   z1(i,:,:)=z1(i+1,:,:)
    ! enddo
    ! print*,' ** x1(0,:,:)=',x1(0,0,0)
    ! print*,' ** x1(in1,:,:)=',x1(in1,0,0)
    !
    zmax=z1(0,0,km1)
    ! 
    print*,' ** zamx=',zmax
    !
    call H5ReadArray(time,'time','outdat/flowfield.h5')
    call H5ReadArray(ro1(in1:im1,:,:),im1-in1,jm1,km1,'ro','outdat/flowfield.h5')
    call H5ReadArray( u1(in1:im1,:,:),im1-in1,jm1,km1,'u1','outdat/flowfield.h5')
    call H5ReadArray( v1(in1:im1,:,:),im1-in1,jm1,km1,'u2','outdat/flowfield.h5')
    call H5ReadArray( w1(in1:im1,:,:),im1-in1,jm1,km1,'u3','outdat/flowfield.h5')
    call H5ReadArray( p1(in1:im1,:,:),im1-in1,jm1,km1, 'p','outdat/flowfield.h5')
    call H5ReadArray( t1(in1:im1,:,:),im1-in1,jm1,km1, 't','outdat/flowfield.h5')
    !
    x1=x1+20.d0
    !
    dim=h5_getdimensio('x','../case2/datin/grid.h5')
    print*,dim-1
    !
    im2=dim(1)-1
    jm2=dim(2)-1
    km2=dim(3)-1
    !
    allocate( x2(0:im2,0:jm2,0:km2),y2(0:im2,0:jm2,0:km2),             &
              z2(0:im2,0:jm2,0:km2),ro2(0:im2,0:jm2,0:km2),            &
              u2(0:im2,0:jm2,0:km2),v2(0:im2,0:jm2,0:km2),             &
              w2(0:im2,0:jm2,0:km2),p2(0:im2,0:jm2,0:km2),             &
              t2(0:im2,0:jm2,0:km2))
    !
    call H5ReadArray(x2,im2,jm2,km2,'x','../case2/datin/grid.h5')
    call H5ReadArray(y2,im2,jm2,km2,'y','../case2/datin/grid.h5')
    call H5ReadArray(z2,im2,jm2,km2,'z','../case2/datin/grid.h5')
    !
    allocate(rom(0:jm2),u1m(0:jm2),u2m(0:jm2),pm(0:jm2),tm(0:jm2))
    !
    open(12,file='../case2/datin/inlet.prof',action='read')    
    read(12,*)
    read(12,*)
    read(12,*)
    read(12,*)
    read(12,"(4(1X,E14.7E2))")(rom(j),u1m(j),u2m(j),tm(j),j=0,jm2)
    close(12)
    print*,' >> ../case2/datin/inlet.prof ... done!'
    !
    call regularlinearinterp(x1(:,:,0),y1(:,:,0),                      &
                             ro1,u1,v1,w1,t1,im1,jm1,                  &
                             x2(:,:,0),y2(:,:,0),                      &
                             ro2,u2,v2,w2,t2,im2,jm2,km2)
    print*,' ** interpolation done'
    !
    ibl=478
    !
    do i=0,ibl
    do j=0,jm2
      !
      ro2(i,j,:)=rom(j)
       u2(i,j,:)=u1m(j)
       v2(i,j,:)=u2m(j)
       w2(i,j,:)=0.d0
       t2(i,j,:)= tm(j)
      !
    enddo
    enddo
    !
    allocate(alpha_up(0:im1))
    !
    jn=240
    j=jn
    xintp=23.469081333894177+20.d0
    yintp=8.5212206574881471
    do i=0,im1
      var1=sqrt((x1(i,j,0)-xintp)**2+(y1(i,j,0)-yintp)**2)
      !
      alpha_up(i)=acos((x1(i,j,0)-xintp)/var1)/pi*180.d0
      !
      if(i==0 .or. i==im1) print*,i,'x',x1(i,j,0),'alpha',alpha_up(i)
      !
    enddo
    !
    do j=0,jm2
      !
      if(y2(0,j,0)>=15.57d0) then
         !
        do i=0,im2
          !
          var1=sqrt((x2(i,j,0)-xintp)**2+(y2(i,j,0)-yintp)**2)
          !
          alpha=acos((x2(i,j,0)-xintp)/var1)/pi*180.d0
          !
          ! print*,i,'x',x2(i,j,0),'alpha',alpha
          !
          do i1=1,im1
            if( alpha<=alpha_up(i1-1) .and. alpha>=alpha_up(i1) ) then
              ro2(i,j,:)=linear1d(alpha_up(i1),alpha_up(i1-1),         &
                                  ro1(i1,jn,:),ro1(i1-1,jn,:),alpha)
              u2(i,j,:)=linear1d(alpha_up(i1),alpha_up(i1-1),          &
                                  u1(i1,jn,:),u1(i1-1,jn,:),alpha)
              v2(i,j,:)=linear1d(alpha_up(i1),alpha_up(i1-1),          &
                                  v1(i1,jn,:),v1(i1-1,jn,:),alpha)
              w2(i,j,:)=linear1d(alpha_up(i1),alpha_up(i1-1),          &
                                  w1(i1,jn,:),w1(i1-1,jn,:),alpha)
              t2(i,j,:)=linear1d(alpha_up(i1),alpha_up(i1-1),          &
                                  t1(i1,jn,:),t1(i1-1,jn,:),alpha)
              exit
            elseif(alpha>alpha_up(0)) then
              ro2(i,j,:)=rom(jm2)
               u2(i,j,:)=u1m(jm2)
               v2(i,j,:)=u2m(jm2)
               w2(i,j,:)=0.d0
               t2(i,j,:)=tm(jm2)
              exit
            elseif(alpha<alpha_up(im1)) then
              ro2(i,j,:)=6.254304479988555
               u2(i,j,:)=0.75633295344776164
               v2(i,j,:)=-0.13359646473414899
               w2(i,j,:)=0.d0
               t2(i,j,:)=3.1127693769152689
            endif
          enddo
          !
        enddo
        !
      endif
      !
    enddo
    !
    ! filename='inflow/islice00000.h5'
    call H5ReadArray(time,'time','inflow/islice00000.h5')
    call H5ReadArray(deltat,'time','inflow/islice00001.h5')
    deltat=deltat-time
    print*,' ** time step:',deltat
    !
    allocate(rof_1(0:jm2,0:km2),u1f_1(0:jm2,0:km2),u2f_1(0:jm2,0:km2), &
             u3f_1(0:jm2,0:km2),pf_1(0:jm2,0:km2),tf_1(0:jm2,0:km2),   &
             rof_2(0:jm2,0:km2),u1f_2(0:jm2,0:km2),u2f_2(0:jm2,0:km2), &
             u3f_2(0:jm2,0:km2),pf_2(0:jm2,0:km2),tf_2(0:jm2,0:km2))
    !
    in1=ibl
    !
    do i=0,in1-1
      !
      lenx=x2(in1,0,0)-x2(i,0,0)
      tvar=lenx/1.d0
      !
      ns1=tvar/deltat
      ns2=ns1+1
      !
      write(fname,'(i5.5)')ns1
      call H5ReadArray(time1,'time','../case2/inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(rof_1,jm2,km2,'ro','../case2/inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(u1f_1,jm2,km2,'u1','../case2/inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(u2f_1,jm2,km2,'u2','../case2/inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(u3f_1,jm2,km2,'u3','../case2/inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray( tf_1,jm2,km2, 't','../case2/inflow/islice'//fname//'.h5',explicit=.false.)
      write(fname,'(i5.5)')ns2
      call H5ReadArray(time2,'time','../case2/inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(rof_2,jm2,km2,'ro','../case2/inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(u1f_2,jm2,km2,'u1','../case2/inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(u2f_2,jm2,km2,'u2','../case2/inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray(u3f_2,jm2,km2,'u3','../case2/inflow/islice'//fname//'.h5',explicit=.false.)
      call H5ReadArray( tf_2,jm2,km2, 't','../case2/inflow/islice'//fname//'.h5',explicit=.false.)
      !
      if(tvar>=time1 .and. tvar<=time2) then
        !
        rof_1=linear1d(time1,time2,rof_1,rof_2,tvar)
        u1f_1=linear1d(time1,time2,u1f_1,u1f_2,tvar)
        u2f_1=linear1d(time1,time2,u2f_1,u2f_2,tvar)
        u3f_1=linear1d(time1,time2,u3f_1,u3f_2,tvar)
         tf_1=linear1d(time1,time2, tf_1, tf_2,tvar)
         !
         do j=0,jm2
         do k=0,km2
           ro2(i,j,k)=ro2(i,j,k)+rof_1(j,k)
            u2(i,j,k)= u2(i,j,k)+u1f_1(j,k)
            v2(i,j,k)= v2(i,j,k)+u2f_1(j,k)
            w2(i,j,k)=           u3f_1(j,k)
            t2(i,j,k)= t2(i,j,k)+ tf_1(j,k)
            !
         enddo
         enddo
         !
      else
        print*,tvar,'|',time1,time2,ns1,ns2
        stop
      endif
      !
      write(*,'(1A1,2(A,I0),$)')char(13),'  ** ... i:',i,'/',in1-1
      !
    enddo
    write(*,*)''
    !
    p2=ro2*t2/const2
    !
    call writetecbin('tecinixy.plt',x2(:,:,km2/2),'x',y2(:,:,km2/2),'y',  &
                      ro2(:,:,km2/2),'ro',u2(:,:,km2/2),'u',           &
                      v2(:,:,km2/2),'v',w2(:,:,km2/2),'w',             &
                      p2(:,:,km2/2),'p',t2(:,:,km2/2),'t',im2,jm2)
    call writetecbin('tecinixz.plt',x2(:,10,:),'x',z2(:,10,:),'z',  &
                      ro2(:,10,:),'ro',u2(:,10,:),'u',           &
                      v2(:,10,:),'v',w2(:,10,:),'w',             &
                      p2(:,10,:),'p',t2(:,10,:),'t',im2,km2)
    call writetecbin('teciniyz.plt',y2(im2/2,:,:),'y',z2(im2/2,:,:),'z',  &
                      ro2(im2/2,:,:),'ro',u2(im2/2,:,:),'u',           &
                      v2(im2/2,:,:),'v',w2(im2/2,:,:),'w',             &
                      p2(im2/2,:,:),'p',t2(im2/2,:,:),'t',jm2,km2)
    !
    call H5WriteArray(time,'time','flowini3d.h5')
    print*,'time=',time
    call H5WriteArray(ro2,im2,jm2,km2,'ro','flowini3d.h5')
    call H5WriteArray( u2,im2,jm2,km2,'u1','flowini3d.h5')
    call H5WriteArray( v2,im2,jm2,km2,'u2','flowini3d.h5')
    call H5WriteArray( w2,im2,jm2,km2,'u3','flowini3d.h5')
    call H5WriteArray( p2,im2,jm2,km2,'p','flowini3d.h5')
    call H5WriteArray( t2,im2,jm2,km2,'t','flowini3d.h5')
    !
  end subroutine flowinterp_swtbli
  !
  subroutine flowinterp_cartesian
    !
    use interpolation
    !
    real(8),allocatable,dimension(:,:,:) :: x1,y1,z1,x2,y2,z2,         &
                                            ro1,u1,v1,w1,p1,t1,        &
                                            ro2,u2,v2,w2,p2,t2
    real(8) :: zmax,time
    integer :: im1,jm1,km1,im2,jm2,km2
    integer :: i,j,k,dim(3)
    !
    dim=h5_getdimensio('x','datin/grid.h5')
    print*,dim-1
    !
    im1=dim(1)-1
    jm1=dim(2)-1
    km1=dim(3)-1
    !
    allocate( x1(0:im1,0:jm1,0:km1),y1(0:im1,0:jm1,0:km1),             &
              z1(0:im1,0:jm1,0:km1),ro1(0:im1,0:jm1,0:km1),            &
              u1(0:im1,0:jm1,0:km1),v1(0:im1,0:jm1,0:km1),             &
              w1(0:im1,0:jm1,0:km1),p1(0:im1,0:jm1,0:km1),             &
              t1(0:im1,0:jm1,0:km1))
    !
    call H5ReadArray(x1,im1,jm1,km1,'x','datin/grid.h5')
    call H5ReadArray(y1,im1,jm1,km1,'y','datin/grid.h5')
    call H5ReadArray(z1,im1,jm1,km1,'z','datin/grid.h5')
    !
    zmax=z1(0,0,km1)
    ! 
    print*,' ** zamx=',zmax
    !
    call H5ReadArray(time,'time','outdat/flowfield.h5')
    call H5ReadArray(ro1,im1,jm1,km1,'ro','outdat/flowfield.h5')
    call H5ReadArray(u1, im1,jm1,km1,'u1','outdat/flowfield.h5')
    call H5ReadArray(v1, im1,jm1,km1,'u2','outdat/flowfield.h5')
    call H5ReadArray(w1, im1,jm1,km1,'u3','outdat/flowfield.h5')
    call H5ReadArray(p1, im1,jm1,km1, 'p','outdat/flowfield.h5')
    call H5ReadArray(t1, im1,jm1,km1, 't','outdat/flowfield.h5')
    !
    dim=h5_getdimensio('x','grid.h5')
    print*,dim-1
    !
    im2=dim(1)-1
    jm2=dim(2)-1
    km2=dim(3)-1
    !
    allocate( x2(0:im2,0:jm2,0:km2),y2(0:im2,0:jm2,0:km2),             &
              z2(0:im2,0:jm2,0:km2),ro2(0:im2,0:jm2,0:km2),            &
              u2(0:im2,0:jm2,0:km2),v2(0:im2,0:jm2,0:km2),             &
              w2(0:im2,0:jm2,0:km2),p2(0:im2,0:jm2,0:km2),             &
              t2(0:im2,0:jm2,0:km2))
    !
    call H5ReadArray(x2,im2,jm2,km2,'x','grid.h5')
    call H5ReadArray(y2,im2,jm2,km2,'y','grid.h5')
    call H5ReadArray(z2,im2,jm2,km2,'z','grid.h5')
    !
    print*,' ** zmax=',z2(0,0,km2)
    ! x2=x2+500.d0
    ! do k=0,km2
    !     !
    !     if(z2(0,0,k)>zmax) then
    !         z2(:,:,k)=z2(:,:,k)-zmax
    !     endif
    !     !
    ! enddo
    !
    ! do j=0,jm2
    ! do i=0,im2
    ! call regularlinearinterp(z1(i,j,:),ro1(i,j,:),u1(i,j,:),v1(i,j,:),  &
    !                                     w1(i,j,:),p1(i,j,:),t1(i,j,:), km1, &
    !                          z2(i,j,:),ro2(i,j,:),u2(i,j,:),v2(i,j,:),  &
    !                                     w2(i,j,:),p2(i,j,:),t2(i,j,:), km2)
    ! enddo
    ! enddo
    call regularlinearinterp(x1,y1,z1,ro1,u1,v1,w1,p1,t1,im1,jm1,km1,    &
                             x2,y2,z2,ro2,u2,v2,w2,p2,t2,im2,jm2,km2)
    !
    ! call H5ReadArray(z2,im2,jm2,km2,'z','grid2.h5')
    !
    call writetecbin('tecinixy.plt',x2(:,:,km2/2),'x',y2(:,:,km2/2),'y',  &
                      ro2(:,:,km2/2),'ro',u2(:,:,km2/2),'u',           &
                      v2(:,:,km2/2),'v',w2(:,:,km2/2),'w',             &
                      p2(:,:,km2/2),'p',t2(:,:,km2/2),'t',im2,jm2)
    call writetecbin('tecinixz.plt',x2(:,10,:),'x',z2(:,10,:),'z',  &
                      ro2(:,10,:),'ro',u2(:,10,:),'u',           &
                      v2(:,10,:),'v',w2(:,10,:),'w',             &
                      p2(:,10,:),'p',t2(:,10,:),'t',im2,km2)
    call writetecbin('teciniyz.plt',y2(im2/2,:,:),'y',z2(im2/2,:,:),'z',  &
                      ro2(im2/2,:,:),'ro',u2(im2/2,:,:),'u',           &
                      v2(im2/2,:,:),'v',w2(im2/2,:,:),'w',             &
                      p2(im2/2,:,:),'p',t2(im2/2,:,:),'t',jm2,km2)
    !
    call H5WriteArray(time,'time','flowini3d.h5')
    print*,'time=',time
    call H5WriteArray(ro2,im2,jm2,km2,'ro','flowini3d.h5')
    call H5WriteArray( u2,im2,jm2,km2,'u1','flowini3d.h5')
    call H5WriteArray( v2,im2,jm2,km2,'u2','flowini3d.h5')
    call H5WriteArray( w2,im2,jm2,km2,'u3','flowini3d.h5')
    call H5WriteArray( p2,im2,jm2,km2,'p','flowini3d.h5')
    call H5WriteArray( t2,im2,jm2,km2,'t','flowini3d.h5')
    !
  end subroutine flowinterp_cartesian
  !
  subroutine flowinterpxy(filename)
    !
    use commvardefine,only: im,jm,km
    use basicfunction,only: argtanh
    use interpolation
    !
    character(len=*),intent(in) :: filename
    !
    real(8),allocatable :: x1(:,:),y1(:,:),x2(:,:),y2(:,:)
    real(8),allocatable,dimension(:,:) :: ro1,u1,v1,w1,p1,t1,ro2,u2,v2,w2,p2,t2
    real(8),allocatable,dimension(:,:) :: xs1,ys1,xs2,ys2
    real(8) :: var1,var2,hstep
    !
    integer :: i,j,k,dim(3),im2,jm2,km2
    !
    k=0
    !
    allocate(x1(0:im,0:jm),y1(0:im,0:jm))
    call h5_read2dfrom3d(x1,im,jm,km,'x','datin/grid.2d',kslice=k)
    call h5_read2dfrom3d(y1,im,jm,km,'y','datin/grid.2d',kslice=k)
    !  
    allocate(ro1(0:im,0:jm),u1(0:im,0:jm),v1(0:im,0:jm),               &
             w1(0:im,0:jm),p1(0:im,0:jm), t1(0:im,0:jm)                )
    !
    call h5_read2dfrom3d(ro1,im,jm,km,'ro',trim(filename),kslice=k)
    call h5_read2dfrom3d(u1,im,jm,km,'u1',trim(filename),kslice=k)
    call h5_read2dfrom3d(v1,im,jm,km,'u2',trim(filename),kslice=k)
    call h5_read2dfrom3d( p1,im,jm,km, 'p',trim(filename),kslice=k)
    call h5_read2dfrom3d( t1,im,jm,km, 't',trim(filename),kslice=k)
    !
    ! call H5ReadArray(ro1,im,jm,'ro','datin/flowini2d.h5')
    ! call H5ReadArray(u1,im,jm,'u1','datin/flowini2d.h5')
    ! call H5ReadArray(v1,im,jm,'u2','datin/flowini2d.h5')
    ! call H5ReadArray(p1,im,jm, 'p','datin/flowini2d.h5')
    ! call H5ReadArray(t1,im,jm, 't','datin/flowini2d.h5')
    !
    dim=h5_getdimensio('x','grid2.h5')
    print*,dim-1
    !
    im2=dim(1)-1
    jm2=dim(2)-1
    km2=dim(3)-1
    !
    allocate(x2(0:im2,0:jm2),y2(0:im2,0:jm2))
    allocate(ro2(0:im2,0:jm2),u2(0:im2,0:jm2),v2(0:im2,0:jm2),         &
              p2(0:im2,0:jm2),t2(0:im2,0:jm2)                          )
    call h5_read2dfrom3d(x2,im2,jm2,km2,'x','grid2.h5',kslice=k)
    call h5_read2dfrom3d(y2,im2,jm2,km2,'y','grid2.h5',kslice=k)
    !
    call regularlinearinterp2d5v(x1,y1,ro1,u1,v1,p1,t1,im,jm,                   &
                                 x2,y2,ro2,u2,v2,p2,t2,im2,jm2)
    ! call InvdisInterp2D(x1,y1,ro1,u1,v1,p1,t1,im,jm,                   &
    !                     x2,y2,ro2,u2,v2,p2,t2,im2,jm2)
    !
    ! !$OMP parallel default(shared) private(i)
    ! !$OMP do
    ! do i=0,im
    !   call regularlinearinterp(y1(i,:),ro1(i,:),u1(i,:),v1(i,:),p1(i,:),t1(i,:),jm,      &
    !                            y2(i,:),ro2(i,:),u2(i,:),v2(i,:),p2(i,:),t2(i,:),jm2)
    ! enddo
    ! !$OMP end do
    ! !$OMP end  parallel
    ! !
    ! !$OMP parallel default(shared) private(i,var1,var2)
    ! !$OMP do
    ! do i=2090,im
    !   var1=y1(i,0)
    !   var2=y2(i,0)
    !   call regularlinearinterp(y1(i,:)-var1,ro1(i,:),u1(i,:),v1(i,:),p1(i,:),t1(i,:),jm,      &
    !                            y2(i,0:180)-var2,ro2(i,0:180),u2(i,0:180),v2(i,0:180),p2(i,0:180),t2(i,0:180),180)
    ! enddo
    ! !$OMP end do
    ! !$OMP end  parallel
    !
    call writetecbin('tecini.plt',x2,'x',y2,'y',ro2,'ro',            &
                                  u2,'u',v2,'v',p2,'p',t2,'t',im2,jm2)
    !
    call H5WriteArray(ro2,im2,jm2,'ro','flowini2d.h5')
    call H5WriteArray( u2,im2,jm2,'u1','flowini2d.h5')
    call H5WriteArray( v2,im2,jm2,'u2','flowini2d.h5')
    call H5WriteArray( t2,im2,jm2,'t','flowini2d.h5')
    call H5WriteArray( p2,im2,jm2,'p','flowini2d.h5')
    !
  end subroutine flowinterpxy
  !
  subroutine flowinterpxy3d
    !
    use commvardefine,only: im,jm,km,outfolder
    use basicfunction,only: argtanh
    use interpolation
    !
    real(8),allocatable :: x1(:,:),y1(:,:),x2(:,:),y2(:,:)
    real(8),allocatable,dimension(:,:,:) :: ro1,u1,v1,w1,p1,t1,ro2,u2,v2,w2,p2,t2
    real(8),allocatable,dimension(:,:) :: xs1,ys1,xs2,ys2
    real(8) :: var1,var2,hstep
    !
    integer :: i,j,k,dim(3),im2,jm2,km2
    character(len=255) :: filename
    !
    allocate(x1(0:im,0:jm),y1(0:im,0:jm))
    call h5_read2dfrom3d(x1,im,jm,km,'x','datin/grid.h5',kslice=1)
    call h5_read2dfrom3d(y1,im,jm,km,'y','datin/grid.h5',kslice=1)
    !  
    allocate(ro1(0:im,0:jm,0:km),u1(0:im,0:jm,0:km),v1(0:im,0:jm,0:km),&
             w1(0:im,0:jm,0:km), t1(0:im,0:jm,0:km) )
    !
    ! filename=outfolder//'flowfield.h5'
    filename='datin/flowini3d.h5'
    !
    call H5ReadArray(ro1,im,jm,km,'ro',trim(filename))
    call H5ReadArray(u1,im,jm,km,'u1',trim(filename))
    call H5ReadArray(v1,im,jm,km,'u2',trim(filename))
    call H5ReadArray(w1,im,jm,km,'u3',trim(filename))
    ! call H5ReadArray(p1,im,jm,km, 'p',trim(filename))
    call H5ReadArray(t1,im,jm,km, 't',trim(filename))
    !
    dim=h5_getdimensio('x','grid.h5')
    print*,dim-1
    !
    im2=dim(1)-1
    jm2=dim(2)-1
    km2=dim(3)-1
    !
    allocate(x2(0:im2,0:jm2),y2(0:im2,0:jm2))
    allocate(ro2(0:im2,0:jm2,0:km2),u2(0:im2,0:jm2,0:km2),             &
              v2(0:im2,0:jm2,0:km2),w2(0:im2,0:jm2,0:km2),             &
              t2(0:im2,0:jm2,0:km2)              )
    call h5_read2dfrom3d(x2,im2,jm2,km2,'x','grid.h5',kslice=1)
    call h5_read2dfrom3d(y2,im2,jm2,km2,'y','grid.h5',kslice=1)
    !
    ! allocate(xs1(0:im,0:jm),ys1(0:im,0:jm))
    ! xs1=x1-15.d0
    ! ys1=y1
    ! hstep=3.1696369630555110
    ! do i=0,im
    !   if(x1(i,0)>13.d0) then
    !     do j=0,jm
    !       var1=-tanh(10.d0*j/jm-6.d0)
    !       var2=0.5d0*(var1+1.d0)*(hstep-1.d0)+1.d0
    !       ! print*,j,var1,var2
    !       xs1(i,j)=(x1(i,j)-15.d0)/var2
    !       ys1(i,j)=y1(i,j)/var2
    !     enddo
    !   endif
    ! enddo
    ! !
    ! allocate(xs2(0:im2,0:jm2),ys2(0:im2,0:jm2))
    ! xs2=x2-15.d0
    ! ys2=y2
    ! hstep=2.1619295881986615
    ! do i=0,im2
    !   if(x2(i,0)>13.d0) then
    !     do j=0,jm2
    !       var1=-tanh(10.d0*j/jm2-6.d0)
    !       var2=0.5d0*(var1+1.d0)*(hstep-1.d0)+1.d0
    !       xs2(i,j)=(x2(i,j)-15.d0)/var2
    !       ys2(i,j)=y2(i,j)/var2
    !     enddo
    !   endif
    ! enddo
    !
    call InvdisInterp2D5_km(x1,y1,ro1,u1,v1,w1,t1,                   &
                            x2,y2,ro2,u2,v2,w2,t2)
    !
    call writetecbin('tecini0.plt',x2,'x',y2,'y',ro2(:,:,0),'ro',    &
                       u2(:,:,0),'u',v2(:,:,0),'v',w2(:,:,0),'w',    &
                                               t2(:,:,0),'t',im2,jm2)
    call writetecbin('tecini.plt',x2,'x',y2,'y',ro2(:,:,km/2),'ro',    &
                u2(:,:,km/2),'u',v2(:,:,km/2),'v',w2(:,:,km/2),'w',    &
                                               t2(:,:,km/2),'t',im2,jm2)
    call writetecbin('tecinim.plt',x2,'x',y2,'y',ro2(:,:,km),'ro',    &
                     u2(:,:,km),'u',v2(:,:,km),'v',w2(:,:,km),'w',    &
                                               t2(:,:,km),'t',im2,jm2)
    !
    call H5WriteArray(ro2,im2,jm2,km2,'ro','flowini3d.h5')
    call H5WriteArray( u2,im2,jm2,km2,'u1','flowini3d.h5')
    call H5WriteArray( v2,im2,jm2,km2,'u2','flowini3d.h5')
    call H5WriteArray( w2,im2,jm2,km2,'u3','flowini3d.h5')
    call H5WriteArray( t2,im2,jm2,km2,'t','flowini3d.h5')
    !
  end subroutine flowinterpxy3d
  !
  subroutine flowinterp3d
    !
    use commvardefine,only: im,jm,km,gridfile,outfolder
    use interpolation
    !
    real(8),allocatable,dimension(:,:,:) :: x2,y2,z2,ro1,u1,v1,w1,p1,t1,ro2,u2,v2,w2,p2,t2
    real(8),allocatable,dimension(:) :: x1,y1,z1
    real(8) :: var1,var2,time,p(3)
    integer :: nstep
    !
    integer :: i,j,k,k1,dim(3),im2,jm2,km2
    character(len=255) :: filename
    !
    allocate(x1(0:im),y1(0:jm),z1(0:km))
    ! call H5ReadArray(x1,im,jm,km,'x',trim(gridfile))
    ! call H5ReadArray(y1,im,jm,km,'y',trim(gridfile))
    ! call H5ReadArray(z1,im,jm,km,'z',trim(gridfile))
    call H5ReadSubset(x1,im,jm,km,'x',trim(gridfile),jslice=0,kslice=0)
    call H5ReadSubset(y1,im,jm,km,'y',trim(gridfile),islice=0,kslice=0)
    call H5ReadSubset(z1,im,jm,km,'z',trim(gridfile),islice=0,jslice=0)
    !
    allocate(ro1(0:im,0:jm,0:km),u1(0:im,0:jm,0:km),v1(0:im,0:jm,0:km),&
             w1(0:im,0:jm,0:km),p1(0:im,0:jm,0:km), t1(0:im,0:jm,0:km) )
    !
    filename=outfolder//'flowfield.h5'
    call H5ReadArray(time,'time',trim(filename))
    call H5ReadArray(nstep,'nstep',trim(filename))
    call H5ReadArray(ro1,im,jm,km,'ro',trim(filename))
    call H5ReadArray(u1,im,jm,km,'u1',trim(filename))
    call H5ReadArray(v1,im,jm,km,'u2',trim(filename))
    call H5ReadArray(w1,im,jm,km,'u3',trim(filename))
    call H5ReadArray( p1,im,jm,km, 'p',trim(filename))
    call H5ReadArray( t1,im,jm,km, 't',trim(filename))
    !
    dim=h5_getdimensio('x','grid2.h5')
    print*,dim-1
    !
    im2=dim(1)-1
    jm2=dim(2)-1
    km2=dim(3)-1
    !
    allocate(x2(0:im2,0:jm2,0:km2),y2(0:im2,0:jm2,0:km2),z2(0:im2,0:jm2,0:km2))
    allocate(ro2(0:im2,0:jm2,0:km2),u2(0:im2,0:jm2,0:km2),v2(0:im2,0:jm2,0:km2),         &
             w2(0:im2,0:jm2,0:km2),p2(0:im2,0:jm2,0:km2),t2(0:im2,0:jm2,0:km2)           )
    call H5ReadArray(x2,im2,jm2,km2,'x','grid2.h5')
    call H5ReadArray(y2,im2,jm2,km2,'y','grid2.h5')
    call H5ReadArray(z2,im2,jm2,km2,'z','grid2.h5')
    ! !
    call regularlinearinterp(x1,y1,z1,ro1,u1,v1,w1,p1,t1,        &
                             x2,y2,z2,ro2,u2,v2,w2,p2,t2)
    !
    ! call writetecbin('tecini.plt',x2,'x',y2,'y',z2,'z',ro2,'ro',       &
    !                             u2,'u',v2,'v',p2,'p',t2,'t',im2,jm2,km2)
    ! !
    ! call H5WriteArray( x2,im2,jm2,km2,'x','flowini3d.h5')
    ! call H5WriteArray( y2,im2,jm2,km2,'y','flowini3d.h5')
    ! call H5WriteArray( z2,im2,jm2,km2,'z','flowini3d.h5')
    !
    call H5WriteArray(time,'time','flowini3d.h5')
    call H5WriteArray(nstep,'nstep','flowini3d.h5')
    call H5WriteArray(ro2,im2,jm2,km2,'ro','flowini3d.h5')
    call H5WriteArray( u2,im2,jm2,km2,'u1','flowini3d.h5')
    call H5WriteArray( v2,im2,jm2,km2,'u2','flowini3d.h5')
    call H5WriteArray( w2,im2,jm2,km2,'u3','flowini3d.h5')
    call H5WriteArray( t2,im2,jm2,km2,'t','flowini3d.h5')
    call H5WriteArray( p2,im2,jm2,km2,'p','flowini3d.h5')
    !
  end subroutine flowinterp3d
  !
  subroutine flowinterpz3d
    !
    use commvardefine,only: im,jm,km,gridfile,outfolder
    use interpolation
    !
    real(8),allocatable,dimension(:,:,:) :: x2,y2,z2,ro1,u1,v1,w1,p1,t1,ro2,u2,v2,w2,p2,t2
    real(8),allocatable,dimension(:) :: x1,y1,z1,z12
    real(8) :: var1,var2,time,p(3),zmax
    integer :: nstep
    !
    integer :: i,j,k,k1,dim(3),im2,jm2,km2
    character(len=255) :: filename
    !
    allocate(x1(0:im),y1(0:jm),z1(0:km))
    ! call H5ReadArray(x1,im,jm,km,'x',trim(gridfile))
    ! call H5ReadArray(y1,im,jm,km,'y',trim(gridfile))
    ! call H5ReadArray(z1,im,jm,km,'z',trim(gridfile))
    call H5ReadSubset(x1,im,jm,km,'x',trim(gridfile),jslice=0,kslice=0)
    call H5ReadSubset(y1,im,jm,km,'y',trim(gridfile),islice=0,kslice=0)
    call H5ReadSubset(z1,im,jm,km,'z',trim(gridfile),islice=0,jslice=0)
    zmax=z1(km)
    !
    allocate(ro1(0:im,0:jm,0:km),u1(0:im,0:jm,0:km),v1(0:im,0:jm,0:km),&
             w1(0:im,0:jm,0:km),p1(0:im,0:jm,0:km), t1(0:im,0:jm,0:km) )
    !
    filename=outfolder//'flowfield.h5'
    call H5ReadArray(time,'time',trim(filename))
    call H5ReadArray(nstep,'nstep',trim(filename))
    call H5ReadArray(ro1,im,jm,km,'ro',trim(filename))
    call H5ReadArray(u1,im,jm,km,'u1',trim(filename))
    call H5ReadArray(v1,im,jm,km,'u2',trim(filename))
    call H5ReadArray(w1,im,jm,km,'u3',trim(filename))
    call H5ReadArray( p1,im,jm,km, 'p',trim(filename))
    call H5ReadArray( t1,im,jm,km, 't',trim(filename))
    !
    dim=h5_getdimensio('x','grid2.h5')
    print*,dim-1
    !
    im2=dim(1)-1
    jm2=dim(2)-1
    km2=dim(3)-1
    !
    allocate(x2(0:im2,0:jm2,0:km2),y2(0:im2,0:jm2,0:km2),              &
             z2(0:im2,0:jm2,0:km2),z12(0:km2))
    allocate(ro2(0:im2,0:jm2,0:km2),u2(0:im2,0:jm2,0:km2),v2(0:im2,0:jm2,0:km2),         &
             w2(0:im2,0:jm2,0:km2),p2(0:im2,0:jm2,0:km2),t2(0:im2,0:jm2,0:km2)           )
    call H5ReadArray(x2,im2,jm2,km2,'x','grid2.h5')
    call H5ReadArray(y2,im2,jm2,km2,'y','grid2.h5')
    call H5ReadArray(z2,im2,jm2,km2,'z','grid2.h5')
    !
    z12(:)=z2(0,0,:)
    do k=0,km2
      do while(z12(k)>=zmax)
        z12(k)=z12(k)-zmax
        print*,k,z12(k),'/',zmax
      enddo
    enddo
    !
    do j=0,jm2
    do i=0,im2
        call regularlinearinterp(z1,ro1(i,j,:),u1(i,j,:),v1(i,j,:),    &
                                     w1(i,j,:),p1(i,j,:),t1(i,j,:),km, &
                                 z12,ro2(i,j,:),u2(i,j,:),v2(i,j,:),   &
                                     w2(i,j,:),p2(i,j,:),t2(i,j,:),km2)
    enddo
      write(*,'(1A1,A20,I3,A1,$)')char(13),' ** Processing ... ',100*j/jm2,'%' 
    enddo
    !
    call writetecbin('tecinixy.plt',x2(:,:,km2/2),'x',y2(:,:,km2/2),'y',  &
                      ro2(:,:,km2/2),'ro',u2(:,:,km2/2),'u',           &
                      v2(:,:,km2/2),'v',w2(:,:,km2/2),'w',             &
                      p2(:,:,km2/2),'p',t2(:,:,km2/2),'t',im2,jm2)
    call writetecbin('tecinixz.plt',x2(:,10,:),'x',z2(:,10,:),'z',  &
                      ro2(:,10,:),'ro',u2(:,10,:),'u',              &
                      v2(:,10,:),'v',w2(:,10,:),'w',                &
                      p2(:,10,:),'p',t2(:,10,:),'t',im2,km2)
    call writetecbin('teciniyz.plt',y2(im2/2,:,:),'y',z2(im2/2,:,:),'z',  &
                      ro2(im2/2,:,:),'ro',u2(im2/2,:,:),'u',           &
                      v2(im2/2,:,:),'v',w2(im2/2,:,:),'w',             &
                      p2(im2/2,:,:),'p',t2(im2/2,:,:),'t',jm2,km2)
    !
    call H5WriteArray(time,'time','flowini3d.h5')
    call H5WriteArray(nstep,'nstep','flowini3d.h5')
    call H5WriteArray(ro2,im2,jm2,km2,'ro','flowini3d.h5')
    call H5WriteArray( u2,im2,jm2,km2,'u1','flowini3d.h5')
    call H5WriteArray( v2,im2,jm2,km2,'u2','flowini3d.h5')
    call H5WriteArray( w2,im2,jm2,km2,'u3','flowini3d.h5')
    call H5WriteArray( t2,im2,jm2,km2,'t','flowini3d.h5')
    call H5WriteArray( p2,im2,jm2,km2,'p','flowini3d.h5')
    !
  end subroutine flowinterpz3d
  !
  subroutine flowdupz(infile)
    !
    use commvardefine,only: im,jm,km
    character(len=*),intent(in) :: infile
    !
    integer :: nstep,im2,jm2,km2,i,j,k,k1
    real(8) :: time
    real(8),allocatable,dimension(:,:,:) :: ro1,u1,v1,w1,p1,t1,ro2,u2,v2,w2,p2,t2
    !
    allocate(ro1(0:im,0:jm,0:km),u1(0:im,0:jm,0:km),v1(0:im,0:jm,0:km),&
             w1(0:im,0:jm,0:km),p1(0:im,0:jm,0:km), t1(0:im,0:jm,0:km) )
    !
    call H5ReadArray(time,'time',trim(infile))
    call H5ReadArray(nstep,'nstep',trim(infile))
    call H5ReadArray(ro1,im,jm,km,'ro',trim(infile))
    call H5ReadArray(u1,im,jm,km,'u1',trim(infile))
    call H5ReadArray(v1,im,jm,km,'u2',trim(infile))
    call H5ReadArray(w1,im,jm,km,'u3',trim(infile))
    call H5ReadArray( p1,im,jm,km, 'p',trim(infile))
    call H5ReadArray( t1,im,jm,km, 't',trim(infile))
    !
    im2=im
    jm2=jm
    km2=2*km
    !
    allocate(ro2(0:im2,0:jm2,0:km2),u2(0:im2,0:jm2,0:km2),v2(0:im2,0:jm2,0:km2),&
             w2(0:im2,0:jm2,0:km2),p2(0:im2,0:jm2,0:km2), t2(0:im2,0:jm2,0:km2) )
    !
    do k=0,km2
      k1=k
      if(k1>km) k1=k1-km
      ro2(:,:,k)=ro1(:,:,k1)
       u2(:,:,k)= u1(:,:,k1)
       v2(:,:,k)= v1(:,:,k1)
       w2(:,:,k)= w1(:,:,k1)
       p2(:,:,k)= p1(:,:,k1)
       t2(:,:,k)= t1(:,:,k1)
    enddo
    !
    call H5WriteArray(time,'time','flowini3d.h5')
    call H5WriteArray(nstep,'nstep','flowini3d.h5')
    call H5WriteArray(ro2,im2,jm2,km2,'ro','flowini3d.h5')
    call H5WriteArray( u2,im2,jm2,km2,'u1','flowini3d.h5')
    call H5WriteArray( v2,im2,jm2,km2,'u2','flowini3d.h5')
    call H5WriteArray( w2,im2,jm2,km2,'u3','flowini3d.h5')
    call H5WriteArray( t2,im2,jm2,km2,'t','flowini3d.h5')
    call H5WriteArray( p2,im2,jm2,km2,'p','flowini3d.h5')
    !
  end subroutine flowdupz
  !
  subroutine datacorase
    !
    use commvardefine, only : outfolder
    !
    integer :: im,jm,km,in,jn,kn,i,j,k
    integer :: dim(3)
    real(8),allocatable,dimension(:,:,:) :: x,y,z,ro,u,v,w,t
    real(8),allocatable,dimension(:,:,:) :: x1,y1,z1,ro1,u1,v1,w1,t1
    !
    dim=h5_getdimensio('x','datin/grid.h5')
    !
    im=dim(1)-1
    jm=dim(2)-1
    km=dim(3)-1
    !
    print*,im,jm,km
    !
    allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km),    &
             ro(0:im,0:jm,0:km),u(0:im,0:jm,0:km),v(0:im,0:jm,0:km),   &
             w(0:im,0:jm,0:km),t(0:im,0:jm,0:km) )
    !
    call H5ReadArray(x,im,jm,km,'x','datin/grid.h5')
    call H5ReadArray(y,im,jm,km,'y','datin/grid.h5')
    call H5ReadArray(z,im,jm,km,'z','datin/grid.h5')
    call H5ReadArray(ro,im,jm,km,'ro',outfolder//'flowfield.h5')
    call H5ReadArray(u,im,jm,km,'u1',outfolder//'flowfield.h5')
    call H5ReadArray(v,im,jm,km,'u2',outfolder//'flowfield.h5')
    call H5ReadArray(w,im,jm,km,'u3',outfolder//'flowfield.h5')
    call H5ReadArray(t,im,jm,km, 't',outfolder//'flowfield.h5')
    !
    in=64
    jn=32
    kn=32
    allocate(x1(0:in,0:jn,0:kn),y1(0:in,0:jn,0:kn),z1(0:in,0:jn,0:kn), &
             ro1(0:in,0:jn,0:kn),u1(0:in,0:jn,0:kn),v1(0:in,0:jn,0:kn),&
             w1(0:in,0:jn,0:kn),t1(0:in,0:jn,0:kn) )
    !
    do k=0,kn
    do j=0,jn
    do i=0,in
      x1(i,j,k) =x(i*4,j*6,k*4)
      y1(i,j,k) =y(i*4,j*6,k*4)
      z1(i,j,k) =z(i*4,j*6,k*4)
      ro1(i,j,k)=ro(i*4,j*6,k*4)
      u1(i,j,k) =u(i*4,j*6,k*4)
      v1(i,j,k) =v(i*4,j*6,k*4)
      w1(i,j,k) =w(i*4,j*6,k*4)
      t1(i,j,k) =t(i*4,j*6,k*4)
    enddo
    enddo
    enddo
    !
    call H5WriteArray(x1,in,jn,kn,'x','flowini3d.h5')
    call H5WriteArray(y1,in,jn,kn,'y','flowini3d.h5')
    call H5WriteArray(z1,in,jn,kn,'z','flowini3d.h5')
    call H5WriteArray(ro1,in,jn,kn,'ro','flowini3d.h5')
    call H5WriteArray(u1,in,jn,kn,'u','flowini3d.h5')
    call H5WriteArray(v1,in,jn,kn,'v','flowini3d.h5')
    call H5WriteArray(w1,in,jn,kn,'w','flowini3d.h5')
    call H5WriteArray(t1,in,jn,kn,'t','flowini3d.h5')
    !
    call writetecbin('tecini.plt',x1,'x',y1,'y',z1,'z',ro1,'ro',       &
                              u1,'u',v1,'v',w1,'w',t1,'t',in,jn,kn)
  end subroutine datacorase
  !
  subroutine flowstat_mon_resum
    !
    integer :: nstep1,nstep2,i,ios,n
    real(8),allocatable :: var(:)
    character(len=2) :: moname
    character(len=14) :: filename
    !
    open(13,file='../flowstate.dat',position='append',action='readwrite')
    backspace(unit=13)
    read(13,*)nstep1
    print*,' ** the file is truncated at nstep: ',nstep1
    !
    allocate(var(7))
    open(14,file='flowstate.dat',action='read')
    read(14,*)
    read(14,*)nstep2
    do while(nstep2.ne.nstep1)
      read(14,*,iostat=ios)nstep2
    enddo
    !
    do while(ios==0)
      read(14,*,iostat=ios)nstep2,(var(i),i=1,7)
      write(13,"(I7,1X,E15.8E2,6(1X,E20.13E2))")nstep2,(var(i),i=1,7)
      print*, '** over write at nstep: ', nstep2
    enddo
    close(13)
    close(14)
    !
    deallocate(var)
    print*,' ** over write flowstate.dat'
    !
    allocate(var(6))
    do n=1,65
      !
      write(moname,'(i2.2)')n
      filename='monitor.'//moname//'.dat'
      !
      print*,' ** resumeing file:',filename
      !
      open(13,file='../'//filename,position='append',action='readwrite')
      backspace(unit=13)
      read(13,*)nstep1
      print*,' ** the file is truncated at nstep: ',nstep1
      !
      open(14,file=filename,action='read')
      read(14,*)
      read(14,*)
      read(14,*)nstep2
      do while(nstep2.ne.nstep1)
        read(14,*,iostat=ios)nstep2
      enddo
      !
      do while(ios==0)
        read(14,*,iostat=ios)nstep2,(var(i),i=1,6)
        write(13,"(I7,1X,E15.8E2,5(1X,E20.13E2))")nstep2,(var(i),i=1,6)
        print*, '** over write at nstep: ', nstep2
      enddo
      close(13)
      close(14)
      !
      print*,' ** overwrite',filename
      !
    enddo
    !
    deallocate(var)
    !
  end subroutine flowstat_mon_resum
  !
  subroutine flowstat_mon_proc
    !
    use interpolation
    !
    integer :: nstep,nstep1,ns,n,i,nline
    real(8),allocatable :: var(:),var1(:),var2(:)
    integer :: ios
    real(8) :: deltat
    character(len=128) :: aline
    character(len=2) :: moname
    character(len=14) :: filename
    !
    ios=0
    allocate(var(5),var1(5),var2(5))
    nstep=0
    var=0.d0
    open(13,file='flowstate.dat')
    read(13,'(A)')aline
    open(14,file='flowstate.pro')
    write(14,'(A)')aline
    nline=0
    do while(ios>=0)
      !
      nline=nline+1
      !
      read(13,*,iostat=ios)nstep1,(var1(i),i=1,5)
      !
      if (ios > 0)  then
        ! print*,'... something wrong ...',nline,ios
        ! backspace(13)
        ! read(13,'(A)')aline
        ! print*,aline
      else if (ios < 0) then
        print*,'... end of file reached ...',nstep1,ios
      else
        !
        var1(1)=nstep1*0.025d0
        !
        if(nstep1==nstep+1) then
          ! normal record
          nstep=nstep1
          deltat=var1(1)-var(1)
          var=var1
          !
          if(deltat<=0.d0) then
            print*,' deltat',deltat
            write(*,"(I7,1X,E15.8E2,4(1X,E20.13E2))")nstep,(var(i),i=1,5)
            stop
          endif
          !
          write(14,"(I7,1X,E15.8E2,4(1X,E20.13E2))")nstep,(var(i),i=1,5)
          !
        else
          !
          if(nstep1==nstep) then
            cycle 
          elseif(nstep1-nstep>1) then
            write(*,"(I7,1X,E15.8E2,4(1X,E20.13E2))")nstep,(var(i),i=1,5)
            write(*,"(I7,1X,E15.8E2,4(1X,E20.13E2))")nstep1,(var1(i),i=1,5)
            print*,' ** interpolating ... '
            !
            ! interpolation
            ns=nstep
            do while(ns<nstep1)
              ns=ns+1
              !
              var2(1)=linear1d(1.d0*nstep,1.d0*nstep1,var(1),var1(1),1.d0*ns)
              var2(2)=linear1d(1.d0*nstep,1.d0*nstep1,var(2),var1(2),1.d0*ns)
              var2(3)=linear1d(1.d0*nstep,1.d0*nstep1,var(3),var1(3),1.d0*ns)
              var2(4)=linear1d(1.d0*nstep,1.d0*nstep1,var(4),var1(4),1.d0*ns)
              var2(5)=linear1d(1.d0*nstep,1.d0*nstep1,var(5),var1(5),1.d0*ns)
              !
              write(*,"(I7,1X,E15.8E2,4(1X,E20.13E2))")ns,(var2(i),i=1,5)
              write(14,"(I7,1X,E15.8E2,4(1X,E20.13E2))")ns,(var2(i),i=1,5)
              !
            enddo
            print*,' ** interpolating ... done'
            !
            nstep=nstep1
            var=var1
            !
          else
            print*,nline,':',nstep,nstep1,var(1),var1(1)
            print*,'*********'
          endif
          ! 
        endif
        !
        
        !
      end if
      !
    enddo
    close(13)
    print*,' >> flowstate.dat'
    close(14)
    print*,' << flowstate.pro'
    !
    print*,' ** checking flowstate.pro ...'
    nstep=0
    var(1)=0.d0
    ios=0
    open(14,file='flowstate.pro')
    read(14,*)
    do while(ios>=0)
      !
      read(14,*,iostat=ios)nstep1,(var1(i),i=1,5)
      !
      if (ios < 0) then
        print*,'... end of file reached ...'
      else
        if(nstep1-nstep==1) then
          !
          if(abs(var1(1)-0.025d0*nstep1)<1.d-10) then
              nstep=nstep1
              var=var1
          else
            print*,var1(1),' time not consistent with nstep'
            stop
          endif
          !
        else
          print*,nstep,nstep1,' nstep not continues'
          stop
        endif
      endif
    enddo
    close(14)
    if(ios<0) print*,' ** flowstate.pro checked'
    !
    deallocate(var,var1,var2)
    !
    allocate(var(6),var1(6),var2(6))
    do n=1,29
      !
      write(moname,'(i2.2)')n
      filename='monitor.'//moname//'.dat'
      !
      ios=0
      nstep=0
      var=0.d0
      open(13,file=filename)
      open(14,file='monitor.'//moname//'.pro')
      read(13,'(A)')aline
      write(14,'(A)')aline
      read(13,'(A)')aline
      write(14,'(A)')aline
      nline=0
      do while(ios>=0)
        !
        nline=nline+1
        !
        read(13,*,iostat=ios)nstep1,(var1(i),i=1,6)
        !
        if (ios > 0)  then
          ! print*,'... something wrong ...',nline,ios
          ! backspace(13)
          ! read(13,'(A)')aline
          ! print*,aline
        else if (ios < 0) then
          print*,'... end of file reached ...',nstep1,ios
        else
          !
          var1(1)=nstep1*0.025d0
          !
          if(nstep1==nstep+1) then
            ! normal record
            nstep=nstep1
            deltat=var1(1)-var(1)
            var=var1
            !
            if(deltat<=0.d0) then
              print*,' deltat',deltat
              write(*,"(I7,1X,E15.8E2,5(1X,E20.13E2))")nstep,(var(i),i=1,6)
              stop
            endif
            !
            write(14,"(I7,1X,E15.8E2,5(1X,E20.13E2))")nstep,(var(i),i=1,6)
            !
          else
            !
            if(nstep1==nstep) then
              cycle 
            elseif(nstep1-nstep>1) then
              write(*,"(I7,1X,E15.8E2,5(1X,E20.13E2))")nstep,(var(i),i=1,6)
              write(*,"(I7,1X,E15.8E2,5(1X,E20.13E2))")nstep1,(var1(i),i=1,6)
              print*,' ** interpolating ... '
              !
              ! interpolation
              ns=nstep
              do while(ns<nstep1)
                ns=ns+1
                !
                var2(1)=linear1d(1.d0*nstep,1.d0*nstep1,var(1),var1(1),1.d0*ns)
                var2(2)=linear1d(1.d0*nstep,1.d0*nstep1,var(2),var1(2),1.d0*ns)
                var2(3)=linear1d(1.d0*nstep,1.d0*nstep1,var(3),var1(3),1.d0*ns)
                var2(4)=linear1d(1.d0*nstep,1.d0*nstep1,var(4),var1(4),1.d0*ns)
                var2(5)=linear1d(1.d0*nstep,1.d0*nstep1,var(5),var1(5),1.d0*ns)
                var2(6)=linear1d(1.d0*nstep,1.d0*nstep1,var(6),var1(6),1.d0*ns)
                !
                write(*,"(I7,1X,E15.8E2,5(1X,E20.13E2))")ns,(var2(i),i=1,6)
                write(14,"(I7,1X,E15.8E2,5(1X,E20.13E2))")ns,(var2(i),i=1,6)
                !
              enddo
              print*,' ** interpolating ... done'
              !
              nstep=nstep1
              var=var1
              !
            else
              print*,nline,':',nstep,nstep1,var(1),var1(1)
              print*,'*********'
            endif
            ! 
          endif
          !
        end if
        !
      enddo
      close(13)
      print*,' >> ',filename
      close(14)
      print*,' << monitor.'//moname//'.pro'
      !
      print*,' ** checking monitor.'//moname//'.pro'
      nstep=0
      var(1)=0.d0
      ios=0
      open(14,file='monitor.'//moname//'.pro')
      read(14,*)
      read(14,*)
      do while(ios>=0)
        !
        read(14,*,iostat=ios)nstep1,(var1(i),i=1,6)
        !
        if (ios < 0) then
          print*,'... end of file reached ...'
        else
          if(nstep1-nstep==1) then
            !
            if(abs(var1(1)-0.025d0*nstep1)<1.d-10) then
                nstep=nstep1
                var=var1
            else
              print*,var1(1),' time not consistent with nstep'
              stop
            endif
            !
          else
            write(*,*)nstep1,(var1(i),i=1,6)
            print*,nstep,nstep1,' nstep not continues'
            stop
          endif
        endif
      enddo
      close(14)
      if(ios<0) print*,' ** monitor.'//moname//'.pro checked.'
      !
    enddo
    !
    !
  end subroutine flowstat_mon_proc
  !
  subroutine monitordatacon(monitor_file_name)
    !
    use commvardefine, only: pinf
    !
    character(len=*),intent(in) :: monitor_file_name
    !
    integer :: n,i,nstep,nsave,nstepmax,ios,record
    character(len=4) :: mfname
    real(8) :: var(15),time,ref_time
    !
    ref_time=0.001d0/951.d0*1.d6
    !
    open(16,file=monitor_file_name,access='direct',action='read',recl=8*4)
    open(18,file=monitor_file_name//'.form',action='write')
    !
    ios=0
    record=0
    do while(ios==0)
      !
      record=record+1
      !
      read(16,rec=record,iostat=ios)nstep,time,(var(i),i=1,2)
      !
      if(ios.ne.0) exit
      !
      write(18,*)time*ref_time,var(1)/pinf,var(2)
      !
    enddo
    !
    close(16)
    print*,' >> ',monitor_file_name
    close(18)
    print*,' << ',monitor_file_name,'.form'
      !
      ! ! checking data
      ! open(16,file='monitor'//mfname//'.dat',access='direct',action='read',recl=8*4)
      ! ios=0
      ! nsave=0
      ! record=0
      ! do while(ios==0)
      !   record=record+1
      !   read(16,rec=record,iostat=ios)nstep,time,(var(i),i=1,15)
      !   !
      !   if(ios.ne.0) exit
      !   !
      !   if(nstep<10 .or. nstep>nstepmax-10) then
      !     print*,record,'|',nstep,time,var(4),ios
      !   endif
      !   !
      !   if(nstep==nsave) then
      !     nsave=nsave+1
      !   else
      !     print*,nstep,time,'-',nsave
      !     stop ' !! data not continues in monitor'//mfname//'.dat' 
      !   endif
      !   !
      ! enddo
      ! close(16)
      ! print*,' monitor.',mfname,'.dat .. checked'
      !
    !
  end subroutine monitordatacon
  !
  !+-------------------------------------------------------------------+
  !| subroutine monitorprocess.                                        |
  !+-------------------------------------------------------------------+
  subroutine monitorprocess(nmon,nss)
    !
    integer,intent(in) :: nmon,nss
    ! nmon: number of mointor files
    ! nss: step of the place statisticall stationary
    integer :: nwindow,nstepend,n,ns,ns1,ns2,nr,nline,nsl,nsline,stat, &
               nmonitor,nsta,nend,i,ios,record
    character(len=4) :: mfname
    real(8),allocatable,dimension(:) :: u1mon,u2mon,u3mon,pmon,tmon,  &
                                        romon,timemon,omegaz
    real(8),allocatable,dimension(:) :: u1f,u2f,u3f,pf,sf,spec
    integer,allocatable,dimension(:) :: nstep
    real(8),allocatable,dimension(:) :: u1ave,u2ave,u3ave,pave,tave,     &
                                        u11ave,u12ave,u22ave,u33ave,     &
                                        ppave,ttave
    real(8) :: var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11

    real(8) :: x0,lref
    real(8),allocatable :: datatemp(:)
    ! 
    x0=0.d0   !10.d0
    lref=1.d0 !0.08105193
    !
    if(nmon==0) then
      ! process all monitor data
      open(12,file='datin/monitor.dat')
      read(12,'()')
      ios=0
      nmonitor=0
      do while(ios==0)
        read(12,*,iostat=ios)n
        !
        if(ios==0) nmonitor=nmonitor+1
        !
      enddo
      close(12)
      print*,nmonitor
      print*,' ** nmonitor',nmonitor
      !
      nsta=1
      nend=nmonitor
    else
      ! process on 1 monitor file
      nsta=nmon
      nend=nmon
    endif
    !
    ! open(12,file='monitor.01.dat',position="append")
    ! backspace(12)
    ! read(12,*)nstepend
    ! close(12)
    ! print*,' ** the last step is: ',nstepend
    !
    allocate(datatemp(1:9))
    open(12,file='monitor/monitor0001.dat',access='direct',recl=8*10,action='read')
    nr=1
    do while(.true.)
      !
      read(12,rec=nr,iostat=stat)n,datatemp
      !
      if(stat==0) then
        nr = nr + 1
      else
        write(*,*) 'Something is wrong with READING operation',stat,nr,datatemp
        exit
      endif
      !
    enddo
    close(12)
    nline=nr-2
    print*,' nstep:',n,'time:',datatemp(1)
    print*,' ** the last record is: ',nline,'the last step is:',n
    !
    allocate(u1mon(nline),u2mon(nline),u3mon(nline),pmon(nline),       &
             tmon(nline),romon(nline),timemon(nline),nstep(nline),omegaz(nline))
    allocate(u1ave(nline),u2ave(nline),u3ave(nline),pave(nline),tave(nline))
    allocate(u11ave(nline),u22ave(nline),u33ave(nline),u12ave(nline),  &
             ppave(nline),ttave(nline))
    !
    print*,' ** The statisticall stationary status begins at nstep=',nss
    !
    nsline=0
    !
    do ns=nsta,nend
      !
      write(mfname,'(i4.4)') ns
      !
      open(16,file='monitor/monitor'//mfname//'.dat',access='direct',recl=8*10,action='read')
      !
      n=1
      do n=1,nline
        !
        read(16,rec=n,iostat=stat)nstep(n),datatemp
        !
        timemon(n)=datatemp(1)
        u1mon(n)  =datatemp(2)
        u2mon(n)  =datatemp(3)
        pmon(n)   =datatemp(4)
        tmon(n)   =datatemp(5)
        !
        omegaz(n) =-datatemp(7)+datatemp(8)
        !
      enddo

      ! read(16,*)
      ! read(16,*)
      ! ns=1
      ! read(16,*)nstep(ns),timemon(ns),u1mon(ns),u2mon(ns),u3mon(ns),pmon(ns),tmon(ns)
      ! do while(nstep(ns)<nstepend)
      !   ns=ns+1
      !   read(16,*)nstep(ns),timemon(ns),u1mon(ns),u2mon(ns),u3mon(ns),pmon(ns),tmon(ns)
      !   if(nsline==0 .and. nstep(ns)>=nss) then
      !     nsline=ns
      !     print*,' ** statistic start at line:',nsline
      !   endif
      ! end do
      close(16)
      print*,' >> monitor',mfname,'.dat'
      !
      open(18,file='monitor.'//mfname//'.dat')
      write(18,"(A10,6(1X,A15))")'nstep','time','u','v','p','t','omegaz'
      do n=1,nline,10
        write(18,"(I10,6(1X,E15.7E3))")nstep(n),timemon(n),u1mon(n),u2mon(n),pmon(n),tmon(n),omegaz(n)
      end do
      close(18)
      print*,' << monitor.',mfname,'.dat'
      !
      stop
      !
      u1ave=0.d0
      u2ave=0.d0
      u3ave=0.d0
       pave=0.d0
       tave=0.d0
      !
      print*,' ** Processing data for monitor.' ,mfname,'.dat'
      var1=0.d0
      var2=0.d0
      var3=0.d0
      var4=0.d0
      var5=0.d0
      !
      var6=0.d0
      var7=0.d0
      var8=0.d0
      var9=0.d0
      var10=0.d0
      var11=0.d0
      do ns1=nsline,nline
        !
        var1=var1+u1mon(ns1)
        var2=var2+u2mon(ns1)
        var3=var3+u3mon(ns1)
        var4=var4+ pmon(ns1)
        var5=var5+ tmon(ns1)
        !
        var6=var6+u1mon(ns1)*u1mon(ns1)
        var7=var7+u2mon(ns1)*u2mon(ns1)
        var8=var8+u3mon(ns1)*u3mon(ns1)
        var9=var9+ pmon(ns1)*pmon(ns1)
        var10=var10+ tmon(ns1)*tmon(ns1)
        var11=var11+ u1mon(ns1)*u2mon(ns1)
        !
        u1ave(ns1)=var1/(1.d0*(ns1-nsline+1))
        u2ave(ns1)=var2/(1.d0*(ns1-nsline+1))
        u3ave(ns1)=var3/(1.d0*(ns1-nsline+1))
         pave(ns1)=var4/(1.d0*(ns1-nsline+1))
         tave(ns1)=var5/(1.d0*(ns1-nsline+1))
        !
        u11ave(ns1)=var6/(1.d0*(ns1-nsline+1))-u1ave(ns1)*u1ave(ns1)
        u22ave(ns1)=var7/(1.d0*(ns1-nsline+1))-u2ave(ns1)*u2ave(ns1)
        u33ave(ns1)=var8/(1.d0*(ns1-nsline+1))-u3ave(ns1)*u3ave(ns1)
         ppave(ns1)=var9/(1.d0*(ns1-nsline+1))-pave(ns1)*pave(ns1)
         ttave(ns1)=var10/(1.d0*(ns1-nsline+1))-tave(ns1)*tave(ns1)
        u12ave(ns1)=var11/(1.d0*(ns1-nsline+1))-u1ave(ns1)*u2ave(ns1)
        !
        !
        write(*,'(1A1,A20,I3,A1,$)')char(13),' ** Processing ... ',100*(ns1-nsline+1)/(nline-nsline+1),'%'
      end do
      write(*,*)
      !
      open(18,file='Mon/mean.monitor.'//mfname//'.dat')
      write(18,"(A10,6(1X,A15))")'nstep','time','u1_aver','u2_aver',   &
                                            'u3_aver','p_aver','t_aver'
      do n=nsline,nline
        write(18,"(I10,6(1X,E15.7E3))")nstep(n),timemon(n),            &
                           u1ave(n),u2ave(n),u3ave(n),pave(n),tave(n)
      end do
      close(18)
      print*,' << mean.monitor.',mfname,'.dat'
      !
      open(18,file='Mon/second.monitor.'//mfname//'.dat')
      write(18,"(A10,7(1X,A15))")'nstep','time','u11','u22','u33','u12','pp','tt'
      do n=nsline,nline
        write(18,"(I10,7(1X,E15.7E3))")nstep(n),timemon(n),          &
                            u11ave(n),u22ave(n),u33ave(ns),          &
                            u12ave(n),ppave(n),ttave(n)
      end do
      close(18)
      print*,' << second.monitor.',mfname,'.dat'
      !
      ! allocate(u1f(nsline:nline),u2f(nsline:nline),u3f(nsline:nline),pf(nsline:nline),sf(nsline:nline))
      ! do ns=nsline,nline
      !   u1f(ns)=u1mon(ns)-u1ave(nline)
      !   u2f(ns)=u2mon(ns)-u2ave(nline)
      !   u3f(ns)=u3mon(ns)-u3ave(nline)
      !    pf(ns)=pmon(ns)-pave(nline)
      !    sf(ns)=sin(2.d0*pi*timemon(ns))
      ! enddo
      ! !
      ! open(18,file='sta.monitor.'//mfname//'.dat')
      ! write(18,"(6(1X,A15))")'time','u1f','u2f','u3f','pf','sf'
      ! do ns=nsline,nline
      !   write(18,"(6(1X,E15.7E3))")timemon(ns),u1f(ns),u2f(ns),u3f(ns),pf(ns),sf(ns)
      ! end do
      ! close(18)
      ! print*,' << sta.monitor.',mfname,'.dat'
      ! !
      ! call freqspetra(pf(nsline:nline),timemon(nline)-timemon(nline-1),spec)
      ! !
      ! open(18,file='spec.pressure.'//mfname//'.dat')
      ! write(18,"(2(1X,A15))")'frequency','psd'
      ! do i=1,size(spec)
      !   write(18,"(2(1X,E15.7E3))")2.d0*pi*i/(timemon(nline)-timemon(nsline))*lref, &
      !                             spec(i)/(0.5d0*lref)
      ! end do
      ! close(18)
      ! print*,' << spec.pressure.',mfname,'.dat'
      ! !
      ! deallocate(u1f,u2f,u3f,pf,sf,spec)
      !
    end do
    !
  end subroutine monitorprocess
  !+-------------------------------------------------------------------+
  !| The end of the subroutine monitorprocess.                         |
  !+-------------------------------------------------------------------+
  subroutine monitor_points
    !
    use commvardefine,only: im,jm,km,gridfile
    use h5readwrite
    !
    integer :: numon,i,j,k
    integer,allocatable :: imon(:),jmon(:),kmon(:)
    integer :: ios
    !
    real(8),allocatable :: x2d(:,:),y2d(:,:)
    !
    numon=1024
    !
    allocate(imon(numon),jmon(numon),kmon(numon))
    !
    i=0
    open(12,file='datin/monitor.dat')
    read(12,*,iostat=ios)
    do while(ios==0)
      i=i+1
      read(12,*,iostat=ios)imon(i),jmon(i),kmon(i)
      !
    enddo
    close(12)
    !
    numon=i-1
    !
    print*,' ** number of monitor points:',numon
    !
    allocate(x2d(0:im,0:jm),y2d(0:im,0:jm))
    call H5ReadSubset(x2d,im,jm,km,'x',gridfile,kslice=0)
    call H5ReadSubset(y2d,im,jm,km,'y',gridfile,kslice=0)
    !
    open(18,file='monitor_coordinates.dat')
    write(18,'(2(1X,A14))')'x','y'
    do i=1,numon
      write(18,'(2(1X,E14.7E2))')x2d(imon(i),jmon(i)),y2d(imon(i),jmon(i))
    enddo
    close(18)
    print*,' << monitor_coordinates.dat'
    !
    !
  end subroutine monitor_points
  !
  subroutine LambChaplyginDipole
    !
    use singleton
    use commvardefine,only: im,jm,km,pi
    use gradsolver,only: grad_xy
    use LinearAlgegra
    !
    integer :: i,j,k,k1,k2,n
    real(8) :: R,uu,radi,theter,ka
    real(8) :: v1(3),v2(3),v3(3)
    real(8),allocatable,dimension(:,:) :: psi,x,y

    real(8),allocatable,dimension(:,:,:) :: u
    !
    allocate( x(0:im,0:jm),y(0:im,0:jm),u(1:2,0:im,0:jm),psi(0:im,0:jm) )
    do j=0,jm
    do i=0,im
      x(i,j)=10.d0/dble(im)*i-5.d0
      y(i,j)=10.d0/dble(jm)*j-5.d0
    enddo
    enddo
    !
    R=2.d0
    uu=1.d0
    ka=3.8317/R
    do j=0,jm
    do i=0,im
      radi=sqrt( x(i,j)**2+ y(i,j)**2)
      !
      if(y(i,j)>=0.d0) then
        theter=acos(x(i,j)/radi)
      else
        theter=-acos(x(i,j)/radi)
      endif
      !
      if(radi>=R) then
        psi(i,j)=uu*(R**2/radi-radi)*(y(i,j)/radi)
      else
        psi(i,j)=-2.d0*uu*BESSEL_JN(1,ka*radi)/(ka*BESSEL_JN(0,ka*R))*(y(i,j)/radi)
      endif
      !
    enddo
    enddo
    !
    u=grad_xy(psi,x,y)
    !
    do j=0,jm
    do i=0,im
        v1(1)=u(1,i,j)
        v1(2)=u(2,i,j)
        v1(3)=0.d0
        v2(1)=0.d0
        v2(2)=0.d0
        v2(3)=-1.d0
        !
        v3=cross_product(v2,v1)
        !
        u(1,i,j)=v3(1)
        u(2,i,j)=v3(2)
    enddo
    enddo
    !
    call writetecbin('tecini.plt',x,'x',y,'y',u(1,:,:),'u',u(2,:,:),'v',psi,'psi',im,jm)
    !
  end subroutine LambChaplyginDipole
  !

  subroutine randomforce
    !
    use singleton
    use commvardefine,only: im,jm,km,x,y,z,u1,u2,u3,pi
    !
    integer :: i,j,k,k1,k2,n
    real(8) :: force(3)
    logical,save :: linit=.true.
    real(8),save :: A(3,3),B(3,3)
    real(8) :: FC,FR,var1,var2,tavg,at,xs,xe,xx,yy,zz,lwave
    real(8) :: xs1,xs2,xs3,xs4,xs5
    integer,allocatable :: seed(:)
    !
    FR=1.d0/16.d0
    FC=0.d0
    !
    var1=sqrt(2.d0/3.d0*FC/1.d0)
    var2=sqrt(1.d0/3.d0*FR/1.d0)
    !
    call random_seed(size=n)
    allocate(seed(n))
    seed = 1    ! putting arbitrary seed to all elements
    call random_seed(put=seed)
    deallocate(seed)
    !
    do j=1,3
    do i=1,3
      call random_number(A(i,j))
      call random_number(B(i,j))
    end do
    end do
    !
    A=sqrt(3.d0)*(2.d0*A-1.d0)
    B=sqrt(3.d0)*(2.d0*B-1.d0)
    !
    do j=1,3
    do i=1,3
      if(i==j) then
        A(i,j)=var1*A(i,j)
        B(i,j)=var1*B(i,j)
      else
        A(i,j)=var2*A(i,j)
        B(i,j)=var2*B(i,j)
      endif
    end do
    end do
    !
    A=1.d0*A
    B=1.d0*B
    !
    allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km), &
             u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),u3(0:im,0:jm,0:km))
    do k=0,km
    do j=0,jm
    do i=0,im
      x(i,j,k)=2.d0*pi/dble(im)*i
      y(i,j,k)=2.d0*pi/dble(jm)*j
      z(i,j,k)=2.d0*pi/dble(km)*k
    enddo
    enddo
    enddo
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      u1(i,j,k)=A(1,1)*sin(x(i,j,k))+B(1,1)*cos(x(i,j,k)) + &
                A(1,2)*sin(y(i,j,k))+B(1,2)*cos(y(i,j,k)) + &
                A(1,3)*sin(z(i,j,k))+B(1,3)*cos(z(i,j,k))
      u2(i,j,k)=A(2,1)*sin(x(i,j,k))+B(2,1)*cos(x(i,j,k)) + &
                A(2,2)*sin(y(i,j,k))+B(2,2)*cos(y(i,j,k)) + &
                A(2,3)*sin(z(i,j,k))+B(2,3)*cos(z(i,j,k))
      u3(i,j,k)=A(3,1)*sin(x(i,j,k))+B(3,1)*cos(x(i,j,k)) + &
                A(3,2)*sin(y(i,j,k))+B(3,2)*cos(y(i,j,k)) + &
                A(3,3)*sin(z(i,j,k))+B(3,3)*cos(z(i,j,k))
    enddo
    enddo
    enddo
    !
    call writetecbin('tecini.plt',x,'x',y,'y',z,'z',u1,'u',u2,'v',u3,'w',im,jm,km)
    !
  end subroutine randomforce
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to gererate the intial flow field for 
  ! the isotropic turbulence.
  ! Ref: Blaisdell, G. A., Numerical simulation of compressible 
  ! homogeneous turbulence, Phd, 1991, Stanford University
  ! Ref3:
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine HomoTurbulenceGen
    !
    use singleton
    use commvardefine,only: im,jm,km,x,y,z,u1,u2,u3,pi,Reynolds,Mach,  &
                            gridfile
    use readwrite, only: H5_ReadGrid
    use gradsolver,only: grad_3d
    !
    integer :: kmi,kmj,kmk,kmax
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! kmi: maximal wavenumber in i direction
    ! kmj: maximal wavenumber in j direction
    ! kmk: maximal wavenumber in k direction
    ! kmax: maximal wavenumber in all direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: wn1,wn2,wn3,wn12,wna
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! wn1: modul of wavenumber in i direction
    ! wn2: modul of wavenumber in j direction
    ! wn3: modul of wavenumber in k direction
    ! wn12: wn12=sqrt(wn1**2+wn2**2)
    ! wna: modul of wavenumber in all direction
    ! (k0*1.d0): the wavenumber at maximum given 
    !     spectrum
    ! Ac: the intensity of given spectrum
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, dimension(:) :: Egv
    real(8), allocatable, dimension(:,:,:) :: u1tp,u2tp,u3tp
    complex(8), allocatable, dimension(:,:,:) :: u1c,u2c,u3c
    complex(8), allocatable, dimension(:,:,:) :: u1ct,u2ct,u3ct,u4ct
    real(8), allocatable, dimension(:,:,:,:) :: du1,du2,du3
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Egv: the given initial energy spectrum
    ! u1c: the spectral velocity in k1 direction
    ! u2c: the spectral velocity in k2 direction
    ! u3c: the spectral velocity in k3 direction
    ! uct: the spectrl variable in (1~*2km)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: Kenergy,Enstropy,ITGscale,LETT,KolmLength,urms,ufmx
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Kenergy: initial volume averaged turbulent 
    !          kinetic energy
    ! Enstropy: initial volume averaged e
    !           nstrophy
    ! ITGscale: initial integral length scale
    ! LETT: initial large-eddy-turnover time
    ! KolmLength: initial Kolmogorov scale
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: k1,k2,k3,ka,k0,i,j,k
    real(8) :: ran1,ran2,ran3,rn1,rn2,rn3,var1,var2,var3,ISEA
    complex(8) :: vac1,vac2,vac3,vac4,crn1,crn2
    real(8) :: dudi,lambda,div
    !
    kmi=im/2
    kmj=jm/2
    kmk=km/2
    kmax=idnint(dsqrt((kmi**2+kmj**2+kmk**2)*1.d0))+1
    !
    allocate(Egv(0:kmax))
    allocate(u1c(-kmi:kmi,-kmj:kmj,-kmk:kmk),                        &
             u2c(-kmi:kmi,-kmj:kmj,-kmk:kmk),                        &
             u3c(-kmi:kmi,-kmj:kmj,-kmk:kmk)                         )
    allocate(u1ct(1:im,1:jm,1:km),u2ct(1:im,1:jm,1:km),              &  
             u3ct(1:im,1:jm,1:km) )
    allocate(u1tp(1:im,1:jm,1:km),u2tp(1:im,1:jm,1:km),              &
             u3tp(1:im,1:jm,1:km),u4ct(1:im,1:jm,1:km)               )
    !
    ! Give the inital energy spectrum.
    ISEA=1.d0/224.7699d0
    k0=4
    !
    do i=0,kmax
      wna=real(i,8)
      Egv(i)=IniEnergDis(ISEA,K0*1.d0,i*1.d0)
    end do
    Egv(0)=0.d0
    !
    ! Generate the random velocity field according the given energy 
    ! spectrum
    ! Blaisdell, G. A. 1991 took E(k)=Integer(ui*uicoj*dA(k)). 
    ! This program takes E(k)=Integer(0.5*ui*uicoj*dA(k)). 
    ! Therefor, we take the Ek as twice of that from Blaisdell.
    print*,' ** Generate the random velocity field according the given energy spectrum'
    do k1=0,kmi
    do k2=-kmj,kmj
    do k3=-kmk,kmk
      !
      call random_number(ran1)
      call random_number(ran2)
      call random_number(ran3)
      ! ran1,ran2,ran3: random number distributied in (0,1)
      !
      rn1=ran1*2.d0*pi
      rn2=ran2*2.d0*pi
      rn3=ran3*2.d0*pi
      !
      ! Calculate the modul of the wavenumber in each direction
      wn1=real(k1,8)
      wn2=real(k2,8)
      wn3=real(k3,8)
      wn12=dsqrt(wn1**2+wn2**2)
      wna=dsqrt(wn1**2+wn2**2+wn3**2)
      !
      ! Calculate the initial energy spectral
      if(k1==0 .and. k2==0 .and. k3==0) then
        var1=0.d0
        var2=0.d0
      else
        var1=IniEnergDis(ISEA*2.d0,K0*1.d0,wna)
        var2=dsqrt(var1/4.d0/pi/wna**2)
        ! var2=1.d0
      end if
      !
      ! Gererate the velocity spectrum in half-wavenumber space.
      crn1=rn1*(0.d0,1.d0)
      crn2=rn2*(0.d0,1.d0)
      !
      vac1=var2*cdexp(crn1)*dcos(rn3)
      vac2=var2*cdexp(crn2)*dsin(rn3)
      !
      if(k1==0 .and. k2==0 .and. k3==0) then
        u1c(k1,k2,k3)=0.d0
        u2c(k1,k2,k3)=0.d0
        u3c(k1,k2,k3)=0.d0
      elseif(k1==0 .and. k2==0) then
        u1c(k1,k2,k3)=vac1
        u2c(k1,k2,k3)=vac2
        u3c(k1,k2,k3)=0.d0
      else
        u1c(k1,k2,k3)=(vac1*wna*wn2+vac2*wn1*wn3)/(wna*wn12)
        u2c(k1,k2,k3)=(vac2*wn2*wn3-vac1*wna*wn1)/(wna*wn12)
        u3c(k1,k2,k3)=-vac2*wn12/wna
      end if
      !
    end do
    end do
    end do
    ! !
    print*,' ** Generate the velocity spectrum in another half-wavenumber space '
    ! Generate the velocity spectrum in another half-wavenumber space
    ! by using conjunction relation
    do k1=-kmi,-1
    do k2=-kmj,kmj
    do k3=-kmk,kmk
      u1c(k1,k2,k3)=dconjg(u1c(-k1,-k2,-k3))
      u2c(k1,k2,k3)=dconjg(u2c(-k1,-k2,-k3))
      u3c(k1,k2,k3)=dconjg(u3c(-k1,-k2,-k3))
    end do
    end do
    end do
    ! !
    ! Transform the spectrum from (-N/2+1,N/2) to (1,N) fo rthe 
    ! convenience of using external FFT subroutine
    print*,' ** Transform the spectrum from (-N/2+1,N/2) to (1,N)  '
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      if(i<=im/2+1) then
        k1=i-1
      else
        k1=i-im-1
      end if
      if(j<=jm/2+1) then
        k2=j-1
      else
        k2=j-jm-1
      end if
      if(k<=km/2+1) then
        k3=k-1
      else
        k3=k-km-1
      end if
      !
      u1ct(i,j,k)=u1c(k1,k2,k3)
      u2ct(i,j,k)=u2c(k1,k2,k3)
      u3ct(i,j,k)=u3c(k1,k2,k3)
    end do
    end do
    end do
    ! !
    u1ct=FFT(u1ct,inv=.true.)
    u2ct=FFT(u2ct,inv=.true.)
    u3ct=FFT(u3ct,inv=.true.)
    !
    print*,' ** project to physical space. '
    !
    allocate(u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),u3(0:im,0:jm,0:km))
    do k=1,km
    do j=1,jm
    do i=1,im
      ! multiply dsqrt(NxNyNz) for return standard FFT
      u1(i,j,k)=real(u1ct(i,j,k),8)*dsqrt(real(im*jm*km,8))
      u2(i,j,k)=real(u2ct(i,j,k),8)*dsqrt(real(im*jm*km,8))
      u3(i,j,k)=real(u3ct(i,j,k),8)*dsqrt(real(im*jm*km,8))
      !
    end do
    end do
    end do
    !
    u1(0,1:jm,1:km)=u1(im,1:jm,1:km)
    u2(0,1:jm,1:km)=u2(im,1:jm,1:km)
    u3(0,1:jm,1:km)=u3(im,1:jm,1:km)
    !
    u1(0:im,0,1:km)=u1(0:im,jm,1:km)
    u2(0:im,0,1:km)=u2(0:im,jm,1:km)
    u3(0:im,0,1:km)=u3(0:im,jm,1:km)
    !
    u1(0:im,0:jm,0)=u1(0:im,0:jm,km)
    u2(0:im,0:jm,0)=u2(0:im,0:jm,km)
    u3(0:im,0:jm,0)=u3(0:im,0:jm,km)
    ! !
    urms=0.d0
    ufmx=0.d0
    Kenergy=0.d0
    do k=1,km
    do j=1,jm
    do i=1,im
      Kenergy=Kenergy+0.5d0*(u1(i,j,k)**2+u2(i,j,k)**2+u3(i,j,k)**2)
      urms=urms+u1(i,j,k)**2+u2(i,j,k)**2+u3(i,j,k)**2
      ufmx=max(ufmx,dabs(u1(i,j,k)),dabs(u2(i,j,k)),dabs(u3(i,j,k)))
    end do
    end do
    end do
    urms=dsqrt(urms/real(im*jm*km,8))
    Kenergy=Kenergy/real(im*jm*km,8)
    !
    u1=u1/urms
    u2=u2/urms
    u3=u3/urms
    Kenergy=Kenergy/urms/urms
    urms=urms/urms
    !
    print*,'Kenergy',Kenergy,'urms',urms,'Mat=',urms*Mach
    !
    call H5WriteArray(u1,im,jm,km,'u1','flowini3d.h5')
    call H5WriteArray(u2,im,jm,km,'u2','flowini3d.h5')
    call H5WriteArray(u3,im,jm,km,'u3','flowini3d.h5')
    call H5WriteArray(k0,'k0','flowini3d.h5')
    call H5WriteArray(ISEA,'ISEA','flowini3d.h5')
    !
    gridfile='./datin/grid.h5'
    call H5_ReadGrid
    !
    allocate(du1(1:3,0:im,0:jm,0:km),du2(1:3,0:im,0:jm,0:km),          &
             du3(1:3,0:im,0:jm,0:km))
    !
    du1=grad_3d(u1,x,y,z)
    du2=grad_3d(u2)
    du3=grad_3d(u3)
    !
    dudi=0.d0
    div=0.d0
    do k=1,km
    do j=1,jm
    do i=1,im
      dudi=dudi+du1(1,i,j,k)**2+du2(2,i,j,k)**2+du3(3,i,j,k)**2
      div=div+du1(1,i,j,k)+du2(2,i,j,k)+du3(3,i,j,k)
    end do
    end do
    end do
    div=div/real(im*jm*km,8)
    print*,' ** divergence:',div
    dudi=dudi/real(im*jm*km,8)
    lambda=urms/sqrt(dudi)
    print*,'lambda',lambda,'Re_lambda',urms/sqrt(3.d0)*lambda*Reynolds
    !
    ! !u1=u1/urms
    ! !u2=u2/urms
    ! !u3=u3/urms
    ! !
    deallocate(Egv)
    deallocate(u1,u2,u3)
    deallocate(u1c,u2c,u3c)
    deallocate(u1ct,u2ct,u3ct)
    deallocate(u1tp,u2tp,u3tp)
    !
  end subroutine HomoTurbulenceGen
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine HomoTurbulenceGen.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This function is used to calcuate the spectral energy at any 
  ! wavenumber.
  ! Ref: S. JAMME, et al. Direct Numerical Simulation of the 
  ! Interaction between a Shock Wave and Various Types of Isotropic 
  ! Turbulence, Flow, Turbulence and Combustion, 2002, 68:227每268.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function IniEnergDis(Ac,k0,wnb)
    !
    real(8) :: k0,Ac,var1,wnb,IniEnergDis
    !
    var1=-2.d0*(wnb/k0)**2
    IniEnergDis=Ac*wnb**4*dexp(var1)
    !IniEnergDis=Ac*wnb**(-5.d0/3.d0)
    !
    return
    !
  end function IniEnergDis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function Ek.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
end module preprocess
!+---------------------------------------------------------------------+
!| The end of the module preprocess.                                   |
!+---------------------------------------------------------------------+