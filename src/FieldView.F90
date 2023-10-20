!+---------------------------------------------------------------------+
!| This module is for viewing flow field.                              |
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2020-04-16                                     |
!+---------------------------------------------------------------------+
module fieldview
  !
  use readwrite
  use writetec
  use WriteVTK
  !
  implicit none
  !
  Interface writeteccs
    module procedure writeteccs_one
    module procedure writeteccs_ser
  end Interface writeteccs
  !
  Interface writetecxy
    module procedure writetecxy_one
    module procedure writetecxy_ser
  end Interface writetecxy
  !
  Interface writetecxz
    module procedure writetecxz_one
  end Interface writetecxz
  !
  Interface writetecyz
    module procedure writetecyz_one
  end Interface writetecyz
  !
  contains
  !
  subroutine techannel
    use gradsolver
    use commvardefine,only: im,jm,km,mach,Reynolds,gridfile,Mach,Prandtl,gamma
    use basicfunction,only: rnsxy,upluscal,miucal
    !
    real(8),allocatable :: x2d(:,:),y2d(:,:)
    real(8),allocatable,dimension(:,:) :: roxy,u1xy,u2xy,u3xy,pxy,txy
    real(8),allocatable :: du1(:,:,:),dt(:,:,:)
    !
    integer :: i,j,k,imm
    real(8) :: ubulk,dy,miu
    !
    imm=im/2
    !
    allocate(x2d(0:im,0:jm),y2d(0:im,0:jm))
    call h5_read2dfrom3d(x2d,im,jm,km,'x',gridfile,kslice=0)
    call h5_read2dfrom3d(y2d,im,jm,km,'y',gridfile,kslice=0)
    !
    call h5_read2dfrom3d(roxy,im,jm,km,'ro','outdat/flowfield.h5',kslice=0)
    call h5_read2dfrom3d(u1xy,im,jm,km,'u1','outdat/flowfield.h5',kslice=0)
    call h5_read2dfrom3d(u2xy,im,jm,km,'u2','outdat/flowfield.h5',kslice=0)
    call h5_read2dfrom3d(u3xy,im,jm,km,'u3','outdat/flowfield.h5',kslice=0)
    call h5_read2dfrom3d( pxy,im,jm,km, 'p','outdat/flowfield.h5',kslice=0)
    call h5_read2dfrom3d( txy,im,jm,km, 't','outdat/flowfield.h5',kslice=0)
    !
    allocate(du1(2,0:im,0:jm),dt(2,0:im,0:jm))
    du1=grad_xy(u1xy,x2d,y2d)
    dt =grad_xy(txy)
    !
    ubulk=0.d0
    !
    do i=1,im
    do j=1,jm
      dy=y2d(i,j)-y2d(i,j-1)
      ubulk=ubulk+0.5d0*(u1xy(i,j)+u1xy(i,j-1))*dy
    enddo
    enddo
    !
    ubulk=ubulk/dble(im)/2.d0
    !
    print*,' ** ubulk=',ubulk
    print*,' ** Re_bulk=',ubulk*Reynolds
    !
    miu=miucal(txy(imm,0))/Reynolds
    !
    open(18,file='profile.dat')
    write(18,"(8(1X,A15))")'y','ro','u','v','p','t','dudy','dtdy'
    write(18,"(8(1X,E15.7E3))")(y2d(imm,j),roxy(imm,j),u1xy(imm,j),     &
                               u2xy(imm,j),pxy(imm,j),txy(imm,j),     &
                               du1(2,imm,j),dt(2,imm,j),j=0,jm)
    !
    write(18,*)'ubulk =',ubulk
    write(18,*)'Re_bulk=',ubulk*Reynolds
    write(18,*)'Cf     =',2.d0*du1(2,imm,0)*miu
    write(18,*)'Ch     =',2.d0*dt(2,imm,0)*miu/(Mach**2*Prandtl*(gamma-1))
    close(18)
    print*,' << profile.dat'
    !
    !
  end subroutine techannel
  !
  subroutine writeprofile(fileinput,fileoutput)
    !
    use gradsolver
    use commvardefine,only: im,jm,km,mach,Reynolds,gridfile
    use basicfunction,only: rnsxy,upluscal,miucal
    !
    character(len=*),intent(in) :: fileinput,fileoutput
    real(8),allocatable :: x(:)
    real(8),allocatable,dimension(:) :: ro,u1,u2,u3,p,t
    logical :: lfex
    integer :: i,j,k
    !
    print*,' ** extract profile: '
    !
    allocate(x(0:im))
    call H5ReadArray(x,im,'x',gridfile)
    !
    allocate(ro(0:im),u1(0:im),u2(0:im),u3(0:im),p(0:im),t(0:im))
    !
    call H5ReadArray(ro,im,'ro',trim(fileinput))
    call H5ReadArray(u1,im,'u1',trim(fileinput))
    call H5ReadArray( p,im, 'p',trim(fileinput))
    call H5ReadArray( t,im, 't',trim(fileinput))
    !
    open(18,file=trim(fileoutput))
    write(18,"(5(1X,A15))")'x','ro','u','p','t'
    write(18,"(5(1X,E15.7E3))")(x(i),ro(i),u1(i),p(i),t(i),i=0,im)
    close(18)
    print*,' << ',trim(fileoutput)
    !
    ! open(18,file='cf.dat')
    ! write(18,"(3(1X,A15))")'x','cf','pw'
    ! write(18,"(3(1X,E15.7E3))")(x2d(i,0),cf(i),pxy(i,0),i=0,im)
    ! close(18)
    ! print*,' << cf.dat ... done. '
    !
  end subroutine writeprofile
  !+-------------------------------------------------------------------+
  !| The subroutine is used to extract a xy slice from a 3D field and  |
  !| to view it with tecplot.                                          |
  !+-------------------------------------------------------------------+
  !| Writen by Jian Fang, 2020-04-17                                   |
  !+-------------------------------------------------------------------+
  subroutine writetecxy_one(fileinput,fileoutput,kslice)
    !
    use gradsolver
    use commvardefine,only: im,jm,km,mach,Reynolds,gridfile
    use basicfunction,only: rnsxy,upluscal,miucal
    !
    character(len=*),intent(in) :: fileinput,fileoutput
    integer,intent(in) :: kslice
    real(8),allocatable :: x2d(:,:),y2d(:,:),uplus(:),yplus(:)
    real(8),allocatable,dimension(:,:) :: roxy,u1xy,u2xy,u3xy,pxy,txy, &
                                          div,rns,ss,ssf,omegaz,tke,omg,miut,crit
    real(8),allocatable :: du1(:,:,:),du2(:,:,:),recount(:,:),ma(:,:),divp(:,:),gradp(:,:,:)
    real(8),allocatable :: cf(:)
    integer,allocatable :: ecount(:,:)
    real(8),allocatable :: spec1(:,:),spec2(:,:)
    logical :: lfex
    integer :: i,j,k
    real(8) :: utaw
    !
    print*,' ** extract an x-y plane at kslice=: ',kslice
    !
    allocate(x2d(0:im,0:jm),y2d(0:im,0:jm))
    call h5_read2dfrom3d(x2d,im,jm,km,'x',gridfile,kslice=kslice)
    call h5_read2dfrom3d(y2d,im,jm,km,'y',gridfile,kslice=kslice)
    !
    allocate(roxy(0:im,0:jm),u1xy(0:im,0:jm),u2xy(0:im,0:jm),          &
             u3xy(0:im,0:jm),pxy(0:im,0:jm), txy(0:im,0:jm)            )
    allocate(ecount(0:im,0:jm),recount(0:im,0:jm),ss(0:im,0:jm),       &
             omegaz(0:im,0:jm),ma(0:im,0:jm),tke(0:im,0:jm),           &
             omg(0:im,0:jm),miut(0:im,0:jm),divp(0:im,0:jm),           &
             ssf(0:im,0:jm))
    allocate(spec1(0:im,0:jm),spec2(0:im,0:jm) )
    !
    call h5_read2dfrom3d(roxy,im,jm,km,'ro',trim(fileinput),kslice=kslice)
    call h5_read2dfrom3d(u1xy,im,jm,km,'u1',trim(fileinput),kslice=kslice)
    call h5_read2dfrom3d(u2xy,im,jm,km,'u2',trim(fileinput),kslice=kslice)
    ! call h5_read2dfrom3d(u3xy,im,jm,km,'u3',trim(fileinput),kslice=kslice)
    call h5_read2dfrom3d( pxy,im,jm,km, 'p',trim(fileinput),kslice=kslice)
    call h5_read2dfrom3d( txy,im,jm,km, 't',trim(fileinput),kslice=kslice)
    !
    !
    ! call h5_read2dfrom3d(crit,im,jm,km,'crit',trim(fileinput),kslice=kslice)
    ! call h5_read2dfrom3d(spec1,im,jm,km,'sp01',trim(fileinput),kslice=kslice)
    ! call h5_read2dfrom3d(spec2,im,jm,km,'sp02',trim(fileinput),kslice=kslice)
    ! !
    ! call h5_read2dfrom3d( tke,im,jm,km, 'k',trim(fileinput),kslice=kslice)
    ! call h5_read2dfrom3d( omg,im,jm,km, 'omega',trim(fileinput),kslice=kslice)
    ! call h5_read2dfrom3d( miut,im,jm,km, 'miut',trim(fileinput),kslice=kslice)
    ! call h5_read2dfrom3d(  ss,im,jm,km, 'crit',trim(fileinput),kslice=kslice)
    call h5_read2dfrom3d( ssf,im,jm,km,'ssf',trim(fileinput),kslice=kslice)
    !vread,dim1,dim2,dim3,vname,fname
    ! call h5_read2dfrom3dint(vread=ecount,dim1=im,dim2=jm,dim3=km,     &
    !                          vname='ecount',fname=trim(fileinput),     &
    !                                                       kslice=kslice)
    ! recount=1.d0*ecount
    ! !
    allocate(du1(2,0:im,0:jm),du2(2,0:im,0:jm))
    du1=grad_xy(u1xy,x2d,y2d)
    du2=grad_xy(u2xy)
    ! ! !
    ! allocate(gradp(2,0:im,0:jm))
    ! gradp=grad_xy(pxy)
    ! ! !
    omegaz=du2(1,:,:)-du1(2,:,:)
    ! ! !
    allocate(div(0:im,0:jm))
    div=du1(1,:,:)+du2(2,:,:)
    ! ! !
    allocate(rns(0:im,0:jm))
    rns=rnsxy(roxy)
    ! ! !
    ! divp=sqrt(gradp(1,:,:)**2+gradp(2,:,:)**2)
    ! !
    ma=sqrt((u1xy**2+u2xy**2)/txy)*mach
    !
    call writetecbin(trim(fileoutput),x2d,'x',y2d,'y',roxy,'ro',       &
                     u1xy,'u',u2xy,'v',ma,'ma',div,'div',omegaz,'omegaz',ssf,'ssf',im,jm)
                     ! spec1,'Y1',spec2,'Y2',im,jm)
                                       
    !
    ! inquire(file='flowini2d.h5',exist=lfex)
    ! if(lfex) call system('mv -v flowini2d.h5 flowini2d.bak')
    ! !
    ! call H5WriteArray(roxy,im,jm,'ro','flowini2d.h5')
    ! call H5WriteArray(u1xy,im,jm,'u1','flowini2d.h5')
    ! call H5WriteArray(u2xy,im,jm,'u2','flowini2d.h5')
    ! call H5WriteArray( txy,im,jm,'t','flowini2d.h5')
    ! call H5WriteArray( pxy,im,jm,'p','flowini2d.h5')
    ! !
    ! open(16,file='inlet.prof')    
    ! write(16,"(A26)")'# parameters of inlet flow'
    ! write(16,"(4(1X,A14))")'delta','delta*','theter','utaw'
    ! write(16,"(4(1X,E14.7E2))")0.d0,0.d0,0.d0,0.d0
    ! write(16,"(4(1X,A14))")'ro','u1','u2','t'
    ! write(16,"(4(1X,E14.7E2))")(roxy(im/2,j),u1xy(im/2,j),u2xy(im/2,j),txy(im/2,j),j=0,jm)
    ! close(16)
    ! print*,' <<< inlet.prof ... done!'
    ! !
    ! open(18,file='profile.dat')
    ! write(18,*)' cf:',2.d0*du1(2,0,0)*miucal(txy(0,0))/Reynolds
    ! write(18,"(5(1X,A15))")'y','u','v','t','p'
    ! write(18,"(5(1X,E15.7E3))")(y2d(0,j),u1xy(0,j),u2xy(0,j),txy(0,j),pxy(0,j),j=0,jm)
    ! close(18)
    ! print*,' << profile.dat'
    ! !
    ! allocate(yplus(0:jm),uplus(0:jm))
    ! call upluscal(uplus=uplus,yplus=yplus,u=u1xy(0,0:jm),            &
    !               y=y2d(0,0:jm),ro=roxy(0,0:jm),tw=txy(0,0),utaw=utaw)
    ! !
    ! print*,' ** utaw=',utaw
    ! !
    ! open(18,file='uplus.dat')
    ! write(18,"(2(1X,A15))")'yplus','uplus'
    ! write(18,"(2(1X,E15.7E3))")(yplus(j),uplus(j),j=0,jm/2)
    ! close(18)
    ! print*,' << uplus.dat ... done. '
    !
    ! allocate(cf(0:im))
    ! do i=0,im
    !   cf(i)=miucal(txy(i,0))*du1(2,i,0)/Reynolds/0.5d0
    ! enddo
    ! open(18,file='cf.dat')
    ! write(18,"(3(1X,A15))")'x','cf','pw'
    ! write(18,"(3(1X,E15.7E3))")(x2d(i,0),cf(i),pxy(i,0),i=0,im)
    ! close(18)
    ! print*,' << cf.dat ... done. '
    !
  end subroutine writetecxy_one
  !
  subroutine writetecxy_ser(nsta,nend,kslice)
    !
    use chem
    use gradsolver
    use commvardefine,only: im,jm,km,mach
    use basicfunction,only: rnsxy
    !
    integer,intent(in) :: nsta,nend,kslice
    !
    real(8),allocatable :: x2d(:,:),y2d(:,:)
    real(8),allocatable,dimension(:,:) :: roxy,u1xy,u2xy,u3xy,pxy,txy, &
                                          spxy,omegaz,div,ma,rns,vtmp,qdot
    real(8),allocatable :: du1(:,:,:),du2(:,:,:)
    real(8),allocatable,dimension(:,:,:) :: spec
    character(len=4) :: fname
    character(len=3) :: spname
    character(len=255) :: filename
    integer :: n,ie,i,j
    logical :: lfex
    real(8) :: qdotmax
    !
#ifdef COMB
    call chemread('datin/Burke12.yaml')
    !
    call thermdyn

    ! call chemrep('./chemrep.dat')
#endif
    !
    print*,' ** extract an x-y plane at kslice=: ',kslice
    !
    allocate(x2d(0:im,0:jm),y2d(0:im,0:jm))
    call h5_read2dfrom3d(x2d,im,jm,km,'x',gridfile,kslice=kslice)
    call h5_read2dfrom3d(y2d,im,jm,km,'y',gridfile,kslice=kslice)
    !
    allocate(roxy(0:im,0:jm),u1xy(0:im,0:jm),u2xy(0:im,0:jm),          &
             u3xy(0:im,0:jm),pxy(0:im,0:jm), txy(0:im,0:jm),           &
             spxy(0:im,0:jm),omegaz(0:im,0:jm),div(0:im,0:jm),         &
             ma(0:im,0:jm),rns(0:im,0:jm),vtmp(0:im,0:jm) )
    !
    allocate(du1(2,0:im,0:jm),du2(2,0:im,0:jm))
    !
    allocate(spec(0:im,0:jm,1:num_species) )
    allocate(qdot(0:im,0:jm) )
    !
    do n=nsta,nend
      !
      write(fname,'(i4.4)') n
      filename=outfolder//'flowfield'//fname//'.h5'
      ! filename='outdat_ready/flowfield'//fname//'.h5'
      !
      call h5_read2dfrom3d(roxy,im,jm,km,'ro',trim(filename),kslice=kslice)
      call h5_read2dfrom3d(u1xy,im,jm,km,'u1',trim(filename),kslice=kslice)
      call h5_read2dfrom3d(u2xy,im,jm,km,'u2',trim(filename),kslice=kslice)
      call h5_read2dfrom3d(u3xy,im,jm,km,'u3',trim(filename),kslice=kslice)
      call h5_read2dfrom3d( pxy,im,jm,km, 'p',trim(filename),kslice=kslice)
      call h5_read2dfrom3d( txy,im,jm,km, 't',trim(filename),kslice=kslice)
      !
#ifdef COMB
      do jsp=1,num_species
        write(spname,'(i3.3)') jsp
        call h5_read2dfrom3d(vtmp,im,jm,km, 'sp'//spname,trim(filename),kslice=kslice)
        spec(:,:,jsp)=vtmp
      enddo
      !
      qdotmax=0.d0
      do j=0,jm
      do i=0,im
        qdot(i,j)=heatrate(roxy(i,j),txy(i,j),spec(i,j,:))
        !
        if(qdot(i,j)>qdotmax)qdotmax=qdot(i,j)
      enddo
      enddo
      !
      print*,' ** qdotmax=',qdotmax
#endif
      !
      du1=grad_xy(u1xy,x2d,y2d)
      du2=grad_xy(u2xy)
      omegaz=du2(1,:,:)-du1(2,:,:)
      ma=sqrt((u1xy**2+u2xy**2)/txy)*mach
      div=du1(1,:,:)+du2(2,:,:)
      !
      rns=rnsxy(roxy)
      !
      ie=im !-100 
      !
      call writetecbin('Results/tecxy'//fname//'.plt',x2d(0:ie,:),'x', &
                                                      y2d(0:ie,:),'y', &
                                                    roxy(0:ie,:),'ro', &
                                                     u1xy(0:ie,:),'u', &
                                                     u2xy(0:ie,:),'v', &
                                                     u3xy(0:ie,:),'w', &
                                                      pxy(0:ie,:),'p', &
                                                      txy(0:ie,:),'t', &
                                              omegaz(0:ie,:),'omegaz', &
                                                    ma(0:ie,:),'ma',   &
                                                  ! rns(0:ie,:),'rns',   &
                                                   ! qdot(0:ie,:),'qdot', &
                                                    div(0:ie,:),'div', &
                                                                   ie,jm)
      ! i=im-100
      ! open(18,file='Results/species_profile'//fname//'.dat')
      ! write(18,"(10(1X,A15))")'y','H2','H','O2','OH','O','H2O','HO2','H2O2','N2'
      ! write(18,"(10(1X,E15.7E3))")(y2d(i,j),(spec(i,j,jsp),jsp=1,num_species),j=0,jm)
      ! close(18)
      ! print*,' << Results/species_profile',fname,'.dat'
      ! !
      ! open(18,file='Results/profile'//fname//'.dat')
      ! write(18,"(6(1X,A15))")'y','ro','u','v','p','T'
      ! write(18,"(6(1X,E15.7E3))")(y2d(i,j),roxy(i,j),u1xy(i,j),u2xy(i,j),pxy(i,j),txy(i,j),j=0,jm)
      ! close(18)
      ! print*,' << Results/profile',fname,'.dat'
      !
    enddo                                 
    !
  end subroutine writetecxy_ser
  !+-------------------------------------------------------------------+
  !| The end of the subroutine writetecxy.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| The subroutine is used to extract a xz slice from a 3D field and  |
  !| to view it with tecplot.                                          |
  !+-------------------------------------------------------------------+
  !| Writen by Jian Fang, 2020-04-17                                   |
  !+-------------------------------------------------------------------+
  subroutine writetecxz_one(fileinput,fileoutput,jslice)
    !
    use commvardefine,only: im,jm,km,gridfile
    !
    character(len=*),intent(in) :: fileinput,fileoutput
    integer,intent(in) :: jslice
    !
    integer :: i,k
    real(8) :: nx,ny
    real(8),allocatable :: alfa(:)
    real(8),allocatable,dimension(:,:) :: x2d,y2d,z2d,us,un
    real(8),allocatable,dimension(:,:) :: roxz,u1xz,u2xz,u3xz,pxz,txz
    real(8),allocatable,dimension(:,:) :: rom,u1m,u2m,u3m,pm,tm
    !
    print*,' ** extract an x-z plane at jslice=: ',jslice
    !
    allocate(x2d(0:im,0:km),y2d(0:im,0:km),z2d(0:im,0:km))
    call h5_read2dfrom3d(x2d,im,jm,km,'x',gridfile,jslice=jslice)
    call h5_read2dfrom3d(y2d,im,jm,km,'y',gridfile,jslice=jslice)
    call h5_read2dfrom3d(z2d,im,jm,km,'z',gridfile,jslice=jslice)
    !
    allocate(roxz(0:im,0:km),u1xz(0:im,0:km),u2xz(0:im,0:km),          &
             u3xz(0:im,0:km), pxz(0:im,0:km), txz(0:im,0:km)           )
    !
    call h5_read2dfrom3d(roxz,im,jm,km,'ro',trim(fileinput),jslice=jslice)
    call h5_read2dfrom3d(u1xz,im,jm,km,'u1',trim(fileinput),jslice=jslice)
    call h5_read2dfrom3d(u2xz,im,jm,km,'u2',trim(fileinput),jslice=jslice)
    call h5_read2dfrom3d(u3xz,im,jm,km,'u3',trim(fileinput),jslice=jslice)
    call h5_read2dfrom3d( pxz,im,jm,km, 'p',trim(fileinput),jslice=jslice)
    call h5_read2dfrom3d( txz,im,jm,km, 't',trim(fileinput),jslice=jslice)
    !
    call writetecbin(trim(fileoutput),x2d,'x',y2d,'y',z2d,'z',roxz,'ro', &
                       u1xz,'u',u2xz,'v',u3xz,'w',pxz,'p',txz,'t',im,km)
    !
    !
  end subroutine writetecxz_one
  !+-------------------------------------------------------------------+
  !| The end of the subroutine writetecxz.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| The subroutine is used to extract a iz slice from a 3D field and  |
  !| to view it with tecplot.                                          |
  !+-------------------------------------------------------------------+
  !| Writen by Jian Fang, 2020-04-17                                   |
  !+-------------------------------------------------------------------+
  subroutine writetecyz_one(fileinput,fileoutput,islice)
    !
    use commvardefine,only: im,jm,km,gridfile
    use basicfunction,only: rnsxy,upluscal,miucal
    !
    character(len=*),intent(in) :: fileinput,fileoutput
    integer,intent(in) :: islice
    real(8),allocatable :: x2d(:,:),y2d(:,:),z2d(:,:)
    real(8),allocatable,dimension(:,:) :: royz,u1yz,u2yz,u3yz,pyz,tyz, &
                                          tke,omega,miut
    real(8),allocatable,dimension(:) :: yplus,uplus
    real(8) :: utaw
    !
    integer :: i,j
    !
    print*,' ** extract an y-z plane at islice=: ',islice
    !
    allocate(x2d(0:jm,0:km),y2d(0:jm,0:km),z2d(0:jm,0:km))
    call h5_read2dfrom3d(x2d,im,jm,km,'x',gridfile,islice=islice)
    call h5_read2dfrom3d(y2d,im,jm,km,'y',gridfile,islice=islice)
    call h5_read2dfrom3d(z2d,im,jm,km,'z',gridfile,islice=islice)
    !
    allocate(royz(0:jm,0:km),u1yz(0:jm,0:km),u2yz(0:jm,0:km),          &
             u3yz(0:jm,0:km), pyz(0:jm,0:km), tyz(0:jm,0:km),          &
             tke(0:jm,0:km),omega(0:jm,0:km),miut(0:jm,0:km)           )
    !
    call h5_read2dfrom3d(royz,im,jm,km,'ro',trim(fileinput),islice=islice)
    call h5_read2dfrom3d(u1yz,im,jm,km,'u1',trim(fileinput),islice=islice)
    call h5_read2dfrom3d(u2yz,im,jm,km,'u2',trim(fileinput),islice=islice)
    call h5_read2dfrom3d(u3yz,im,jm,km,'u3',trim(fileinput),islice=islice)
    call h5_read2dfrom3d( pyz,im,jm,km, 'p',trim(fileinput),islice=islice)
    call h5_read2dfrom3d( tyz,im,jm,km, 't',trim(fileinput),islice=islice)
    !
    call h5_read2dfrom3d(tke,im,jm,km,'k',    trim(fileinput),islice=islice)
    call h5_read2dfrom3d(omega,im,jm,km,'omega',trim(fileinput),islice=islice)
    call h5_read2dfrom3d(miut,im,jm,km,'miut', trim(fileinput),islice=islice)
    !
    call writetecbin(trim(fileoutput),x2d,'x',y2d,'y',z2d,'z',         &
                         royz,'ro',u1yz,'u',u2yz,'v',u3yz,'w',         &
                                                  pyz,'p',tyz,'t',jm,km)
    !
    open(18,file='profile.dat')
    write(18,"(9(1X,A15))")'y','ro','u','v','p','T','k','omega','miut'
    write(18,"(9(1X,E15.7E3))")(y2d(j,0),royz(j,0),u1yz(j,0),          &
                                u2yz(j,0),pyz(j,0),tyz(j,0),tke(j,0),  &
                                omega(j,0),miut(j,0),j=0,jm)
    close(18)
    print*,' << profile.dat ... done. '
    !
    allocate(yplus(0:jm/2),uplus(0:jm/2))
    !
    call upluscal(uplus=uplus,yplus=yplus,u=u1yz(0:jm/2,0),            &
                  y=y2d(0:jm/2,0),ro=royz(0:jm/2,0),tw=tyz(0,0),utaw=utaw)
    !
    print*,' ** utaw=',utaw
    !
    open(18,file='uplus.dat')
    write(18,"(2(1X,A15))")'yplus','uplus'
    write(18,"(2(1X,E15.7E3))")(yplus(j),uplus(j),j=0,jm/2)
    close(18)
    print*,' << uplus.dat ... done. '
    !
  end subroutine writetecyz_one
  !+-------------------------------------------------------------------+
  !| The end of the subroutine writetecyz.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| The subroutine is used to view 3D flow field.                     |
  !+-------------------------------------------------------------------+
  !| Writen by Jian Fang, 2020-09-10                                   |
  !+-------------------------------------------------------------------+
  subroutine writetec3d(fileinput,fileoutput)
    !
    use commvardefine,only: im,jm,km,u1,u2,u3,ro,p,t,mach
    use basicfunction,only: rns3d
    use gradsolver
    !
    character(len=*),intent(in) :: fileinput,fileoutput
    !
    real(8),allocatable :: ma(:,:,:),dro(:,:,:,:),rns(:,:,:),rns2d(:,:)
    logical :: lfex
    integer :: i,j,k
    !
    print*,' ** to view 3d flow field structures.'
    !
    call H5_ReadGrid
    !
    allocate(u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),                    &
             u3(0:im,0:jm,0:km),ro(0:im,0:jm,0:km),p(0:im,0:jm,0:km),  &
             t(0:im,0:jm,0:km),ma(0:im,0:jm,0:km))  
    !
    call H5ReadArray(ro,im,jm,km,'ro',trim(fileinput))
    call H5ReadArray(u1,im,jm,km,'u1',trim(fileinput))
    call H5ReadArray(u2,im,jm,km,'u2',trim(fileinput))
    call H5ReadArray(u3,im,jm,km,'u3',trim(fileinput))
    call H5ReadArray(p, im,jm,km,'p',trim(fileinput))
    call H5ReadArray(t, im,jm,km,'t',trim(fileinput))
    !
    ! allocate(dro(1:3,0:im,0:jm,0:km),rns(0:im,0:jm,0:km),rns2d(0:im,0:jm))
    ! !
    ! dro=grad_3d(ro,x,y,z)
    ! !
    ! rns=rns3d(ro,dro)
    ! !
    ! rns2d=0.d0
    ! !
    ! do j=0,jm
    ! do i=0,im
    !   do k=1,km
    !     rns2d(i,j)=rns2d(i,j)+rns(i,j,k)
    !   enddo
    ! enddo
    ! enddo
    ! rns2d=rns2d/dble(km)
    !
    ! inquire(file='flowini3d.h5',exist=lfex)
    ! if(lfex) call system('mv -v flowini3d.h5 flowini3d.bak')
    ! !
    ! call H5WriteArray(ro(0:im,0:jm,0:km/2),im,jm,km/2,'ro','flowini3d.h5')
    ! call H5WriteArray(u1(0:im,0:jm,0:km/2),im,jm,km/2,'u1','flowini3d.h5')
    ! call H5WriteArray(u2(0:im,0:jm,0:km/2),im,jm,km/2,'u2','flowini3d.h5')
    ! call H5WriteArray(u3(0:im,0:jm,0:km/2),im,jm,km/2,'u3','flowini3d.h5')
    ! call H5WriteArray( t(0:im,0:jm,0:km/2),im,jm,km/2,'t','flowini3d.h5')
    ! call H5WriteArray( p(0:im,0:jm,0:km/2),im,jm,km/2,'p','flowini3d.h5')
    !
    ma=sqrt((u1**2+u2**2+u3**2)/t)*mach
    ! !
    call writetecbin(trim(fileoutput),x,'x',y,'y',z,'z',ro,'ro',u1,'u',&
                         u2,'v',u3,'w',p,'p',t,'t',ma,'Ma',im,jm,km)
    ! call writetecbin(trim(fileoutput),x(:,:,km/2),'x',y(:,:,km/2),'y', &
    !                  ro(:,:,km/2),'ro',u1(:,:,km/2),'u',u2(:,:,km/2),'v', &
    !                  u3(:,:,km/2),'w',p(:,:,km/2),'p',t(:,:,km/2),'t',    &
    !                  rns(:,:,km/2),'rns',rns2d,'rns_avg',im,jm)
    !
  end subroutine writetec3d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine writetecyz.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| The subroutine is used to view 3D turbulent coherent structures.  |
  !+-------------------------------------------------------------------+
  !| Writen by Jian Fang, 2020-04-17                                   |
  !+-------------------------------------------------------------------+
  subroutine writeteccs_one(fileinput,fileoutput)
    !
    use commvardefine,only: im,jm,km,u1,u2,u3,x,y,z,ro,p,t,mach,gridfile,gamma
    use gradsolver
    use basicfunction,only:lambdaci,omega,rns3d,qcrition
    !
    character(len=*),intent(in) :: fileinput,fileoutput
    !
    ! local data
    real(8),allocatable,dimension(:,:,:,:) :: du1,du2,du3,dp
    real(8),allocatable,dimension(:,:,:) :: lambda_ci,omegax,omegay,   &
                                            omegaz,omegas,us,rns,ma,gradp,qcri
    real(8) :: alfa,nx,ny
    logical :: lfex
    integer :: i,j,is,ie,js,je,n,ia,ib
    character(len=2) :: tempname
    character(len=16) :: datadir,filename,meshfile
    !
    print*,' ** to view coherent structures.'
    !
    allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
    call H5ReadArray(x,im,jm,km,'x',gridfile)
    call H5ReadArray(y,im,jm,km,'y',gridfile)
    call H5ReadArray(z,im,jm,km,'z',gridfile)
    !
    allocate(u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),                    &
             u3(0:im,0:jm,0:km) )
    allocate(du1(1:3,0:im,0:jm,0:km),du2(1:3,0:im,0:jm,0:km),          &
             du3(1:3,0:im,0:jm,0:km) )
    !
    call H5ReadArray(u1,im,jm,km,'u1',trim(fileinput))
    call H5ReadArray(u2,im,jm,km,'u2',trim(fileinput))
    call H5ReadArray(u3,im,jm,km,'u3',trim(fileinput)) 
    !
    ! p=p*gamma*mach*mach
    !
    du1=grad_3d(u1,x,y,z)
    du2=grad_3d(u2)
    du3=grad_3d(u3)
    !
    ! dp=grad_3d(p)
    ! !
    ! gradp=sqrt(dp(1,:,:,:)**2+dp(2,:,:,:)**2+dp(3,:,:,:)**2)
    !
    allocate(lambda_ci(0:im,0:jm,0:km),omegax(0:im,0:jm,0:km) )
    !
    lambda_ci=lambdaci(du1,du2,du3)
    !
    ! allocate(qcri(0:im,0:jm,0:km))
    ! qcri=qcrition(du1,du2,du3)
    !
    omegax=omega(du1,du2,du3,'x')
    ! omegaz=omega(du1,du2,du3,'z')
    !
    ! ma=sqrt((u1**2+u2**2+u3**2)/t)*mach
    !
    ! rns=rns3d(ro)
    !
    ! omegay=omega(du1,du2,du3,'y')
    ! ! !
    ! do j=0,jm
    ! do i=0,im
    !   alfa=asin(y(i,j,0)/sqrt((x(i,j,0)-1.d0)**2+y(i,j,0)**2))
    !   nx=-cos(alfa)
    !   ny= sin(alfa)
    !   !
    !   us(i,j,:)=    ny*u1(i,j,:)-nx*u2(i,j,:)
    !   ! un(i,j,:)=    nx*u1(i,j,:)+ny*u2(i,j,:)
    !   omegas(i,j,:)=ny*omegax(i,j,:)-nx*omegay(i,j,:)
    !   !
    ! enddo
    ! enddo
    !
    is=0;    ie=im
    js=0;    je=jm

      datadir='./Results/'
      filename=fileoutput
      meshfile='mesh3d.h5'
      !
      inquire(file=trim(datadir)//trim(meshfile),exist=lfex)
      !
      if(.not. lfex) then
        call H5WriteArray(x(is:ie,js:je,0:km),ie-is,je-js,km,'x',trim(datadir)//trim(meshfile))
        call H5WriteArray(y(is:ie,js:je,0:km),ie-is,je-js,km,'y',trim(datadir)//trim(meshfile))
        call H5WriteArray(z(is:ie,js:je,0:km),ie-is,je-js,km,'z',trim(datadir)//trim(meshfile))
      endif
      !
      call xdmfwriter(dir=trim(datadir),filename=trim(filename),gridfile=trim(meshfile), &
                      var1=lambda_ci(is:ie,js:je,0:km), var1name='lambda_ci', &
                      var2=omegax(is:ie,js:je,0:km),    var2name='omegax',     &
                      var4=u1(is:ie,js:je,0:km),      var4name='u',      &
                        im=ie-is,jm=je-js,km=km)
      !
    ! call writetecbin(trim(fileoutput), x(is:ie,js:je,:),'x',           &
    !                                    y(is:ie,js:je,:),'y',           &
    !                                    z(is:ie,js:je,:),'z',           &
    !                            lambda_ci(is:ie,js:je,:),'lambda_ci',   &
    !                                   ! t(is:ie,js:je,:),'T',   &
    !                                 ! qcri(is:ie,js:je,:),'q',           &
    !                                   u1(is:ie,js:je,:),'u',           &
    !                                   u2(is:ie,js:je,:),'v',           &
    !                                   u3(is:ie,js:je,:),'w',           &
    !                                   ! ma(is:ie,js:je,:),'ma',          &
    !                                    p(is:ie,js:je,:),'p',       &
    !                               omegax(is:ie,js:je,:),'omegax',ie-is,je-js,km)
    ! is=0; ie=im
    ! js=0; je=250
    ! ! is=550; ie=1250
    ! !
    ! do n=1,10
    !   !
    !   is=800/10*(n-1)+550
    !   ie=800/10*n+550
    !   write(tempname,"(I2.2)")n
    !   !
    !   call writetecbin('teccs'//tempname//'.plt',x(is:ie,js:je,:),'x', &
    !                                   y(is:ie,js:je,:),'y',            &
    !                                   z(is:ie,js:je,:),'z',            &
    !                            lambda_ci(is:ie,js:je,:),'lambda_ci',   &
    !                                  omegax(is:ie,js:je,:),'omegax',   &
    !                                      ma(is:ie,js:je,:),'ma',       &
    !                                   gradp(is:ie,js:je,:),'dp',       &
    !                                       p(is:ie,js:je,:),'p',        &
    !                                                      ie-is,je-js,km)
    ! enddo
    !
    !
  end subroutine writeteccs_one
  !
  subroutine writeteccs_ser(nsta,nend)
    !
    use chem
    use commvardefine,only: im,jm,km,u1,u2,u3,x,y,z,Reynolds,gridfile
    use gradsolver
    use basicfunction,only:lambdaci,omega,rns3d
    !
    integer,intent(in) :: nsta,nend
    !
    ! local data
    real(8),allocatable,dimension(:,:,:,:) :: du1,du2,du3,sp,dp
    real(8),allocatable,dimension(:,:,:) :: lambda_ci,omegax,omegay,   &
                                            omegas,us,ma,div,rns,   &
                                            gradp,vtmp,qdot
    real(8),allocatable,dimension(:,:) :: cf
    real(8),allocatable,dimension(:,:,:) :: x4,y4,z4
    real(8) :: alfa,nx,ny,qdotmax,zmax,ymax
    integer :: i,j,k,is,ie,js,je,n
    character(len=4) :: fname
    character(len=3) :: spname
    character(len=255) :: filename,datadir,meshfile
    logical :: lfex
    logical,save :: firstcal=.true.
    !
    print*,' ** to view coherent structures.'
    !
#ifdef COMB
    call chemread('datin/Burke12.yaml')
    !
    call thermdyn
#endif
    !
    allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
    call H5ReadArray(x,im,jm,km,'x',gridfile)
    call H5ReadArray(y,im,jm,km,'y',gridfile)
    call H5ReadArray(z,im,jm,km,'z',gridfile)
    !
    allocate(x4(0:im,0:jm,0:km),y4(0:im,0:jm,0:km),z4(0:im,0:jm,0:km))
    x4=x
    y4=y
    z4=z
    !
    ymax=y(0,jm,0)
    zmax=z(0,0,km)
    !
    allocate(u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),                    &
             u3(0:im,0:jm,0:km),p(0:im,0:jm,0:km),                     &
             t(0:im,0:jm,0:km),ro(0:im,0:jm,0:km)                      )
    allocate(du1(1:3,0:im,0:jm,0:km),du2(1:3,0:im,0:jm,0:km),          &
             du3(1:3,0:im,0:jm,0:km))
    allocate(dp(1:3,0:im,0:jm,0:km) )

    !
    allocate(lambda_ci(0:im,0:jm,0:km),omegax(0:im,0:jm,0:km),         &
             ma(0:im,0:jm,0:km),div(0:im,0:jm,0:km),rns(0:im,0:jm,0:km), &
             qdot(0:im,0:jm,0:km),gradp(0:im,0:jm,0:km) )
    !
    allocate(sp(0:im,0:jm,0:km,1:num_species))
    allocate(vtmp(0:im,0:jm,0:km))
    !
    do n=nsta,nend
      !
      write(fname,'(i4.4)') n
      filename=outfolder//'flowfield'//fname//'.h5'
      ! filename='outdat/flowfield'//fname//'.h5'
      !
      call H5ReadArray(u1,im,jm,km,'u1',trim(filename))
      call H5ReadArray(u2,im,jm,km,'u2',trim(filename))
      call H5ReadArray(u3,im,jm,km,'u3',trim(filename))
      call H5ReadArray(p,im,jm,km,'p',trim(filename))
      call H5ReadArray(ro,im,jm,km,'ro',trim(filename))
      call H5ReadArray(t,im,jm,km,'t',trim(filename))
      !
#ifdef COMB
      do jsp=1,num_species
        write(spname,'(i3.3)') jsp
        call H5ReadArray(vtmp,im,jm,km, 'sp'//spname,trim(filename))
        sp(:,:,:,jsp)=vtmp
      enddo
      !
      qdotmax=0.d0
      !
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
      !$OMP DO
      do k=0,km
      do j=0,jm
      do i=0,im
        qdot(i,j,k)=heatrate(ro(i,j,k),t(i,j,k),sp(i,j,k,:))
        !
        if(qdot(i,j,k)>qdotmax)qdotmax=qdot(i,j,k)
      enddo
      enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      print*,' ** qdotmax=',qdotmax
#endif
      !
      if(n==nsta) then
        du1=grad_3d(u1,x,y,z)
      else
        du1=grad_3d(u1)
      endif
      du2=grad_3d(u2)
      du3=grad_3d(u3)
      !
      lambda_ci=lambdaci(du1,du2,du3)
      !
      omegax=omega(du1,du2,du3,'m')
      !
      dp=grad_3d(p)
      !
      gradp=sqrt(dp(1,:,:,:)**2+dp(2,:,:,:)**2+dp(3,:,:,:)**2)
      !
      div=du1(1,:,:,:)+du2(2,:,:,:)+du3(3,:,:,:)
      ! du3=grad_3d(p)
      ! dp=sqrt(du3(1,:,:,:)**2+du3(2,:,:,:)**2+du3(3,:,:,:)**2)
      ! !
      ! ma=sqrt((u1**2)/t)*mach
      ! !
      ! rns=rns3d(ro)
      !
      is=650; ie=1250
      js=0; je=jm
      !
      ! call writetecbin('teccfwall',x(is:ie,0,:),'x',            &
      !                                 y(is:ie,0,:),'y',            &
      !                                 z(is:ie,0,:),'z',            &
      !                                 cf(is:ie,:),'cf', &
      !                                p(is:ie,0,:),'p', ie-is,km)
      ! call writetecbin('tecslice',x(is:ie,js:je,0),'x',            &
      !                                 y(is:ie,js:je,0),'y',            &
      !                                 z(is:ie,js:je,0),'z',            &
      !                               rns(is:ie,js:je,0),'rns', &
      !                                p(is:ie,js:je,0),'p',       &
      !                                ma(is:ie,js:je,0),'mach',ie-is,je-js)
      !
      datadir='./Results/'
      filename='teccs'//fname
      meshfile='mesh3d.h5'
      !
      inquire(file=trim(datadir)//trim(meshfile),exist=lfex)
      !
      if(.not. lfex) then
        call H5WriteArray(x4(is:ie,js:je,0:km),ie-is,je-js,km,'x',trim(datadir)//trim(meshfile))
        call H5WriteArray(y4(is:ie,js:je,0:km),ie-is,je-js,km,'y',trim(datadir)//trim(meshfile))
        call H5WriteArray(z4(is:ie,js:je,0:km),ie-is,je-js,km,'z',trim(datadir)//trim(meshfile))
      endif
      !
      call xdmfwriter(dir=trim(datadir),filename=trim(filename),gridfile=trim(meshfile), &
                      var1=lambda_ci(is:ie,js:je,0:km), var1name='lambda_ci', &
                      var2=omegax(is:ie,js:je,0:km),    var2name='omega',     &
                      var3=t(is:ie,js:je,0:km),         var3name='T',         &
                      var4=qdot(is:ie,js:je,0:km),      var4name='qdot',      &
                        im=ie-is,jm=je-js,km=km)
      !

      filename='tecxy'//fname
      meshfile='meshxy.h5'
      !
      inquire(file=trim(datadir)//trim(meshfile),exist=lfex)
      if(.not. lfex) then
        call H5WriteArray(x4(is:ie,js:je,0:0),ie-is,je-js,0,'x',trim(datadir)//trim(meshfile))
        call H5WriteArray(y4(is:ie,js:je,0:0),ie-is,je-js,0,'y',trim(datadir)//trim(meshfile))
        call H5WriteArray(z4(is:ie,js:je,0:0)-zmax,ie-is,je-js,0,'z',trim(datadir)//trim(meshfile))
      endif
      !
      call xdmfwriter(trim(datadir),trim(filename),trim(meshfile), &
                          u1(is:ie,js:je,0:0),    'u1', &
                           t(is:ie,js:je,0:0),     'T', &
                         div(is:ie,js:je,0:0),   'div', &
                       gradp(is:ie,js:je,0:0), 'gradp', &
                      omegax(is:ie,js:je,0:0), 'omega', &
                        qdot(is:ie,js:je,0:0),  'qdot', &
                        im=ie-is,jm=je-js,km=0)
      !
      filename='tecxz'//fname
      meshfile='meshxz.h5'
      !
      inquire(file=trim(datadir)//trim(meshfile),exist=lfex)
      if(.not. lfex) then
        call H5WriteArray(x4(is:ie,0:0,0:km),ie-is,0,km,'x',trim(datadir)//trim(meshfile))
        call H5WriteArray(y4(is:ie,0:0,0:km)-ymax,ie-is,0,km,'y',trim(datadir)//trim(meshfile))
        call H5WriteArray(z4(is:ie,0:0,0:km),ie-is,0,km,'z',trim(datadir)//trim(meshfile))
      endif
      !
      call xdmfwriter(trim(datadir),trim(filename),trim(meshfile), &
                          u1(is:ie,0:0,0:km),    'u1',    &
                           t(is:ie,0:0,0:km),     'T',    &
                         div(is:ie,0:0,0:km),   'div',    &
                       gradp(is:ie,0:0,0:km), 'gradp',    &
                      omegax(is:ie,0:0,0:km),  'omega',   &
                        qdot(is:ie,0:0,0:km),   'qdot',   &
                        im=ie-is,jm=0,km=km)
      ! call writeprvbin(trim(filename),x(is:ie,js:je,:),'x',            &
      !                                 y(is:ie,js:je,:),'y',            &
      !                                 z(is:ie,js:je,:),'z',            &
      !                            lambda_ci(is:ie,js:je,:),'lambda_ci', &
      !                            omegax(is:ie,js:je,:),'omega',        &
      !                               gradp(is:ie,js:je,:),'gradp',      &
      !                               ie-is,je-js,km)

      ! call writetecbin(trim(filename),x(is:ie,js:je,:),'x',            &
      !                                 y(is:ie,js:je,:),'y',            &
      !                                 z(is:ie,js:je,:),'z',            &
      !                            lambda_ci(is:ie,js:je,:),'lambda_ci', &
      !                            omegax(is:ie,js:je,:),'omega',        &
      !                                u1(is:ie,js:je,:),'u',            &
      !                               gradp(is:ie,js:je,:),'gradp',      &
      !                              !    t(is:ie,js:je,:),'T',            &
      !                              ! qdot(is:ie,js:je,:),'qdot',         &
      !                            !     ma(is:ie,js:je,:),'Mach',         &
      !                            !     dp(is:ie,js:je,:),'div_p',        &
      !                                ! rns(is:ie,js:je,:),'rns',         &
      !                                ! div(is:ie,js:je,:),'div',         &
      !                                                    ie-is,je-js,km)
    enddo
    !
    !
  end subroutine writeteccs_ser
  !+-------------------------------------------------------------------+
  !| The end of the subroutine writetecyz.                             |
  !+-------------------------------------------------------------------+
  !!
  subroutine ksliceview(nsta,nend)
    !
    use gradsolver
    use commvardefine,only: im,jm,km,const2
    use basicfunction,only: rnsxy
    !
    ! argument
    integer,intent(in) :: nsta,nend
    !
    ! local data
    real(8),allocatable :: x2d(:,:),y2d(:,:)
    real(8),allocatable,dimension(:,:) :: ro,u1,u2,u3,p,t,ma,divp, &
                                          div,rns,omegax,omegay,omegaz
    real(8),allocatable :: du1(:,:,:),du2(:,:,:),du3(:,:,:),dp(:,:,:)
    real(8) :: time
    integer :: n,ix,jy
    character(len=5) :: fname
    character(len=255) :: filename
    !
    allocate(x2d(0:im,0:jm),y2d(0:im,0:jm))
    call h5_read2dfrom3d(x2d,im,jm,km,'x','datin/grid.h5',kslice=0)
    call h5_read2dfrom3d(y2d,im,jm,km,'y','datin/grid.h5',kslice=0)
    !
    allocate(ro(0:im,0:jm),u1(0:im,0:jm),u2(0:im,0:jm),          &
             u3(0:im,0:jm),p(0:im,0:jm), t(0:im,0:jm)            )
    !
    allocate(du1(2,0:im,0:jm),du2(2,0:im,0:jm))
    !
    allocate(dp(1:2,0:im,0:jm),divp(0:im,0:jm))
    allocate(omegaz(0:im,0:jm))
    allocate(div(0:im,0:jm))
    allocate(rns(0:im,0:jm))
    allocate(ma(0:im,0:jm))
    !
    do n=nsta,nend
      !
      if(n>=0.and.n<=9)then
        write(fname,'(4h0000,i1)') n
      else if(n>=10.and.n<=99)then
        write(fname,'(3h000,i2)') n
      else if(n>=100.and.n<=999)then
        write(fname,'(2h00,i3)') n
      else if(n>=1000.and.n<=9999)then
        write(fname,'(1h0,i4)') n
      else if(n>=10000.and.n<=99999)then
        write(fname,'(i5)') n
      else
        print *, ' Error: file number not in the range [1,9999]'
        stop
      end if
      !
      filename='kslice/'//'kslice'//fname//'.h5'
      !
      call H5ReadArray(time,'time',trim(filename))
      call H5ReadArray(p,im,jm,'p',trim(filename))
      call H5ReadArray(u1,im,jm,'u1',trim(filename))
      call H5ReadArray(u2,im,jm,'u2',trim(filename))
      call H5ReadArray(u3,im,jm,'u3',trim(filename))
      call H5ReadArray(t, im,jm, 't',trim(filename))
      !
      call H5ReadArray(du1(1,:,:),im,jm,'dudx',trim(filename))
      call H5ReadArray(du1(2,:,:),im,jm,'dudy',trim(filename))
      ! call H5ReadArray(du1(3,:,:),im,jm,'du1dz',trim(filename))
      call H5ReadArray(du2(1,:,:),im,jm,'dvdx',trim(filename))
      call H5ReadArray(du2(2,:,:),im,jm,'dvdy',trim(filename))
      ! call H5ReadArray(du2(3,:,:),im,jm,'du2dz',trim(filename))
      ! call H5ReadArray(du3(1,:,:),im,jm,'du3dx',trim(filename))
      ! call H5ReadArray(du3(2,:,:),im,jm,'du3dy',trim(filename))
      ! call H5ReadArray(du3(3,:,:),im,jm,'du3dz',trim(filename))
      !
      ro=p/t*const2
      ! p=ro*t/const2
      !
      ! du1=grad_xy(u1,x2d,y2d)
      ! du2=grad_xy(u2)
      !
      dp=grad_xy(p,x2d,y2d)
      divp=sqrt(dp(1,:,:)**2+dp(2,:,:)**2)
      !
      omegaz=du2(1,:,:)-du1(2,:,:)
      !
      div=du1(1,:,:)+du2(2,:,:)
      !
      rns=rnsxy(ro)
      !
      ma=sqrt((u1**2+u2**2)/t)*mach
      !
      ! div=du1(1,:,:)+du2(2,:,:)+du3(3,:,:)
      !
      ! omegax=du3(2,:,:)-du2(3,:,:)
      ! omegay=du1(3,:,:)-du3(1,:,:)
      ! omegaz=du2(1,:,:)-du1(2,:,:)
      !
      rns=rnsxy(ro)
      !
      ix=im-50
      jy=jm-20
      !
      call writetecbin('Results/tecxy'//fname//'.plt',                 &
                              x2d(0:ix,0:jy),'x',y2d(0:ix,0:jy),'y',   &
                              ro(0:ix,0:jy),'ro',u1(0:ix,0:jy),'u',    &
                              u2(0:ix,0:jy),'v',  p(0:ix,0:jy),'p',    &
                              t(0:ix,0:jy),'t',                        &
                              div(0:ix,0:jy),'div',                    &
                              rns(0:ix,0:jy),'rns',                    &
                              omegaz(0:ix,0:jy),'omegaz',              &
                              divp(0:ix,0:jy),'divp',                  &
                              ma(0:ix,0:jy),'ma',                      &
                                                                  ix,jy)
      !
    enddo
    !
  end subroutine ksliceview
  !+-------------------------------------------------------------------+
  !| The end of the subroutine ksliceview.                             |
  !+-------------------------------------------------------------------+
  !
  subroutine isliceview(nsta,nend,islice)
    !
    use gradsolver
    use commvardefine,only: im,jm,km,const2,gridfile
    use basicfunction,only: rnsxy
    !
    ! argument
    integer,intent(in) :: nsta,nend,islice
    !
    ! local data
    real(8),allocatable :: x2d(:,:),y2d(:,:),z2d(:,:)
    real(8),allocatable,dimension(:,:) :: royz,u1yz,u2yz,u3yz,pyz,tyz, &
                                          div,rns
    real(8),allocatable :: du1(:,:,:),du2(:,:,:)
    real(8) :: time
    integer :: n,ix,jy
    character(len=5) :: fname
    character(len=255) :: filename
    !
    allocate(x2d(0:jm,0:km),y2d(0:jm,0:km),z2d(0:jm,0:km))
    call h5_read2dfrom3d(x2d,im,jm,km,'x',trim(gridfile),islice=islice)
    call h5_read2dfrom3d(y2d,im,jm,km,'y',trim(gridfile),islice=islice)
    call h5_read2dfrom3d(z2d,im,jm,km,'z',trim(gridfile),islice=islice)
    !
    allocate(royz(0:jm,0:km),u1yz(0:jm,0:km),u2yz(0:jm,0:km),          &
             u3yz(0:jm,0:km),pyz(0:jm,0:km), tyz(0:jm,0:km)            )
    !
    do n=nsta,nend
      !
      write(fname,'(i5.5)')n
      !
      filename='islice/islice'//fname//'.h5'
      print*,' >> data:',trim(filename)
      !
      call H5ReadArray(time,'time',trim(filename))
      print*,' ** time=',time
      call H5ReadArray(royz,jm,km,'ro',trim(filename))
      call H5ReadArray(u1yz,jm,km,'u1',trim(filename))
      call H5ReadArray(u2yz,jm,km,'u2',trim(filename))
      call H5ReadArray(u3yz,jm,km,'u3',trim(filename))
      call H5ReadArray(tyz, jm,km, 't',trim(filename))
      !
      pyz=royz*tyz/const2
      !
      call writetecbin('./tecyz'//fname//'.plt',                       &
                                x2d,'x',y2d,'y',z2d,'z',royz,'ro',     &
                                u1yz,'u',u2yz,'v',u3yz,'w',pyz,'p',    &
                                tyz,'t',jm,km)
      !
    enddo
    !
  end subroutine isliceview
  !+-------------------------------------------------------------------+
  !| The end of the subroutine isliceview.                             |
  !+-------------------------------------------------------------------+
  !
  subroutine jsliceview(nsta,nend,jslice)
    !
    use gradsolver
    use commvardefine,only: im,jm,km,const2
    use basicfunction,only: rnsxy
    !
    ! argument
    integer,intent(in) :: nsta,nend,jslice
    !
    ! local data
    real(8),allocatable :: x2d(:,:),y2d(:,:),z2d(:,:)
    real(8),allocatable,dimension(:,:) :: roxz,u1xz,u2xz,u3xz,pxz,txz, &
                                          us,un
    real(8),allocatable :: dudy(:,:),dwdy(:,:)
    real(8),allocatable,dimension(:) :: rom,u1m,u2m,u3m,pm,tm,alfa
    real(8) :: time,nx,ny
    integer :: n,ix,jy,i,j,k
    character(len=5) :: fname
    character(len=255) :: filename
    !
    allocate(x2d(0:im,0:km),y2d(0:im,0:km),z2d(0:im,0:km))
    call h5_read2dfrom3d(x2d,im,im,km,'x','datin/grid.h5',jslice=jslice)
    call h5_read2dfrom3d(y2d,im,im,km,'y','datin/grid.h5',jslice=jslice)
    call h5_read2dfrom3d(z2d,im,im,km,'z','datin/grid.h5',jslice=jslice)
    !
    allocate(roxz(0:im,0:km),u1xz(0:im,0:km),u2xz(0:im,0:km),          &
             u3xz(0:im,0:km),pxz(0:im,0:km), txz(0:im,0:km),           &
             dudy(0:im,0:km),dwdy(0:im,0:km)          )
    !
    ! allocate(rom(0:im),u1m(0:im),u2m(0:im),u3m(0:im),pm(0:im),tm(0:im))
    ! rom=0.d0
    ! u1m=0.d0
    ! u2m=0.d0
    ! u3m=0.d0
    !  pm=0.d0
    !  tm=0.d0
    ! do n=nsta,nend
    !   !
    !   write(fname,'(i5.5)')n
    !   !
    !   filename='wall/jslice'//fname//'.h5'
    !   print*,' >> data:',trim(filename)
    !   !
    !   call H5ReadArray(time,'time',trim(filename))
    !   print*,' ** time=',time
    !   call H5ReadArray(roxz,im,km,'ro',trim(filename))
    !   call H5ReadArray(u1xz,im,km,'u1',trim(filename))
    !   call H5ReadArray(u2xz,im,km,'u2',trim(filename))
    !   call H5ReadArray(u3xz,im,km,'u3',trim(filename))
    !   call H5ReadArray(txz, im,km, 't',trim(filename))
    !   pxz=roxz*txz/const2
    !   !
    !   do k=1,km
    !     rom(:)=rom(:)+roxz(:,k)
    !     u1m(:)=u1m(:)+roxz(:,k)*u1xz(:,k)
    !     u2m(:)=u2m(:)+roxz(:,k)*u2xz(:,k)
    !     u3m(:)=u3m(:)+roxz(:,k)*u3xz(:,k)
    !      tm(:)=tm(:) +roxz(:,k)*txz(:,k)
    !      pm(:)=pm(:) +pxz(:,k)
    !   enddo
    !   !
    ! enddo
    ! !
    ! u1m=u1m/rom
    ! u2m=u2m/rom
    ! u3m=u3m/rom
    !  tm= tm/rom
    !  pm= pm/real(nend-nsta+1,8)
    ! rom=rom/real(nend-nsta+1,8)
    ! !
    ! print*,' ** mean flow calculated.'
    !
    do n=nsta,nend
      !
      write(fname,'(i5.5)')n
      !
      filename='jslice/jslice'//fname//'.h5'
      print*,' >> data:',trim(filename)
      !
      call H5ReadArray(time,'time',trim(filename))
      print*,' ** time=',time
      call H5ReadArray(roxz,im,km,'ro',trim(filename))
      call H5ReadArray(u1xz,im,km,'u1',trim(filename))
      call H5ReadArray(u2xz,im,km,'u2',trim(filename))
      call H5ReadArray(u3xz,im,km,'u3',trim(filename))
      call H5ReadArray(pxz, im,km, 'p',trim(filename))
      call H5ReadArray(txz, im,km, 't',trim(filename))
      ! call H5ReadArray(dudy,im,km,'du1dy',trim(filename))
      ! call H5ReadArray(dwdy,im,km,'du3dy',trim(filename))
      !
      ! pxz=roxz*txz/const2
      !
      ! do k=1,km
      !   roxz(:,k)=roxz(:,k)-rom(:)
      !   u1xz(:,k)=u1xz(:,k)-u1m(:)
      !   u2xz(:,k)=u2xz(:,k)-u2m(:)
      !   u3xz(:,k)=u3xz(:,k)-u3m(:)
      !    txz(:,k)= txz(:,k) -tm(:)
      !    pxz(:,k)= pxz(:,k) -pm(:)
      ! enddo
      !

      call writetecbin('./Results/tecxz'//fname//'.plt',               &
                                x2d,'x',y2d,'y',z2d,'z',roxz,'ro"',    &
                                u1xz,'u',u2xz,'v',u3xz,'w',txz,'T',    &
                                pxz,'p',im,km)

      ! allocate(us(0:im,0:km),un(0:im,0:km),alfa(0:im) )
      ! !
      ! do i=0,im
      !   alfa(i)=asin(y2d(i,0)/sqrt((x2d(i,0)-1.d0)**2+y2d(i,0)**2))
      !   nx=-cos(alfa(i))
      !   ny= sin(alfa(i))
      !   do k=0,km
      !     us(i,k)= ny*u1xz(i,k)-nx*u2xz(i,k)
      !     un(i,k)= nx*u1xz(i,k)+ny*u2xz(i,k)
      !   enddo
      ! enddo
      ! !
      ! call writetecbin('./Results/tecxz'//fname//'.plt',               &
      !                           x2d,'x',y2d,'y',z2d,'z',roxz,'ro"',    &
      !                           us,'us"',un,'un"',u3xz,'w"',pxz,"p'",  &
      !                           txz,'t"',im,km)
      !
    enddo
    !
  end subroutine jsliceview
  !+-------------------------------------------------------------------+
  !| The end of the subroutine jsliceview.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write a xdmf head file for             |
  !| visulisation of flow field.                                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 01-07-2022  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine xdmfwriter(dir,filename,gridfile,var1,var1name,var2,var2name, &
                                              var3,var3name,var4,var4name, &
                                              var5,var5name,var6,var6name, im,jm,km)
    !
    ! arguments
    character(len=*),intent(in) :: dir,filename,gridfile
    real(8),intent(in),dimension(:,:,:) :: var1
    character(len=*),intent(in) :: var1name
    integer :: im,jm,km
    real(8),intent(in),dimension(:,:,:),optional :: var2,var3,var4,var5,var6
    character(len=*),intent(in),optional :: var2name,var3name,var4name,var5name,var6name
    !
    real(4),allocatable :: bufr4(:,:,:)
    !
    ! local data
    integer :: fh,i
    !
    allocate(bufr4(0:im,0:jm,0:km))
    !
    fh=18
    !
    ! write the head of xdmf and grid
    open(fh,file=dir//filename//'.xdmf',form='formatted')
    write(fh,'(A)')'<?xml version="1.0" ?>'
    write(fh,'(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(fh,'(A)')'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(fh,'(A)')'  <Domain>'
    !
    write(fh,'(A,3(1X,I0),A)')'    <Topology name="topo" TopologyType="3DSMESH" Dimensions="',  &
                              km+1,jm+1,im+1,'"> </Topology>'
    write(fh,'(A)')'    <Geometry name="geo" Type="X_Y_Z">'
    write(fh,'(A)')'      <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,3(1X,I0),3(A))')'                Dimensions="',km+1,jm+1,im+1,'"> ',gridfile,':x </DataItem>'
    write(fh,'(A)')'      <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,3(1X,I0),3(A))')'                Dimensions="',km+1,jm+1,im+1,'"> ',gridfile,':y </DataItem>'
    write(fh,'(A)')'      <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,3(1X,I0),3(A))')'                Dimensions="',km+1,jm+1,im+1,'"> ',gridfile,':z </DataItem>'
    write(fh,'(A)')'    </Geometry>'
    !
    write(fh,'(A)')'    <Grid Name="001" GridType="Uniform">'
    write(fh,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    write(fh,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    !
    write(fh,'(3(A))')'      <Attribute Name="',var1name,'" Center="Node">'
    write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,3(1X,I0),5(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'.h5',':',var1name,'</DataItem>'
    write(fh,'(A)')'      </Attribute>'
    bufr4=var1
    call H5WriteArray(bufr4,im,jm,km,var1name,dir//filename//'.h5')

    if(present(var2)) then
      write(fh,'(3(A))')'      <Attribute Name="',var2name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),5(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'.h5',':',var2name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var2
      call H5WriteArray(bufr4,im,jm,km,var2name,dir//filename//'.h5')
    endif
    if(present(var3)) then
      write(fh,'(3(A))')'      <Attribute Name="',var3name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),5(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'.h5',':',var3name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var3
      call H5WriteArray(bufr4,im,jm,km,var3name,dir//filename//'.h5')
    endif
    if(present(var4)) then
      write(fh,'(3(A))')'      <Attribute Name="',var4name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),5(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'.h5',':',var4name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var4
      call H5WriteArray(bufr4,im,jm,km,var4name,dir//filename//'.h5')
    endif
    if(present(var5)) then
      write(fh,'(3(A))')'      <Attribute Name="',var5name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),5(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'.h5',':',var5name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var5
      call H5WriteArray(bufr4,im,jm,km,var5name,dir//filename//'.h5')
    endif
    if(present(var6)) then
      write(fh,'(3(A))')'      <Attribute Name="',var6name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),5(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'.h5',':',var6name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var6
      call H5WriteArray(bufr4,im,jm,km,var6name,dir//filename//'.h5')
    endif

    write(fh,'(A)')'     </Grid>'
    !
    write(fh,'(A)')'  </Domain>'
    write(fh,'(A)')'</Xdmf>'
    !
    close(fh)
    !
    !
  end subroutine xdmfwriter
  !+-------------------------------------------------------------------+
  !| The end of the subroutine xdmfwriter.                             |
  !+-------------------------------------------------------------------+
  !
end module fieldview
!+---------------------------------------------------------------------+
!| The end of the module fieldview.                                    |
!+---------------------------------------------------------------------+
