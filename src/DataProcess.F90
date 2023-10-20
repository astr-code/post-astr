!+---------------------------------------------------------------------+
!| This module contains subroutines of data processing.                |
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2020-04-16                                     |
!+---------------------------------------------------------------------+
module dataprocess
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to convered slice* to hdf5.               |
  !+-------------------------------------------------------------------+
  subroutine kslice2h5(nsta,nend)
    !
    use commvardefine, only: ink,jnk,knk,imp,jmp,kmp,iop,jop,kop,      &
                             im,jm,km,rankmax
    use h5readwrite
    !
    integer,intent(in) :: nsta,nend
    !
    integer :: nrank,ni,nj,nk,ima,jma,kma,krks
    integer :: i,j,k,n,irk,jrk,krk
    real(8) :: time
    logical :: lfilalive
    character(len=5) :: fname
    character(len=8) :: mpirankname
    character(len=255) :: filename
    real(8) :: ro(0:im,0:jm),u1(0:im,0:jm),u2(0:im,0:jm),u3(0:im,0:jm),&
               t(0:im,0:jm),p(0:im,0:jm),du1(3,0:im,0:jm),             &
               du2(3,0:im,0:jm),du3(3,0:im,0:jm)
    real(8),allocatable :: var(:,:,:),dvar(:,:,:,:)
    integer :: error(15)
    !
    !$ save var,dvar
    !
    !$OMP threadprivate(var,dvar)
    !
    print*,' ** convert kslice data from ',nsta,' to ',nend
    !determin krks
    n=nsta
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
    do nrank=0,rankmax
      !
        if(nrank>=0 .and. nrank<=9)then
        write(mpirankname,'(7h.000000,i1)') nrank
      elseif(nrank>=10 .and. nrank<=99)then
        write(mpirankname,'(6h.00000,i2)') nrank
      elseif(nrank>=100 .and. nrank<=999)then
        write(mpirankname,'(5h.0000,i3)') nrank
      elseif(nrank>=1000 .and. nrank<=9999)then
        write(mpirankname,'(4h.000,i4)') nrank
      elseif(nrank>=10000 .and. nrank<=99999)then
        write(mpirankname,'(3h.00,i5)') nrank
      elseif(nrank>=100000 .and. nrank<=999999)then
        write(mpirankname,'(2h.0,i6)') nrank
      elseif(nrank>=1000000 .and. nrank<=9999999)then
        write(mpirankname,'(1h.,i7)') nrank
      else
        print *, ' !! Error: rank number not in the range [0,9999999]'
        stop
      end if
      !
      filename='Outdat/kslice'//fname//mpirankname
      !
      irk=ink(nrank)
      jrk=jnk(nrank)
      krk=knk(nrank)
      !
      inquire(file=trim(filename),exist=lfilalive)
      !
      ! print*,' ** checking file:',trim(filename),' -- ',lfilalive
      !
      if(lfilalive) then
        !
        krks=krk
        !
        exit
        !
      endif
      !
    end do
    !
    if(.not. lfilalive) then
      print*,' !! cant find anyfile for Outdat/kslice',fname
      stop
    else
      print*,' ** slices were extracted at the k_rank=',krks
    endif
    !
    call sleep(3)
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
      !$OMP parallel default(shared) private(nrank,i,j,k,irk,jrk,krk,  &
      !$OMP                ni,nj,nk,ima,jma,kma,mpirankname,filename)
      !$OMP do
      do nrank=0,rankmax
        !
        if(nrank>=0 .and. nrank<=9)then
          write(mpirankname,'(7h.000000,i1)') nrank
        elseif(nrank>=10 .and. nrank<=99)then
          write(mpirankname,'(6h.00000,i2)') nrank
        elseif(nrank>=100 .and. nrank<=999)then
          write(mpirankname,'(5h.0000,i3)') nrank
        elseif(nrank>=1000 .and. nrank<=9999)then
          write(mpirankname,'(4h.000,i4)') nrank
        elseif(nrank>=10000 .and. nrank<=99999)then
          write(mpirankname,'(3h.00,i5)') nrank
        elseif(nrank>=100000 .and. nrank<=999999)then
          write(mpirankname,'(2h.0,i6)') nrank
        elseif(nrank>=1000000 .and. nrank<=9999999)then
          write(mpirankname,'(1h.,i7)') nrank
        else
          print *, ' !! Error: rank number not in the range [0,9999999]'
          stop
        end if
        !
        irk=ink(nrank)
        jrk=jnk(nrank)
        krk=knk(nrank)
        !
        ni=imp(nrank)
        nj=jmp(nrank)
        nk=kmp(nrank)
        !
        ima=iop(nrank)
        jma=jop(nrank)
        kma=kop(nrank)
        !
        filename='Outdat/kslice'//fname//mpirankname
        !
        if(krk==krks) then
          !
          allocate( var(0:ni,0:nj,5),dvar(3,0:ni,0:nj,3) )
          !
          open(12+nrank,file=trim(filename),form='unformatted',        &
                                            status='old',action ='read')
          read(12+nrank)time
          read(12+nrank)var(:,:,1),var(:,:,2),var(:,:,3),var(:,:,4),   &
                        var(:,:,5)
          ! read(12+nrank)dvar(:,:,:,1),dvar(:,:,:,2),dvar(:,:,:,3)
          close(12+nrank)
          !
          write(*,'(1A1,A,A,$)')char(13),'  >> ',trim(filename)
          !
          do j=0,nj
          do i=0,ni
            ro(i+ima,j+jma)=var(i,j,1)
            u1(i+ima,j+jma)=var(i,j,2)
            u2(i+ima,j+jma)=var(i,j,3)
            u3(i+ima,j+jma)=var(i,j,4)
             t(i+ima,j+jma)=var(i,j,5)
            !
            ! du1(1,i+ima,j+jma)=dvar(1,i,j,1)
            ! du1(2,i+ima,j+jma)=dvar(2,i,j,1)
            ! du1(3,i+ima,j+jma)=dvar(3,i,j,1)
            ! du2(1,i+ima,j+jma)=dvar(1,i,j,2)
            ! du2(2,i+ima,j+jma)=dvar(2,i,j,2)
            ! du2(3,i+ima,j+jma)=dvar(3,i,j,2)
            ! du3(1,i+ima,j+jma)=dvar(1,i,j,3)
            ! du3(2,i+ima,j+jma)=dvar(2,i,j,3)
            ! du3(3,i+ima,j+jma)=dvar(3,i,j,3)
            !
          end do
          end do
          !
          deallocate(var,dvar)
          !
        end if
        !
      end do
      !$OMP end do
      !$OMP end parallel
      !
      call H5WriteArray(time,'time','Outdat/kslice'//fname//'.h5'  ,ierr=error(1))
      call H5WriteArray(ro,im,jm,'ro','Outdat/kslice'//fname//'.h5',ierr=error(2))
      call H5WriteArray(u1,im,jm,'u1','Outdat/kslice'//fname//'.h5',ierr=error(3))
      call H5WriteArray(u2,im,jm,'u2','Outdat/kslice'//fname//'.h5',ierr=error(4))
      call H5WriteArray(u3,im,jm,'u3','Outdat/kslice'//fname//'.h5',ierr=error(5))
      call H5WriteArray(t,im,jm,  't','Outdat/kslice'//fname//'.h5',ierr=error(6))
      !
      ! call H5WriteArray(du1(1,:,:),im,jm,'du1dx','Outdat/kslice'//fname//'.h5',ierr=error(7))
      ! call H5WriteArray(du1(2,:,:),im,jm,'du1dy','Outdat/kslice'//fname//'.h5',ierr=error(8))
      ! call H5WriteArray(du1(3,:,:),im,jm,'du1dz','Outdat/kslice'//fname//'.h5',ierr=error(9))
      ! call H5WriteArray(du2(1,:,:),im,jm,'du2dx','Outdat/kslice'//fname//'.h5',ierr=error(10))
      ! call H5WriteArray(du2(2,:,:),im,jm,'du2dy','Outdat/kslice'//fname//'.h5',ierr=error(11))
      ! call H5WriteArray(du2(3,:,:),im,jm,'du2dz','Outdat/kslice'//fname//'.h5',ierr=error(12))
      ! call H5WriteArray(du3(1,:,:),im,jm,'du3dx','Outdat/kslice'//fname//'.h5',ierr=error(13))
      ! call H5WriteArray(du3(2,:,:),im,jm,'du3dy','Outdat/kslice'//fname//'.h5',ierr=error(14))
      ! call H5WriteArray(du3(3,:,:),im,jm,'du3dz','Outdat/kslice'//fname//'.h5',ierr=error(15))
      !
      if(any(error .ne. 0)) then
        print*,error
        print*,'Outdat/kslice',fname,'.h5 not generated correctly'
        stop ' !! kslice2h5 !!'
      else
        !
        filename='Outdat/kslice'//fname//'.0*'
        call system('rm '//trim(filename))
        print*,' ** ',trim(filename),'deleted'
         !
      endif
      !
    end do
    !
  end subroutine kslice2h5
  !
  subroutine inflowpro(nsta,nend)
    !
    use commvardefine, only: im,jm,km,x,y,z
    use interpolation
    use h5readwrite
    use writetec
    !
    integer,intent(in) :: nsta,nend
    !
    ! local data
    real(8),allocatable,dimension(:,:) :: ro,u1,u2,u3,p,t
    real(8),allocatable :: x2dt(:,:),z2dt(:,:),y1dt(:)
    real(8),allocatable,dimension(:,:) :: ros,u1s,u2s,u3s,ps,ts
    real(8),allocatable,dimension(:,:,:) :: xs,ys,zs
    real(8) :: time
    integer :: i,j,k,jm2,km2,n
    character(len=4) :: fname1
    character(len=5) :: fname2
    logical :: lfilalive
    !
    allocate( x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km),   &
              ro(0:jm,0:km),u1(0:jm,0:km),u2(0:jm,0:km),u3(0:jm,0:km), &
              p(0:jm,0:km),t(0:jm,0:km) )
    !
    allocate( x2dt(0:im,0:km),z2dt(0:im,0:km),y1dt(0:jm) )
    open(12,file='datin/grid.dat',form='unformatted')
    read(12)x2dt
    read(12)y1dt
    read(12)z2dt
    close(12)
    print*,' >> grid.dat'
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      x(i,j,k)=x2dt(i,k)
      y(i,j,k)=y1dt(j)
      z(i,j,k)=z2dt(i,k)
    end do
    end do
    end do
    !
    jm2=jm
    km2=128
    allocate( xs(0:im,0:jm2,0:km2),ys(0:im,0:jm2,0:km2),zs(0:im,0:jm2,0:km2), &
              ros(0:jm2,0:km2),u1s(0:jm2,0:km2),u2s(0:jm2,0:km2),      &
              u3s(0:jm2,0:km2),ps(0:jm2,0:km2),ts(0:jm2,0:km2) )
    do j=0,jm2
      ys(0,j,:)=y1dt(j)
    enddo
    !
    do k=0,km2
      zs(0,:,k)=10.d0/dble(km2)*dble(k)
    enddo
    !
    do i=0,im
      xs(i,:,:)=10.d0/dble(im)*dble(i)
      ys(i,:,:)=ys(0,:,:)
      zs(i,:,:)=zs(0,:,:)
    enddo
    deallocate(x2dt,z2dt,y1dt)
    !
    ! call H5WriteArray(xs,im,jm2,km2,'x','datin/grid.h5')
    ! call H5WriteArray(ys,im,jm2,km2,'y','datin/grid.h5')
    ! call H5WriteArray(zs,im,jm2,km2,'z','datin/grid.h5')
    ! !
    do n=nsta,nend
      !
      write(fname1,'(I4.4)')n
      write(fname2,'(I5.5)')n
      !
      inquire(file='inflow/indata'//fname1,exist=lfilalive)
      !
      if(lfilalive) then
        open(12,file='inflow/indata'//fname1,form='unformatted')
        read(12)time
        read(12)ro,u1,u2,u3,t
        close(12)
        print*,' >> inflow/indata',fname1
      else
        print*,' ** inflow/indata',fname1,' is missing.'
        cycle
      endif
      !
      call regularlinearinterp(y(0,:,:),z(0,:,:),ro,u1,u2,u3,t,jm,km,  &
                             ys(0,:,:),zs(0,:,:),ros,u1s,u2s,u3s,ts,jm2,km2 )
    
      call H5WriteArray(time,'time','inflow2/islice'//fname2//'.h5'    )
      call H5WriteArray(ros,jm2,km2,'ro','inflow2/islice'//fname2//'.h5')
      call H5WriteArray(u1s,jm2,km2,'u1','inflow2/islice'//fname2//'.h5')
      call H5WriteArray(u2s,jm2,km2,'u2','inflow2/islice'//fname2//'.h5')
      call H5WriteArray(u3s,jm2,km2,'u3','inflow2/islice'//fname2//'.h5')
      call H5WriteArray(ts, jm2,km2, 't','inflow2/islice'//fname2//'.h5')
      !
    enddo

    call writetecbin( 'tecinflow.plt',ys(0,:,:),'y',zs(0,:,:),'z', &
                      ros,'ro',u1s,'u',u2s,'v',u3s,'w',ts,'t',jm2,km2 )
    stop
    !
  end subroutine inflowpro
  !
  subroutine isdatabase(nsta,nend,islice)
    !
    use commvardefine,only: im,jm,km,Reynolds,Mach,pinf,pi,const2,gridfile
    use h5readwrite
    use writetec
    use interpolation
    use basicfunction
    !
    integer,intent(in) :: nsta,nend,islice
    !
    real(8),allocatable,dimension(:) :: rom,u1m,u2m,u3m,pm,tm,yplus,uplus
    real(8),allocatable,dimension(:,:) :: ro,u1,u2,u3,p,t,x,y,z
    !
    character(len=5) :: fname
    character(len=3) :: jname
    character(len=255) :: filename
    !
    integer :: n,i,j,k,k1,k2,js(3),ju
    real(8) :: time,var1,var2,bl_delta99,bl_dstar,bl_theta,utaw,  &
               lvis,nsamples,shapefac,scale
    logical :: lfilalive
    !
    allocate(x(0:jm,0:km),y(0:jm,0:km),z(0:jm,0:km))
    call h5_read2dfrom3d(x,im,jm,km,'x',trim(gridfile),islice=islice)
    call h5_read2dfrom3d(y,im,jm,km,'y',trim(gridfile),islice=islice)
    call h5_read2dfrom3d(z,im,jm,km,'z',trim(gridfile),islice=islice)
    !
    ju=jm ! 2.5delta
    !
    scale=1.d0
    ! scale=1.d0/1.1694605d0
    ! scale=1.d0/7.11d0
    ! scale=1.d0/7.33863d0
    print*,' ** scale:',scale
    !
    x=x*scale
    y=y*scale
    z=z*scale
    !
    inquire(file='database/grid.h5',exist=lfilalive)
    if(.not. lfilalive) then
      call H5WriteArray(x(0:ju,:),ju,km,'x','database/grid.h5')
      call H5WriteArray(y(0:ju,:),ju,km,'y','database/grid.h5')
      call H5WriteArray(z(0:ju,:),ju,km,'z','database/grid.h5')
    endif
    !
    inquire(file='database/meanprofile.h5',exist=lfilalive)
    if(.not. lfilalive) then
      !
      allocate(rom(0:jm),u1m(0:jm),u2m(0:jm),u3m(0:jm),pm(0:jm),tm(0:jm))
      !
      call H5ReadArray(rom,jm,'ro','Results/meanprofile.h5')
      call H5ReadArray(u1m,jm,'u1','Results/meanprofile.h5')
      call H5ReadArray(u2m,jm,'u2','Results/meanprofile.h5')
      call H5ReadArray(u3m,jm,'u3','Results/meanprofile.h5')
      call H5ReadArray( tm,jm, 't','Results/meanprofile.h5')
      call H5ReadArray( pm,jm, 'p','Results/meanprofile.h5')
      !
      call H5WriteArray(rom(0:ju),ju,'ro','database/meanprofile.h5')
      call H5WriteArray(u1m(0:ju),ju,'u1','database/meanprofile.h5')
      call H5WriteArray(u2m(0:ju),ju,'u2','database/meanprofile.h5')
      call H5WriteArray(u3m(0:ju),ju,'u3','database/meanprofile.h5')
      call H5WriteArray( tm(0:ju),ju, 't','database/meanprofile.h5')
      call H5WriteArray( pm(0:ju),ju, 'p','database/meanprofile.h5')
      !
    endif
    !
    allocate(ro(0:jm,0:km),u1(0:jm,0:km),u2(0:jm,0:km),u3(0:jm,0:km),  &
              p(0:jm,0:km), t(0:jm,0:km))
    !
    do n=nsta,nend
      !
      write(fname,'(I5.5)')n
      !
      filename='islice/islice'//fname//'.h5'
      !
      inquire(file=trim(filename),exist=lfilalive)
      !
      if(lfilalive) then
        !
        call H5ReadArray(time,'time',trim(filename))
        ! call H5ReadArray(ro,jm,km,'ro',trim(filename))
        call H5ReadArray( p,jm,km, 'p',trim(filename))
        call H5ReadArray(u1,jm,km,'u1',trim(filename))
        call H5ReadArray(u2,jm,km,'u2',trim(filename))
        call H5ReadArray(u3,jm,km,'u3',trim(filename))
        call H5ReadArray( t,jm,km, 't',trim(filename))
        !
        ! p=ro*t/const2
        ro=p/t*const2
        !
        nsamples=nsamples+1.d0
      else
        print*,' ** ',trim(filename),' is missing.'
        cycle
      endif
      !
      time=time*scale
      !
      filename='database/islice/islice'//fname//'.h5'
      call H5WriteArray(time,'time',trim(filename))
      call H5WriteArray(ro(0:ju,:),ju,km,'ro',trim(filename))
      call H5WriteArray(u1(0:ju,:),ju,km,'u1',trim(filename))
      call H5WriteArray(u2(0:ju,:),ju,km,'u2',trim(filename))
      call H5WriteArray(u3(0:ju,:),ju,km,'u3',trim(filename))
      call H5WriteArray( t(0:ju,:),ju,km, 't',trim(filename))
      !
    enddo
    !
  end subroutine isdatabase
  !
  subroutine islice2h5(nsta,nend)
    !
    use commvardefine, only: ink,jnk,knk,imp,jmp,kmp,iop,jop,kop,      &
                             im,jm,km,rankmax
    use h5readwrite
    !
    integer,intent(in) :: nsta,nend
    !
    integer :: nrank,ni,nj,nk,ima,jma,kma,irks
    integer :: i,j,k,n,irk,jrk,krk
    real(8) :: time
    character(len=5) :: fname
    character(len=8) :: mpirankname
    character(len=255) :: filename,cmd
    real(8) :: royz(0:jm,0:km),u1yz(0:jm,0:km),u2yz(0:jm,0:km),        &
               u3yz(0:jm,0:km),tyz(0:jm,0:km),pyz(0:jm,0:km)
    real(8) :: du1yz(3,0:jm,0:km),du2yz(3,0:jm,0:km),du3yz(3,0:jm,0:km)
    real(8),allocatable,dimension(:,:,:) :: var
    real(8),allocatable,dimension(:,:,:,:) :: dvar
    integer :: error(15)
    real(8) :: time0
    logical :: lfilalive
    !
    time0=0.d0
    !
    n=nsta
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
    do nrank=0,rankmax
      !
      if(nrank>=0 .and. nrank<=9)then
        write(mpirankname,'(7h.000000,i1)') nrank
      elseif(nrank>=10 .and. nrank<=99)then
        write(mpirankname,'(6h.00000,i2)') nrank
      elseif(nrank>=100 .and. nrank<=999)then
        write(mpirankname,'(5h.0000,i3)') nrank
      elseif(nrank>=1000 .and. nrank<=9999)then
        write(mpirankname,'(4h.000,i4)') nrank
      elseif(nrank>=10000 .and. nrank<=99999)then
        write(mpirankname,'(3h.00,i5)') nrank
      elseif(nrank>=100000 .and. nrank<=999999)then
        write(mpirankname,'(2h.0,i6)') nrank
      elseif(nrank>=1000000 .and. nrank<=9999999)then
        write(mpirankname,'(1h.,i7)') nrank
      else
        print *, ' !! Error: rank number not in the range [0,9999999]'
        stop
      end if
      !
      irk=ink(nrank)
      jrk=jnk(nrank)
      krk=knk(nrank)
      !
      filename='Outdat/islice'//fname//mpirankname
      !
      inquire(file=trim(filename),exist=lfilalive)
      !
      if(lfilalive) then
        !
        irks=irk
        ! 
        print*,' ** file alive for',trim(filename)
        !
        exit
        !
      endif
      !
    end do
    !
    if(.not. lfilalive) then
      print*,' No file can be found for islice',fname,'*******'
      stop
    end if
    !
    do n=nsta,nend
      !
      write(fname,'(i5.5)')n
      !
      do nrank=0,rankmax
        !
        if(nrank>=0 .and. nrank<=9)then
          write(mpirankname,'(7h.000000,i1)') nrank
        elseif(nrank>=10 .and. nrank<=99)then
          write(mpirankname,'(6h.00000,i2)') nrank
        elseif(nrank>=100 .and. nrank<=999)then
          write(mpirankname,'(5h.0000,i3)') nrank
        elseif(nrank>=1000 .and. nrank<=9999)then
          write(mpirankname,'(4h.000,i4)') nrank
        elseif(nrank>=10000 .and. nrank<=99999)then
          write(mpirankname,'(3h.00,i5)') nrank
        elseif(nrank>=100000 .and. nrank<=999999)then
          write(mpirankname,'(2h.0,i6)') nrank
        elseif(nrank>=1000000 .and. nrank<=9999999)then
          write(mpirankname,'(1h.,i7)') nrank
        else
          print *, ' !! Error: rank number not in the range [0,9999999]'
          stop
        end if
        !
        irk=ink(nrank)
        jrk=jnk(nrank)
        krk=knk(nrank)
        !
        ni=imp(nrank)
        nj=jmp(nrank)
        nk=kmp(nrank)
        !
        ima=iop(nrank)
        jma=jop(nrank)
        kma=kop(nrank)
        !
        filename='Outdat/islice'//fname//mpirankname
        !
        if(irk==irks) then
            !
            allocate( var(5,0:nj,0:nk),dvar(3,3,0:nj,0:nk) )
            !
            open(12,file=trim(filename),form='unformatted')
            read(12)time
            read(12)var(1,:,:),var(2,:,:),var(3,:,:),var(4,:,:),var(5,:,:)
            read(12)dvar(1,:,:,:),dvar(2,:,:,:),dvar(3,:,:,:)
            close(12)
            do k=0,nk
            do j=0,nj
              royz(j+jma,k+kma)=var(1,j,k)
              u1yz(j+jma,k+kma)=var(2,j,k)
              u2yz(j+jma,k+kma)=var(3,j,k)
              u3yz(j+jma,k+kma)=var(4,j,k)
               tyz(j+jma,k+kma)=var(5,j,k)
              !
              du1yz(1,j+jma,k+kma)=dvar(1,1,j,k)
              du1yz(2,j+jma,k+kma)=dvar(1,2,j,k)
              du1yz(3,j+jma,k+kma)=dvar(1,3,j,k)
              !
              du2yz(1,j+jma,k+kma)=dvar(2,1,j,k)
              du2yz(2,j+jma,k+kma)=dvar(2,2,j,k)
              du2yz(3,j+jma,k+kma)=dvar(2,3,j,k)
              !
              du3yz(1,j+jma,k+kma)=dvar(3,1,j,k)
              du3yz(2,j+jma,k+kma)=dvar(3,2,j,k)
              du3yz(3,j+jma,k+kma)=dvar(3,3,j,k)
              !
            end do
            end do
            !
            deallocate(var,dvar)
            print*,' >> ',trim(filename),'...done.'
            !
            !if(n==nsta) time0=time(1)
            !
         end if
         !
       end do
       !
       write(fname,'(i5.5)')n+1
       call H5WriteArray(time,'time','inflow/islice'//fname//'.h5'    ,ierr=error(1))
       call H5WriteArray(royz,jm,km,'ro','inflow/islice'//fname//'.h5',ierr=error(2))
       call H5WriteArray(u1yz,jm,km,'u1','inflow/islice'//fname//'.h5',ierr=error(3))
       call H5WriteArray(u2yz,jm,km,'u2','inflow/islice'//fname//'.h5',ierr=error(4))
       call H5WriteArray(u3yz,jm,km,'u3','inflow/islice'//fname//'.h5',ierr=error(5))
       call H5WriteArray(tyz,jm,km,  't','inflow/islice'//fname//'.h5',ierr=error(6))
       !
       call H5WriteArray(du1yz(1,:,:),jm,km,'du1dx','inflow/islice'//fname//'.h5',ierr=error(7))
       call H5WriteArray(du1yz(2,:,:),jm,km,'du1dy','inflow/islice'//fname//'.h5',ierr=error(8))
       call H5WriteArray(du1yz(3,:,:),jm,km,'du1dz','inflow/islice'//fname//'.h5',ierr=error(9))
       call H5WriteArray(du2yz(1,:,:),jm,km,'du2dx','inflow/islice'//fname//'.h5',ierr=error(10))
       call H5WriteArray(du2yz(2,:,:),jm,km,'du2dy','inflow/islice'//fname//'.h5',ierr=error(11))
       call H5WriteArray(du2yz(3,:,:),jm,km,'du2dz','inflow/islice'//fname//'.h5',ierr=error(12))
       call H5WriteArray(du3yz(1,:,:),jm,km,'du3dx','inflow/islice'//fname//'.h5',ierr=error(13))
       call H5WriteArray(du3yz(2,:,:),jm,km,'du3dy','inflow/islice'//fname//'.h5',ierr=error(14))
       call H5WriteArray(du3yz(3,:,:),jm,km,'du3dz','inflow/islice'//fname//'.h5',ierr=error(15))
       !
       if(any(error .ne. 0)) then
         print*,error
         print*,'inflow/islice',fname,'.h5 not generated correctly'
         stop ' !! islice2h5 !!'
       else
         write(fname,'(i5.5)')n
         filename='Outdat/islice'//fname//'.0*'
         call system('rm -v '//trim(filename))
       endif
       !
       ! cmd='rm -v Outdat/islice'//fname//'.*'
       ! call system(trim(cmd))
       !
    end do
    !
  end subroutine islice2h5
  !
  subroutine jslice2h5(nsta,nend)
    !
    use commvardefine, only: ink,jnk,knk,imp,jmp,kmp,iop,jop,kop,      &
                             im,jm,km,rankmax
    use h5readwrite
    !
    integer,intent(in) :: nsta,nend
    !
    integer :: nrank,ni,nj,nk,ima,jma,kma,jrks
    integer :: i,j,k,n,irk,jrk,krk
    real(8) :: time
    character(len=5) :: fname
    character(len=8) :: mpirankname
    character(len=255) :: filename,cmd
    real(8) :: ro(0:im,0:km),u1(0:im,0:km),u2(0:im,0:km),              &
               u3(0:im,0:km),t(0:im,0:km),p(0:im,0:km),                &
               du1(3,0:im,0:km),du2(3,0:im,0:km),du3(3,0:im,0:km)
    real(8),allocatable :: var(:,:,:),dvar(:,:,:,:)
    integer :: error(15)
    real(8) :: time0
    logical :: lfilalive
    !
    !$ save var,dvar
    !
    !$OMP threadprivate(var,dvar)
    !
    time0=0.d0
    !
    n=nsta
    !
    write(fname,'(i5.5)')n
    !
    do nrank=0,rankmax
      !
      write(mpirankname,'(1h.,i7.7)')nrank
      !
      irk=ink(nrank)
      jrk=jnk(nrank)
      krk=knk(nrank)
      !
      filename='Outdat/jslice'//fname//mpirankname
      !
      inquire(file=trim(filename),exist=lfilalive)
      !
      if(lfilalive) then
        !
        jrks=jrk
        ! 
        print*,' ** file alive for',trim(filename)
        !
        exit
        !
      endif
      !
    end do
    !
    if(.not. lfilalive) then
      print*,' No file can be found for jslice',fname,'*******'
      stop
    end if
    !
    do n=nsta,nend
      !
      write(fname,'(i5.5)')n
      !
      !$OMP parallel default(shared) private(nrank,i,j,k,irk,jrk,krk,  &
      !$OMP                ni,nj,nk,ima,jma,kma,mpirankname,filename)
      !$OMP do
      do nrank=0,rankmax
        !
        write(mpirankname,'(1h.,i7.7)')nrank
        !
        irk=ink(nrank)
        jrk=jnk(nrank)
        krk=knk(nrank)
        !
        ni=imp(nrank)
        nj=jmp(nrank)
        nk=kmp(nrank)
        !
        ima=iop(nrank)
        jma=jop(nrank)
        kma=kop(nrank)
        !
        filename='Outdat/jslice'//fname//mpirankname
        !
        if(jrk==jrks) then
           !
           allocate( var(0:ni,0:nk,5),dvar(3,0:ni,0:nk,3) )
           !
           open(12+nrank,file=trim(filename),form='unformatted')
           read(12+nrank)time
           read(12+nrank)var(:,:,1),var(:,:,2),var(:,:,3),var(:,:,4),  &
                         var(:,:,5)
           read(12+nrank)dvar(:,:,:,1),dvar(:,:,:,2),dvar(:,:,:,3)
           close(12+nrank)
           !
           write(*,'(1A1,A,A,$)')char(13),'  >> ',trim(filename)

           !
           do k=0,nk
           do i=0,ni
             ro(i+ima,k+kma)=var(i,k,1)
             u1(i+ima,k+kma)=var(i,k,2)
             u2(i+ima,k+kma)=var(i,k,3)
             u3(i+ima,k+kma)=var(i,k,4)
              t(i+ima,k+kma)=var(i,k,5)
             !
             du1(1,i+ima,k+kma)=dvar(1,i,k,1)
             du1(2,i+ima,k+kma)=dvar(2,i,k,1)
             du1(3,i+ima,k+kma)=dvar(3,i,k,1)
             du2(1,i+ima,k+kma)=dvar(1,i,k,2)
             du2(2,i+ima,k+kma)=dvar(2,i,k,2)
             du2(3,i+ima,k+kma)=dvar(3,i,k,2)
             du3(1,i+ima,k+kma)=dvar(1,i,k,3)
             du3(2,i+ima,k+kma)=dvar(2,i,k,3)
             du3(3,i+ima,k+kma)=dvar(3,i,k,3)
             !
           end do
           end do
           !
           deallocate(var,dvar)
           !
           !if(n==nsta) time0=time(1)
           !
        end if
        !
      end do
      !$OMP end do
      !$OMP end parallel
      !
      call H5WriteArray(time,'time','Outdat/jslice'//fname//'.h5',ierr=error(1))
      !
      call H5WriteArray(ro,im,km,'ro','Outdat/jslice'//fname//'.h5',ierr=error(2))
      call H5WriteArray(u1,im,km,'u1','Outdat/jslice'//fname//'.h5',ierr=error(3))
      call H5WriteArray(u2,im,km,'u2','Outdat/jslice'//fname//'.h5',ierr=error(4))
      call H5WriteArray(u3,im,km,'u3','Outdat/jslice'//fname//'.h5',ierr=error(5))
      call H5WriteArray(t, im,km,'t','Outdat/jslice'//fname//'.h5',ierr=error(6))
      !
      call H5WriteArray(du1(1,:,:),im,km,'du1dx','Outdat/jslice'//fname//'.h5',ierr=error(7))
      call H5WriteArray(du1(2,:,:),im,km,'du1dy','Outdat/jslice'//fname//'.h5',ierr=error(8))
      call H5WriteArray(du1(3,:,:),im,km,'du1dz','Outdat/jslice'//fname//'.h5',ierr=error(9))
      call H5WriteArray(du2(1,:,:),im,km,'du2dx','Outdat/jslice'//fname//'.h5',ierr=error(10))
      call H5WriteArray(du2(2,:,:),im,km,'du2dy','Outdat/jslice'//fname//'.h5',ierr=error(11))
      call H5WriteArray(du2(3,:,:),im,km,'du2dz','Outdat/jslice'//fname//'.h5',ierr=error(12))
      call H5WriteArray(du3(1,:,:),im,km,'du3dx','Outdat/jslice'//fname//'.h5',ierr=error(13))
      call H5WriteArray(du3(2,:,:),im,km,'du3dy','Outdat/jslice'//fname//'.h5',ierr=error(14))
      call H5WriteArray(du3(3,:,:),im,km,'du3dz','Outdat/jslice'//fname//'.h5',ierr=error(15))
      !
      if(any(error .ne. 0)) then
        print*,error
        print*,'wall/jslice',fname,'.h5 not generated correctly'
        stop ' !! jslice2h5 !!'
      else
        write(fname,'(i5.5)')n
        filename='Outdat/jslice'//fname//'.0*'
        call system('rm '//trim(filename))
        print*,' ** ',trim(filename),'deleted'
      endif
      !
    end do
    !
  end subroutine jslice2h5
  !+-------------------------------------------------------------------+
  !| The end of the subroutine slice2h5.                               |
  !+-------------------------------------------------------------------+
  !+-------------------------------------------------------------------+
  !| This subroutine is used to combine all meanflow.****** and
  !| budget.****** to one HDF5 file.
  !+-------------------------------------------------------------------+
  subroutine meanflow_datacon
    !
    use commvardefine, only: ink,jnk,knk,imp,jmp,kmp,iop,jop,kop,      &
                             im,jm,km,rankmax
    use h5readwrite
    !
    logical :: lfex,lprint
    integer :: nsamples,nstep,nrank
    integer :: ni,nj,nk,ima,jma,kma,i,j,k,irk,jrk,krk
    character(len=8) :: mpirankname
    real(8) :: timesb,timese
    integer,allocatable :: SampStep(:)
    real(8),allocatable,dimension(:,:,:) :: rom,u1m,u2m,u3m,pm,tm,u11m,&
                                            u22m,u33m,u12m,u13m,u23m,  &
                                            ppm,ttm,tu1m,tu2m,tu3m,    &
                                            u111m,u222m,u333m,u112m,   &
                                            u113m,u122m,u133m,u223m,   &
                                            u233m,u123m,               &
                                            pu1m,pu2m,pu3m,u1rem,u2rem,&
                                            u3rem,sgmam11m,sgmam22m,   &
                                            sgmam33m,sgmam12m,sgmam13m,&
                                            sgmam23m,disspam,predilm,  &
                                            visdifm1,visdifm2,visdifm3
    !
    real(8),allocatable,dimension(:,:,:) :: rot,u1t,u2t,u3t,pt,tt,sgst
    real(8),allocatable,dimension(:,:,:) :: u11t,u22t,u33t,u12t,u13t,  &
                 u23t,ttt,ppt,tu1t,tu2t,tu3t,u111t,u222t,u333t,u112t,  &
                                 u113t,u122t,u133t,u223t,u233t,u123t
    real(8),allocatable,dimension(:,:,:) :: u1remt,u2remt,u3remt,pu1t, &
                                 pu2t,pu3t,sgmam11t,sgmam22t,sgmam33t, &
                              sgmam12t,sgmam13t,sgmam23t,disspat,predilt
    real(8), allocatable, dimension(:,:,:,:) :: visdift
    real(8), allocatable, dimension(:,:,:,:) :: du1m1,du2m2,du3m3,     &
                                                du1m1t,du2m2t,du3m3t
    !
    !
    !$ save rot,u1t,u2t,u3t,pt,tt,u11t,u22t,u33t,u12t,u13t,u23t,ttt,   &
    !$ ppt,tu1t,tu2t,tu3t,u111t,u222t,u333t,u112t,u113t,u122t,u133t,   &
    !$ u223t,u233t,u123t,u1remt,u2remt,u3remt,pu1t,pu2t,pu3t,          &
    !$ sgmam11t,sgmam22t,sgmam33t,sgmam12t,sgmam13t,sgmam23t,disspat,  &
    !$ predilt,visdift
    !
    !$OMP threadprivate(rot,u1t,u2t,u3t,pt,tt,u11t,u22t,u33t,u12t,     &
    !$OMP u13t,u23t,ttt,ppt,tu1t,tu2t,tu3t,u111t,u222t,u333t,u112t,    &
    !$OMP u113t,u122t,u133t,u223t,u233t,u123t,u1remt,u2remt,u3remt,    &
    !$OMP pu1t,pu2t,pu3t,sgmam11t,sgmam22t,sgmam33t,sgmam12t,sgmam13t, &
    !$OMP sgmam23t,disspat,predilt,visdift)
    !
    !+---------------------+
    !| processe meanflow   |
    !+---------------------+
    !
    lprint=.true.
    !
    allocate(rom(0:im,0:jm,0:km),u1m(0:im,0:jm,0:km),                  &
             u2m(0:im,0:jm,0:km),u3m(0:im,0:jm,0:km),                  &
             pm(0:im,0:jm,0:km),tm(0:im,0:jm,0:km)                     )
    allocate(u11m(0:im,0:jm,0:km),u22m(0:im,0:jm,0:km),                &
             u33m(0:im,0:jm,0:km),u12m(0:im,0:jm,0:km),                &
             u13m(0:im,0:jm,0:km),u23m(0:im,0:jm,0:km),                &
             ppm(0:im,0:jm,0:km),ttm(0:im,0:jm,0:km),                  &
             tu1m(0:im,0:jm,0:km),tu2m(0:im,0:jm,0:km),                &
             tu3m(0:im,0:jm,0:km)                                      )
    allocate(u111m(0:im,0:jm,0:km),u222m(0:im,0:jm,0:km),              &
             u333m(0:im,0:jm,0:km),u112m(0:im,0:jm,0:km),              &
             u113m(0:im,0:jm,0:km),u122m(0:im,0:jm,0:km),              &
             u133m(0:im,0:jm,0:km),u223m(0:im,0:jm,0:km),              &
             u233m(0:im,0:jm,0:km),u123m(0:im,0:jm,0:km)               )
    !
    ima=0
    jma=0
    kma=0
    ! 
    !$OMP parallel default(shared) private(nrank,i,j,k,irk,jrk,krk,    &
    !$OMP                              ni,nj,nk,ima,jma,kma,mpirankname)
    !$OMP do
    do nrank=0,rankmax
      !
      write(mpirankname,'(i7.7)')nrank
      !
      irk=ink(nrank)
      jrk=jnk(nrank)
      krk=knk(nrank)
      !
      ni=imp(nrank)
      nj=jmp(nrank)
      nk=kmp(nrank)
      !
      ima=iop(nrank)
      jma=jop(nrank)
      kma=kop(nrank)
      !
      allocate( rot(0:ni,0:nj,0:nk),u1t(0:ni,0:nj,0:nk),               &
                u2t(0:ni,0:nj,0:nk),u3t(0:ni,0:nj,0:nk),               &
                pt(0:ni,0:nj,0:nk),  tt(0:ni,0:nj,0:nk)                )
      allocate(u11t(0:ni,0:nj,0:nk),u22t(0:ni,0:nj,0:nk),              &
               u33t(0:ni,0:nj,0:nk),u12t(0:ni,0:nj,0:nk),              &
               u13t(0:ni,0:nj,0:nk),u23t(0:ni,0:nj,0:nk),              &
                ttt(0:ni,0:nj,0:nk), ppt(0:ni,0:nj,0:nk),              &
                tu1t(0:ni,0:nj,0:nk),tu2t(0:ni,0:nj,0:nk),             &
                tu3t(0:ni,0:nj,0:nk),                                  &
                u111t(0:ni,0:nj,0:nk),u222t(0:ni,0:nj,0:nk),           &
                u333t(0:ni,0:nj,0:nk),u112t(0:ni,0:nj,0:nk),           &
                u113t(0:ni,0:nj,0:nk),u122t(0:ni,0:nj,0:nk),           &
                u133t(0:ni,0:nj,0:nk),u223t(0:ni,0:nj,0:nk),           &
                u233t(0:ni,0:nj,0:nk),u123t(0:ni,0:nj,0:nk)            )
      !
      open(1000+nrank,file='Outdat/meanflow.'//mpirankname,            &
                       form='unformatted',status='old',action ='read', &
                               access='sequential')
      read(1000+nrank)nsamples,nstep
      read(1000+nrank)timesb,timese
      !
      if(lprint) then
        write(*,'(A)')'-------------------------------------------------'
        write(*,'(1x,2(1x,A8),2(1x,A14))')'nsamples','nstep','timesb','timese'
        write(*,'(1x,2(1x,I8),2(1x,F14.7))')nsamples,nstep,timesb,timese
        write(*,'(A)')'-------------------------------------------------'
        lprint=.false.
      endif
      !
      read(1000+nrank)rot,u1t,u2t,u3t,pt,tt
      read(1000+nrank)u11t,u22t,u33t,u12t,u13t,u23t,ttt,ppt,tu1t,tu2t,tu3t
      read(1000+nrank)u111t,u222t,u333t,u112t,u113t,u122t,u133t,u223t,u233t,u123t
      !
      if(.not. allocated(SampStep)) then
        allocate(SampStep(1:nsamples))
        read(1000+nrank)SampStep(1:nsamples)
      endif
      !
      close(1000+nrank)
      !
      print*,' >> file: meanflow.',mpirankname
      !
      do k=0,kmp(nrank)
      do j=0,jmp(nrank)
      do i=0,imp(nrank)
        !
        rom(i+ima,j+jma,k+kma)=rot(i,j,k)
        u1m(i+ima,j+jma,k+kma)=u1t(i,j,k)
        u2m(i+ima,j+jma,k+kma)=u2t(i,j,k)
        u3m(i+ima,j+jma,k+kma)=u3t(i,j,k)
         pm(i+ima,j+jma,k+kma)= pt(i,j,k)
         tm(i+ima,j+jma,k+kma)= tt(i,j,k)
        !
        u11m(i+ima,j+jma,k+kma)=u11t(i,j,k)
        u22m(i+ima,j+jma,k+kma)=u22t(i,j,k)
        u33m(i+ima,j+jma,k+kma)=u33t(i,j,k)
        u12m(i+ima,j+jma,k+kma)=u12t(i,j,k)
        u13m(i+ima,j+jma,k+kma)=u13t(i,j,k)
        u23m(i+ima,j+jma,k+kma)=u23t(i,j,k)
         ppm(i+ima,j+jma,k+kma)= ppt(i,j,k)
        ttm(i+ima,j+jma,k+kma)= ttt(i,j,k)
        tu1m(i+ima,j+jma,k+kma)=tu1t(i,j,k)
        tu2m(i+ima,j+jma,k+kma)=tu2t(i,j,k)
        tu3m(i+ima,j+jma,k+kma)=tu3t(i,j,k)
        !
        u111m(i+ima,j+jma,k+kma)=u111t(i,j,k)
        u222m(i+ima,j+jma,k+kma)=u222t(i,j,k)
        u333m(i+ima,j+jma,k+kma)=u333t(i,j,k)
        u112m(i+ima,j+jma,k+kma)=u112t(i,j,k)
        u113m(i+ima,j+jma,k+kma)=u113t(i,j,k)
        u122m(i+ima,j+jma,k+kma)=u122t(i,j,k)
        u133m(i+ima,j+jma,k+kma)=u133t(i,j,k)
        u223m(i+ima,j+jma,k+kma)=u223t(i,j,k)
        u233m(i+ima,j+jma,k+kma)=u233t(i,j,k)
        u123m(i+ima,j+jma,k+kma)=u123t(i,j,k)
        !
      end do
      end do
      end do
      !
      deallocate( rot,u1t,u2t,u3t,pt,tt )
      deallocate( u11t,u22t,u33t,u12t,u13t,u23t,ttt,ppt,tu1t,tu2t,tu3t,&
                  u111t,u222t,u333t,u112t,u113t,u122t,u133t,u223t,     &
                  u233t,u123t)
      !
    end do
    !$OMP end do
    !$OMP end parallel
    !
    inquire(file='Outdat/meanflow.h5',exist=lfex)
    if(lfex) call system('mv -v Outdat/meanflow.h5 Outdat/meanflow.h5.bak')
    !
    call H5WriteArray(nsamples,'nsamples','Outdat/meanflow.h5')
    call H5WriteArray(nstep,'nstep','Outdat/meanflow.h5')
    !
    call H5WriteArray(timesb,'timesb','Outdat/meanflow.h5')
    call H5WriteArray(timese,'timese','Outdat/meanflow.h5')
    !
    call H5WriteArray(rom,im,jm,km,'rom','Outdat/meanflow.h5')
    call H5WriteArray(u1m,im,jm,km,'u1m','Outdat/meanflow.h5')
    call H5WriteArray(u2m,im,jm,km,'u2m','Outdat/meanflow.h5')
    call H5WriteArray(u3m,im,jm,km,'u3m','Outdat/meanflow.h5')
    call H5WriteArray( pm,im,jm,km, 'pm','Outdat/meanflow.h5')
    call H5WriteArray( tm,im,jm,km, 'tm','Outdat/meanflow.h5')
    !
    inquire(file='Outdat/2ndsta.h5',exist=lfex)
    if(lfex) call system('mv -v Outdat/2ndsta.h5 Outdat/2ndsta.h5.bak')
    !
    call H5WriteArray(nsamples,'nsamples','Outdat/2ndsta.h5')
    call H5WriteArray(nstep,'nstep','Outdat/2ndsta.h5')
    call H5WriteArray(u11m,im,jm,km,'u11','Outdat/2ndsta.h5')
    call H5WriteArray(u22m,im,jm,km,'u22','Outdat/2ndsta.h5')
    call H5WriteArray(u33m,im,jm,km,'u33','Outdat/2ndsta.h5')
    call H5WriteArray(u12m,im,jm,km,'u12','Outdat/2ndsta.h5')
    call H5WriteArray(u13m,im,jm,km,'u13','Outdat/2ndsta.h5')
    call H5WriteArray(u23m,im,jm,km,'u23','Outdat/2ndsta.h5')
    call H5WriteArray( ppm,im,jm,km, 'pp','Outdat/2ndsta.h5')
    call H5WriteArray( ttm,im,jm,km, 'tt','Outdat/2ndsta.h5')
    call H5WriteArray(tu1m,im,jm,km,'tu1','Outdat/2ndsta.h5')
    call H5WriteArray(tu2m,im,jm,km,'tu2','Outdat/2ndsta.h5')
    call H5WriteArray(tu3m,im,jm,km,'tu3','Outdat/2ndsta.h5')
    !
    inquire(file='Outdat/3rdsta.h5',exist=lfex)
    if(lfex) call system('mv -v Outdat/3rdsta.h5 Outdat/3rdsta.h5.bak')
    !
    call H5WriteArray(nsamples,'nsamples','Outdat/3rdsta.h5')
    call H5WriteArray(nstep,'nstep','Outdat/3rdsta.h5')
    call H5WriteArray(u111m,im,jm,km,'u111','Outdat/3rdsta.h5')
    call H5WriteArray(u222m,im,jm,km,'u222','Outdat/3rdsta.h5')
    call H5WriteArray(u333m,im,jm,km,'u333','Outdat/3rdsta.h5')
    call H5WriteArray(u112m,im,jm,km,'u112','Outdat/3rdsta.h5')
    call H5WriteArray(u113m,im,jm,km,'u113','Outdat/3rdsta.h5')
    call H5WriteArray(u122m,im,jm,km,'u122','Outdat/3rdsta.h5')
    call H5WriteArray(u133m,im,jm,km,'u133','Outdat/3rdsta.h5')
    call H5WriteArray(u223m,im,jm,km,'u223','Outdat/3rdsta.h5')
    call H5WriteArray(u233m,im,jm,km,'u233','Outdat/3rdsta.h5')
    call H5WriteArray(u123m,im,jm,km,'u123','Outdat/3rdsta.h5')
    !
    call H5WriteArray(SampStep,nsamples-1,'SampStep','Outdat/meanflow.h5')
    !
    deallocate(rom,u1m,u2m,u3m,pm,tm,u11m,u22m,u33m,u12m,u13m,         &
               u23m,ppm,ttm,tu1m,tu2m,tu3m,u111m,u222m,u333m,u112m,    &
               u113m,u122m,u133m,u223m,u233m,u123m,SampStep            )
    !
    ! !+---------------------+
    ! !| processe budget     |
    ! !+---------------------+
    ! lprint=.true.
    ! !
    ! allocate(pu1m(0:im,0:jm,0:km),pu2m(0:im,0:jm,0:km),                &
    !          pu3m(0:im,0:jm,0:km),u1rem(0:im,0:jm,0:km),               &
    !          u2rem(0:im,0:jm,0:km),u3rem(0:im,0:jm,0:km),              &
    !          sgmam11m(0:im,0:jm,0:km),sgmam22m(0:im,0:jm,0:km),        &
    !          sgmam33m(0:im,0:jm,0:km),sgmam12m(0:im,0:jm,0:km),        &
    !          sgmam13m(0:im,0:jm,0:km),sgmam23m(0:im,0:jm,0:km),        &
    !          disspam(0:im,0:jm,0:km),predilm(0:im,0:jm,0:km),          &
    !          visdifm1(0:im,0:jm,0:km),visdifm2(0:im,0:jm,0:km),        &
    !          visdifm3(0:im,0:jm,0:km)                                 )
    ! !
    ! !$OMP parallel default(shared) private(nrank,i,j,k,irk,jrk,krk,    &
    ! !$OMP                              ni,nj,nk,ima,jma,kma,mpirankname)
    ! !$OMP do
    ! do nrank=0,rankmax
    !   !
    !   write(mpirankname,'(i7.7)')nrank
    !   !
    !   irk=ink(nrank)
    !   jrk=jnk(nrank)
    !   krk=knk(nrank)
    !   !
    !   ni=imp(nrank)
    !   nj=jmp(nrank)
    !   nk=kmp(nrank)
    !   !
    !   ima=iop(nrank)
    !   jma=jop(nrank)
    !   kma=kop(nrank)
    !   !
    !   allocate(u1remt(0:ni,0:nj,0:nk),u2remt(0:ni,0:nj,0:nk),          &
    !            u3remt(0:ni,0:nj,0:nk),  pu1t(0:ni,0:nj,0:nk),          &
    !              pu2t(0:ni,0:nj,0:nk),  pu3t(0:ni,0:nj,0:nk),          &
    !            sgmam11t(0:ni,0:nj,0:nk),sgmam22t(0:ni,0:nj,0:nk),      &
    !            sgmam33t(0:ni,0:nj,0:nk),sgmam12t(0:ni,0:nj,0:nk),      &
    !            sgmam13t(0:ni,0:nj,0:nk),sgmam23t(0:ni,0:nj,0:nk),      &
    !             disspat(0:ni,0:nj,0:nk), predilt(0:ni,0:nj,0:nk)       )
    !   allocate( visdift(3,0:ni,0:nj,0:nk)                              )
    !   !
    !   open(1000+nrank,file='Outdat/budget.'//mpirankname,              &
    !                    form='unformatted',status='old',action ='read', &
    !                            access='sequential')
    !   read(1000+nrank)nsamples,nstep
    !   read(1000+nrank)timesb,timese
    !   !
    !   if(lprint) then
    !     write(*,'(A)')'-------------------------------------------------'
    !     write(*,'(1x,2(1x,A8),2(1x,A14))')'nsamples','nstep','timesb','timese'
    !     write(*,'(1x,2(1x,I8),2(1x,F14.7))')nsamples,nstep,timesb,timese
    !     write(*,'(A)')'-------------------------------------------------'
    !     lprint=.false.
    !   endif
    !   !
    !   read(1000+nrank)u1remt,u2remt,u3remt,pu1t,pu2t,pu3t
    !   read(1000+nrank)sgmam11t,sgmam22t,sgmam33t,sgmam12t,sgmam13t,    &
    !                   sgmam23t,disspat,predilt
    !   read(1000+nrank)visdift
    !   !
    !   if(.not. allocated(SampStep)) then
    !     allocate(SampStep(1:nsamples))
    !     read(1000+nrank)SampStep(1:nsamples)
    !   endif
    !   !
    !   close(1000+nrank)
    !   !
    !   !
    !   print*,' >> file: budget.',mpirankname
    !   !
    !   do k=0,kmp(nrank)
    !   do j=0,jmp(nrank)
    !   do i=0,imp(nrank)
    !     !
    !     u1rem(i+ima,j+jma,k+kma)=u1remt(i,j,k)
    !     u2rem(i+ima,j+jma,k+kma)=u2remt(i,j,k)
    !     u3rem(i+ima,j+jma,k+kma)=u3remt(i,j,k)
    !      pu1m(i+ima,j+jma,k+kma)=pu1t(i,j,k)
    !      pu2m(i+ima,j+jma,k+kma)=pu2t(i,j,k)
    !      pu3m(i+ima,j+jma,k+kma)=pu3t(i,j,k)
    !     !
    !     sgmam11m(i+ima,j+jma,k+kma)=sgmam11t(i,j,k)
    !     sgmam22m(i+ima,j+jma,k+kma)=sgmam22t(i,j,k)
    !     sgmam33m(i+ima,j+jma,k+kma)=sgmam33t(i,j,k)
    !     sgmam12m(i+ima,j+jma,k+kma)=sgmam12t(i,j,k)
    !     sgmam13m(i+ima,j+jma,k+kma)=sgmam13t(i,j,k)
    !     sgmam23m(i+ima,j+jma,k+kma)=sgmam23t(i,j,k)
    !      disspam(i+ima,j+jma,k+kma)= disspat(i,j,k)
    !      predilm(i+ima,j+jma,k+kma)= predilt(i,j,k)
    !     !
    !     visdifm1(i+ima,j+jma,k+kma)=visdift(1,i,j,k)
    !     visdifm2(i+ima,j+jma,k+kma)=visdift(2,i,j,k)
    !     visdifm3(i+ima,j+jma,k+kma)=visdift(3,i,j,k)
    !     !
    !   end do
    !   end do
    !   end do
    !   !
    !   deallocate( u1remt,u2remt,u3remt,pu1t,pu2t,pu3t,                 &
    !               sgmam11t,sgmam22t,sgmam33t,sgmam12t,sgmam13t,        &
    !               sgmam23t,disspat,predilt                             )
    !   deallocate(visdift)
    !   !
    ! end do
    ! !$OMP end do
    ! !$OMP end parallel
    ! !
    ! inquire(file='Outdat/budget.h5',exist=lfex)
    ! if(lfex) call system('mv -v Outdat/budget.h5 Outdat/budget.h5.bak')
    ! !
    ! call H5WriteArray(nsamples,'nsamples','Outdat/budget.h5')
    ! call H5WriteArray(nstep,'nstep','Outdat/budget.h5')
    ! !
    ! call H5WriteArray(timesb,'timesb','Outdat/budget.h5')
    ! call H5WriteArray(timese,'timese','Outdat/budget.h5')
    ! !
    ! call H5WriteArray(u1rem,im,jm,km,'u1rem','Outdat/budget.h5')
    ! call H5WriteArray(u2rem,im,jm,km,'u2rem','Outdat/budget.h5')
    ! call H5WriteArray(u3rem,im,jm,km,'u3rem','Outdat/budget.h5')
    ! call H5WriteArray(pu1m,im,jm,km,'pu1','Outdat/budget.h5')
    ! call H5WriteArray(pu2m,im,jm,km,'pu2','Outdat/budget.h5')
    ! call H5WriteArray(pu3m,im,jm,km,'pu3','Outdat/budget.h5')
    ! !
    ! call H5WriteArray(sgmam11m,im,jm,km,'sgmam11','Outdat/budget.h5')
    ! call H5WriteArray(sgmam22m,im,jm,km,'sgmam22','Outdat/budget.h5')
    ! call H5WriteArray(sgmam33m,im,jm,km,'sgmam33','Outdat/budget.h5')
    ! call H5WriteArray(sgmam12m,im,jm,km,'sgmam12','Outdat/budget.h5')
    ! call H5WriteArray(sgmam13m,im,jm,km,'sgmam13','Outdat/budget.h5')
    ! call H5WriteArray(sgmam23m,im,jm,km,'sgmam23','Outdat/budget.h5')
    ! !
    ! call H5WriteArray(disspam,im,jm,km,'disspa','Outdat/budget.h5')
    ! call H5WriteArray(predilm,im,jm,km,'predil','Outdat/budget.h5')
    ! call H5WriteArray(visdifm1,im,jm,km,'visdif1','Outdat/budget.h5')
    ! call H5WriteArray(visdifm2,im,jm,km,'visdif2','Outdat/budget.h5')
    ! call H5WriteArray(visdifm3,im,jm,km,'visdif3','Outdat/budget.h5')
    ! !
    ! call H5WriteArray(SampStep,nsamples-1,'SampStep','Outdat/budget.h5')
    ! !
    ! deallocate(pu1m,pu2m,pu3m,u1rem,u2rem,u3rem,disspam,predilm,       &
    !            sgmam11m,sgmam22m,sgmam33m,sgmam12m,sgmam13m,sgmam23m,  &
    !            visdifm1,visdifm2,visdifm3,SampStep   )
    !
  end subroutine meanflow_datacon
  !+-------------------------------------------------------------------+
  !|the end of the subroutine meanflow_datacon.                        |
  !+-------------------------------------------------------------------+
  subroutine budget_datacon
    !
    use commvardefine, only: ink,jnk,knk,imp,jmp,kmp,iop,jop,kop,      &
                             im,jm,km,rankmax
    use h5readwrite
    !
    logical :: lfex,lprint
    integer :: nsamples,nstep,nrank
    integer :: ni,nj,nk,ima,jma,kma,i,j,k,irk,jrk,krk
    character(len=8) :: mpirankname
    real(8) :: timesb,timese
    integer,allocatable :: SampStep(:)
    real(8),allocatable,dimension(:,:,:) :: u1rem,u2rem,u3rem,         &
                                            pu1,pu2,pu3,               &
                                            sgmam11,sgmam22,sgmam33,   &
                                            sgmam12,sgmam13,sgmam23,   &
                                            u1rem_t,u2rem_t,u3rem_t,         &
                                            pu1_t,pu2_t,pu3_t,               &
                                            sgmam11_t,sgmam22_t,sgmam33_t,   &
                                            sgmam12_t,sgmam13_t,sgmam23_t
    real(8),allocatable,dimension(:,:,:,:) :: sigmau,sigmadu,pdu,      &
                                              sigmau_t,sigmadu_t,pdu_t
    !
    !$ save u1rem_t,u2rem_t,u3rem_t,pu1_t,pu2_t,pu3_t,                 &
    !$      sgmam11_t,sgmam22_t,sgmam33_t,sgmam12_t,sgmam13_t,         &
    !$      sgmam23_t,sigmau_t,sigmadu_t,pdu_t
    ! 
    !$OMP threadprivate(u1rem_t,u2rem_t,u3rem_t,pu1_t,pu2_t,pu3_t,     &
    !$OMP       sgmam11_t,sgmam22_t,sgmam33_t,sgmam12_t,sgmam13_t,     &
    !$OMP       sgmam23_t,sigmau_t,sigmadu_t,pdu_t )
    !
    !+---------------------+
    !| processe budget     |
    !+---------------------+
    lprint=.true.
    !
    allocate( u1rem(0:im,0:jm,0:km),u2rem(0:im,0:jm,0:km),             &
              u3rem(0:im,0:jm,0:km),pu1(0:im,0:jm,0:km),               &
              pu2(0:im,0:jm,0:km),  pu3(0:im,0:jm,0:km),               &
              sgmam11(0:im,0:jm,0:km),sgmam22(0:im,0:jm,0:km),         &
              sgmam33(0:im,0:jm,0:km),sgmam12(0:im,0:jm,0:km),         &
              sgmam13(0:im,0:jm,0:km),sgmam23(0:im,0:jm,0:km)  )
    allocate( pdu(1:6,0:im,0:jm,0:km),sigmau(18,0:im,0:jm,0:km),       &
              sigmadu(9,0:im,0:jm,0:km) )
    !
    !$OMP parallel default(shared) private(nrank,i,j,k,irk,jrk,krk,    &
    !$OMP                              ni,nj,nk,ima,jma,kma,mpirankname)
    !$OMP do
    do nrank=0,rankmax
      !
      write(mpirankname,'(i7.7)')nrank
      !
      irk=ink(nrank)
      jrk=jnk(nrank)
      krk=knk(nrank)
      !
      ni=imp(nrank)
      nj=jmp(nrank)
      nk=kmp(nrank)
      !
      ima=iop(nrank)
      jma=jop(nrank)
      kma=kop(nrank)
      !
      allocate( u1rem_t(0:ni,0:nj,0:nk),u2rem_t(0:ni,0:nj,0:nk),       &
              u3rem_t(0:ni,0:nj,0:nk),pu1_t(0:ni,0:nj,0:nk),           &
              pu2_t(0:ni,0:nj,0:nk),  pu3_t(0:ni,0:nj,0:nk),           &
              sgmam11_t(0:ni,0:nj,0:nk),sgmam22_t(0:ni,0:nj,0:nk),     &
              sgmam33_t(0:ni,0:nj,0:nk),sgmam12_t(0:ni,0:nj,0:nk),     &
              sgmam13_t(0:ni,0:nj,0:nk),sgmam23_t(0:ni,0:nj,0:nk) )
      allocate( pdu_t(1:6,0:ni,0:nj,0:nk),sigmau_t(18,0:ni,0:nj,0:nk), &
                sigmadu_t(9,0:ni,0:nj,0:nk) )
      !
      open(1000+nrank,file='Outdat/budget.'//mpirankname,              &
                       form='unformatted',status='old',action ='read', &
                               access='sequential')
      read(1000+nrank)nsamples,nstep
      read(1000+nrank)timesb,timese
      !
      if(lprint) then
        write(*,'(A)')'-------------------------------------------------'
        write(*,'(1x,2(1x,A8),2(1x,A14))')'nsamples','nstep','timesb','timese'
        write(*,'(1x,2(1x,I8),2(1x,F14.7))')nsamples,nstep,timesb,timese
        write(*,'(A)')'-------------------------------------------------'
        lprint=.false.
      endif
      !
      read(1000+nrank)u1rem_t,u2rem_t,u3rem_t
      read(1000+nrank)pu1_t,pu2_t,pu3_t
      read(1000+nrank)pdu_t
      read(1000+nrank)sgmam11_t,sgmam22_t,sgmam33_t,                 &
                        sgmam12_t,sgmam13_t,sgmam23_t
      read(1000+nrank)sigmau_t
      read(1000+nrank)sigmadu_t
      !
      if(.not. allocated(SampStep)) then
        allocate(SampStep(1:nsamples))
        read(1000+nrank)SampStep(1:nsamples)
      endif
      !
      close(1000+nrank)
      !
      !
      print*,' >> file: budget.',mpirankname
      !
      do k=0,kmp(nrank)
      do j=0,jmp(nrank)
      do i=0,imp(nrank)
        !
        u1rem(i+ima,j+jma,k+kma)    = u1rem_t(i,j,k)
        u2rem(i+ima,j+jma,k+kma)    = u2rem_t(i,j,k)
        u3rem(i+ima,j+jma,k+kma)    = u3rem_t(i,j,k)
        pu1(i+ima,j+jma,k+kma)      = pu1_t(i,j,k)
        pu2(i+ima,j+jma,k+kma)      = pu2_t(i,j,k)
        pu3(i+ima,j+jma,k+kma)      = pu3_t(i,j,k)
        pdu(:,i+ima,j+jma,k+kma)    = pdu_t(:,i,j,k)
        sgmam11(i+ima,j+jma,k+kma)  = sgmam11_t(i,j,k)
        sgmam22(i+ima,j+jma,k+kma)  = sgmam22_t(i,j,k)
        sgmam33(i+ima,j+jma,k+kma)  = sgmam33_t(i,j,k)
        sgmam12(i+ima,j+jma,k+kma)  = sgmam12_t(i,j,k)
        sgmam13(i+ima,j+jma,k+kma)  = sgmam13_t(i,j,k)
        sgmam23(i+ima,j+jma,k+kma)  = sgmam23_t(i,j,k)
        sigmau(:,i+ima,j+jma,k+kma) = sigmau_t(:,i,j,k)
        sigmadu(:,i+ima,j+jma,k+kma)= sigmadu_t(:,i,j,k)
        !
      end do
      end do
      end do
      !
      deallocate( u1rem_t,u2rem_t,u3rem_t,pu1_t,pu2_t,pu3_t,pdu_t,     &
                  sgmam11_t,sgmam22_t,sgmam33_t,sgmam12_t,sgmam13_t,   &
                  sgmam23_t,sigmau_t,sigmadu_t )
      !
    end do
    !$OMP end do
    !$OMP end parallel
    !
    inquire(file='Outdat/budget.h5',exist=lfex)
    if(lfex) call system('mv -v Outdat/budget.h5 Outdat/budget.h5.bak')
    !
    call H5WriteArray(nsamples,'nsamples','Outdat/budget.h5')
    call H5WriteArray(nstep,'nstep','Outdat/budget.h5')
    !
    call H5WriteArray(timesb,'timesb','Outdat/budget.h5')
    call H5WriteArray(timese,'timese','Outdat/budget.h5')
    !
    call H5WriteArray(u1rem,im,jm,km,'u1rem','Outdat/budget.h5')
    call H5WriteArray(u2rem,im,jm,km,'u2rem','Outdat/budget.h5')
    call H5WriteArray(u3rem,im,jm,km,'u3rem','Outdat/budget.h5')
    !
    call H5WriteArray(pu1,im,jm,km,'pu1','Outdat/budget.h5')
    call H5WriteArray(pu2,im,jm,km,'pu2','Outdat/budget.h5')
    call H5WriteArray(pu3,im,jm,km,'pu3','Outdat/budget.h5')
    !
    call H5WriteArray(sgmam11,im,jm,km,'sgmam11','Outdat/budget.h5')
    call H5WriteArray(sgmam22,im,jm,km,'sgmam22','Outdat/budget.h5')
    call H5WriteArray(sgmam33,im,jm,km,'sgmam33','Outdat/budget.h5')
    call H5WriteArray(sgmam12,im,jm,km,'sgmam12','Outdat/budget.h5')
    call H5WriteArray(sgmam13,im,jm,km,'sgmam13','Outdat/budget.h5')
    call H5WriteArray(sgmam23,im,jm,km,'sgmam23','Outdat/budget.h5')
    !
    call H5WriteArray(pdu(1,:,:,:),im,jm,km,'pdu1','Outdat/budget.h5')
    call H5WriteArray(pdu(2,:,:,:),im,jm,km,'pdu2','Outdat/budget.h5')
    call H5WriteArray(pdu(3,:,:,:),im,jm,km,'pdu3','Outdat/budget.h5')
    call H5WriteArray(pdu(4,:,:,:),im,jm,km,'pdu4','Outdat/budget.h5')
    call H5WriteArray(pdu(5,:,:,:),im,jm,km,'pdu5','Outdat/budget.h5')
    call H5WriteArray(pdu(6,:,:,:),im,jm,km,'pdu6','Outdat/budget.h5')
    !
    call H5WriteArray(sigmau(1,:,:,:) ,im,jm,km,'sigmau1' ,'Outdat/budget.h5')
    call H5WriteArray(sigmau(2,:,:,:) ,im,jm,km,'sigmau2' ,'Outdat/budget.h5')
    call H5WriteArray(sigmau(3,:,:,:) ,im,jm,km,'sigmau3' ,'Outdat/budget.h5')
    call H5WriteArray(sigmau(4,:,:,:) ,im,jm,km,'sigmau4' ,'Outdat/budget.h5')
    call H5WriteArray(sigmau(5,:,:,:) ,im,jm,km,'sigmau5' ,'Outdat/budget.h5')
    call H5WriteArray(sigmau(6,:,:,:) ,im,jm,km,'sigmau6' ,'Outdat/budget.h5')
    call H5WriteArray(sigmau(7,:,:,:) ,im,jm,km,'sigmau7' ,'Outdat/budget.h5')
    call H5WriteArray(sigmau(8,:,:,:) ,im,jm,km,'sigmau8' ,'Outdat/budget.h5')
    call H5WriteArray(sigmau(9,:,:,:) ,im,jm,km,'sigmau9' ,'Outdat/budget.h5')
    call H5WriteArray(sigmau(10,:,:,:),im,jm,km,'sigmau10','Outdat/budget.h5')
    call H5WriteArray(sigmau(11,:,:,:),im,jm,km,'sigmau11','Outdat/budget.h5')
    call H5WriteArray(sigmau(12,:,:,:),im,jm,km,'sigmau12','Outdat/budget.h5')
    call H5WriteArray(sigmau(13,:,:,:),im,jm,km,'sigmau13','Outdat/budget.h5')
    call H5WriteArray(sigmau(14,:,:,:),im,jm,km,'sigmau14','Outdat/budget.h5')
    call H5WriteArray(sigmau(15,:,:,:),im,jm,km,'sigmau15','Outdat/budget.h5')
    call H5WriteArray(sigmau(16,:,:,:),im,jm,km,'sigmau16','Outdat/budget.h5')
    call H5WriteArray(sigmau(17,:,:,:),im,jm,km,'sigmau17','Outdat/budget.h5')
    call H5WriteArray(sigmau(18,:,:,:),im,jm,km,'sigmau18','Outdat/budget.h5')
    !
    call H5WriteArray(sigmadu(1,:,:,:),im,jm,km,'sigmadu1','Outdat/budget.h5')
    call H5WriteArray(sigmadu(2,:,:,:),im,jm,km,'sigmadu2','Outdat/budget.h5')
    call H5WriteArray(sigmadu(3,:,:,:),im,jm,km,'sigmadu3','Outdat/budget.h5')
    call H5WriteArray(sigmadu(4,:,:,:),im,jm,km,'sigmadu4','Outdat/budget.h5')
    call H5WriteArray(sigmadu(5,:,:,:),im,jm,km,'sigmadu5','Outdat/budget.h5')
    call H5WriteArray(sigmadu(6,:,:,:),im,jm,km,'sigmadu6','Outdat/budget.h5')
    call H5WriteArray(sigmadu(7,:,:,:),im,jm,km,'sigmadu7','Outdat/budget.h5')
    call H5WriteArray(sigmadu(8,:,:,:),im,jm,km,'sigmadu8','Outdat/budget.h5')
    call H5WriteArray(sigmadu(9,:,:,:),im,jm,km,'sigmadu9','Outdat/budget.h5')
    !
    !
  end subroutine budget_datacon
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to interpolate islice data                |
  !+-------------------------------------------------------------------+
  subroutine isliceinterp(nsta,nend,database)
    !
    use interpolation
    use commvardefine,only: im,jm,km,gridfile,const2
    !
    use h5readwrite
    use writetec
    !
    integer,intent(in) :: nsta,nend
    character(len=*),intent(in) :: database
    ! 
    integer :: n,i,j,k,k1,j1
    integer :: im2,jm2,km2,nmax
    character(len=5) :: fname
    character(len=255) :: filename
    !
    integer :: dims(3)
    real(8),allocatable,dimension(:,:) :: x,y,z
    real(8),allocatable,dimension(:,:)  :: ro,u1,u2,u3,t,p
    real(8),allocatable,dimension(:) :: rom,u1m,u2m,tm,rom2,u1m2,u2m2,tm2
    real(8),allocatable,dimension(:,:) :: royz2,u1yz2,u2yz2,u3yz2,tyz2,pyz2
    real(8),allocatable,dimension(:,:) :: x2,y2,z2
    real(8),allocatable,dimension(:) :: y1
    real(8) :: time,time0,timemax
    real(8) :: delta,dstar,theter,utaw,zmax,var1
    !
    nmax=3787
    !
    dims=h5_getdimensio('x',database//'datin/grid.h5')
    jm=dims(1)-1
    km=dims(2)-1
    print*,' ** dimension of the input data: ',jm,'x',km
    !
    print*,' ** reading mean profile data... '
    allocate(rom(0:jm),u1m(0:jm),u2m(0:jm),tm(0:jm))
    open(12,file=database//'Results/meanprofile.dat')
    read(12,*)
    read(12,*)delta,dstar,theter,utaw
    read(12,*)
    do j=0,jm
      read(12,*)var1,rom(j),u1m(j),u2m(j),tm(j)
    enddo
    close(12)
    print*,' << ',database,'meanprofile.dat'
    !
    ! call H5ReadArray(rom,jm,'ro',database//'Results/meanprofile.h5')
    ! call H5ReadArray(u1m,jm,'u1',database//'Results/meanprofile.h5')
    ! call H5ReadArray(u2m,jm,'u2',database//'Results/meanprofile.h5')
    ! call H5ReadArray( tm,jm,'t',database//'Results/meanprofile.h5')
    ! !
    filename=database//'islice/islice00000.h5'
    call H5ReadArray(time0,'time',trim(filename))
    print*,' ** time0 : ',time0
    !
    write(fname,'(i5.5)')nmax
    filename=database//'islice/islice'//fname//'.h5'
    call H5ReadArray(timemax,'time',trim(filename))
    print*,' ** timemax : ',timemax
    !
    print*,' ** reading input grid ... '
    !
    allocate(y(0:jm,0:km),z(0:jm,0:km))
    ! call h5_read2dfrom3d(x2,im,jm,km,'x','datin/grid.h5',islice=0)
    call H5ReadArray(y,jm,km,'y',database//'datin/grid.h5')
    call H5ReadArray(z,jm,km,'z',database//'datin/grid.h5')
    !
    zmax=z(0,km)
    print*,' ** zmax=',zmax
    !
    dims=h5_getdimensio('x','datin/grid.h5')
    print*,dims-1
    !
    im2=dims(1)-1
    jm2=dims(2)-1
    km2=dims(3)-1
    print*,' ** dimension of grid: ',jm2,'x',km2
    !
    ! !
    ! print*,' ** reading  grid ... '
    allocate(y2(0:jm2,0:km2),z2(0:jm2,0:km2))
    ! ! call h5_read2dfrom3d(x,im,jm,km,'x','datin/grid.h5',islice=0)
    call H5ReadSubset(y2,im2,jm2,km2,'y','datin/grid.h5',islice=0)
    call H5ReadSubset(z2,im2,jm2,km2,'z','datin/grid.h5',islice=0)
    !
    allocate(rom2(0:jm2),u1m2(0:jm2),u2m2(0:jm2),tm2(0:jm2))
    !
    call regularlinearinterp(y(:,0),rom,u1m,u2m,tm,jm,               &
                             y2(:,0),rom2,u1m2,u2m2,tm2,jm2  )
    !
    open(16,file='inlet.prof')    
    write(16,"(A26)")'# parameters of inlet flow'
    write(16,"(4(1X,A14))")'delta','delta*','theter','utaw'
    write(16,"(4(1X,E14.7E2))")delta,dstar,theter,utaw
    write(16,"(4(1X,A14))")'ro','u1','u2','t'
    write(16,"(4(1X,E14.7E2))")(rom2(j),u1m2(j),u2m2(j),tm2(j),j=0,jm2)
    close(16)
    print*,' << inlet.prof ... done!'
    ! !
    ! open(16,file='meanprofile.dat')  
    ! write(16,"(5(1X,A14))")'y','ro','u1','u2','t'
    ! write(16,"(5(1X,E14.7E2))")(y2(j,0),rom(j),u1m(j),u2m(j),tm(j),j=0,jm2)
    ! close(16)
    ! print*,' << meanprofile.dat ... done!'
    ! !
    ! open(16,file='inlet.dat')  
    ! write(16,"(5(1X,A14))")'y','ro','u1','u2','t'
    ! write(16,"(5(1X,E14.7E2))")(y(j,0),rom2(j),u1m2(j),u2m2(j),tm2(j),j=0,jm)
    ! close(16)
    ! print*,' << inlet.dat ... done!'
    !
    ! allocate(y2(0:jm2,0:km2),z2(0:jm2,0:km2))
    ! call h5_read2dfrom3d(y2,im2,jm2,km2,'y','datin/grid.h5',islice=0)
    ! call h5_read2dfrom3d(z2,im2,jm2,km2,'z','datin/grid.h5',islice=0)
    ! !
    ! print*,' ** reading mean profile data... '
    ! allocate(rom(0:jm2),u1m(0:jm2),u2m(0:jm2),tm(0:jm2))
    ! open(12,file='datin/inlet.prof')    
    ! read(12,*)
    ! read(12,*)
    ! read(12,*)
    ! read(12,*)
    ! read(12,*)(rom(j),u1m(j),u2m(j),tm(j),j=0,jm2)
    ! close(12)
    ! print*,' >> inlet.prof ... done!'
    ! !
    ! allocate(y1(0:jm))
    ! open(12,file='gridy.dat')
    ! do j=0,jm
    !   read(12,*)j1,y1(j)
    ! end do
    ! close(12)
    ! print*,' >> gridy.dat ... done.'
    ! !
    ! do k=0,km
    ! do j=0,jm
    !     y(j,k)=y1(j)
    !     z(j,k)=z2(0,k)
    ! enddo
    ! enddo
    ! !
    ! allocate(rom2(0:jm),u1m2(0:jm),u2m2(0:jm),tm2(0:jm))
    ! !
    ! call regularlinearinterp(y2(:,0),rom,u1m,u2m,tm,jm2,               &
    !                          y(:,0),rom2,u1m2,u2m2,tm2,jm  )
    ! !
    ! open(16,file='inlet.prof')    
    ! write(16,"(A26)")'# parameters of inlet flow'
    ! write(16,"(4(1X,A14))")'delta','delta*','theter','utaw'
    ! write(16,"(4(1X,A14))")'ro','u1','u2','t'
    ! write(16,"(4(1X,E14.7E2))")(rom2(j),u1m2(j),u2m2(j),tm2(j),j=0,jm)
    ! close(16)
    ! print*,' << inlet.prof ... done!'
    ! !
    ! open(16,file='meanprofile.dat')  
    ! write(16,"(5(1X,A14))")'y','ro','u1','u2','t'
    ! write(16,"(5(1X,E14.7E2))")(y2(j,0),rom(j),u1m(j),u2m(j),tm(j),j=0,jm2)
    ! close(16)
    ! print*,' << meanprofile.dat ... done!'
    !
    allocate(ro(0:jm,0:km),u1(0:jm,0:km),u2(0:jm,0:km),u3(0:jm,0:km),t(0:jm,0:km))
    allocate(royz2(0:jm2,0:km2),u1yz2(0:jm2,0:km2),u2yz2(0:jm2,0:km2), &
             u3yz2(0:jm2,0:km2), tyz2(0:jm2,0:km2), pyz2(0:jm2,0:km2)  )
    !
    do n=nsta,nend
      !
      if(n>nmax) then
        write(fname,'(i5.5)')n-nmax
      else
        write(fname,'(i5.5)')n
      endif
      !
      filename=database//'islice/islice'//fname//'.h5'
      !
      call H5ReadArray(time,'time',  trim(filename),explicit=.false.)
      call H5ReadArray(ro,jm,km,'ro',trim(filename),explicit=.false.)
      call H5ReadArray(u1,jm,km,'u1',trim(filename),explicit=.false.)
      call H5ReadArray(u2,jm,km,'u2',trim(filename),explicit=.false.)
      call H5ReadArray(u3,jm,km,'u3',trim(filename),explicit=.false.)
      call H5ReadArray( t,jm,km, 't',trim(filename),explicit=.false.)
      !
      if(n>nmax) then
        time=time+timemax-time0
      endif
      !
      do j=0,jm
        ro(j,:)=ro(j,:)-rom(j)
        u1(j,:)=u1(j,:)-u1m(j)
        u2(j,:)=u2(j,:)-u2m(j)
         t(j,:)= t(j,:)- tm(j)
      enddo
      !
      write(*,'(1A1,2(A),$)')char(13),'  >> ',trim(filename)
      !
      call regularlinearinterp(y,z,ro,u1,u2,u3,t,jm,km,                   &
                               y2,z2,royz2,u1yz2,u2yz2,u3yz2,tyz2,jm2,km2,&
                               homoy=.true.,progress=.false.)
      !
      write(fname,'(i5.5)')n
      filename='inflow/islice'//fname//'.h5'
      !
      time=time-time0
      !
      call H5WriteArray(time,'time',trim(filename),explicit=.false.)
      call H5WriteArray(royz2,jm2,km2,'ro',trim(filename),explicit=.false.)
      call H5WriteArray(u1yz2,jm2,km2,'u1',trim(filename),explicit=.false.)
      call H5WriteArray(u2yz2,jm2,km2,'u2',trim(filename),explicit=.false.)
      call H5WriteArray(u3yz2,jm2,km2,'u3',trim(filename),explicit=.false.)
      call H5WriteArray( tyz2,jm2,km2,'t',trim(filename),explicit=.false.)
      !
      write(*,'(1A1,3(A),I0,A,$)')char(13),'  << ',trim(filename),   &
                                   ' ... ',(n-nsta)*100/(nend-nsta),' %'
      !
    end do
    !
    call writetecbin('tecro2.plt',y2,'y',z2,'z',royz2,'ro',u1yz2,'u',  &
                                  u2yz2,'v',u3yz2,'w',tyz2,'t',jm2,km2)
    do j=0,jm2
      royz2(j,:)=royz2(j,:)+rom2(j)
      u1yz2(j,:)=u1yz2(j,:)+u1m2(j)
      u2yz2(j,:)=u2yz2(j,:)+u2m2(j)
       tyz2(j,:)= tyz2(j,:)+ tm2(j)
      !
    enddo
    !
    pyz2=royz2*tyz2/const2
    !
    call writetecbin('tecinflow.plt',y2,'y',z2,'z',royz2,'ro',u1yz2,'u',  &
                          u2yz2,'v',u3yz2,'w',tyz2,'t',pyz2,'p',jm2,km2)
      
  end subroutine isliceinterp
  !+-------------------------------------------------------------------+
  !| The end of the subroutine isliceinterp.                           |
  !+-------------------------------------------------------------------+
  !
  subroutine isliceinterp_chan(nsta,nend,database)
    !
    use interpolation
    use commvardefine,only: im,jm,km,gridfile,const2
    !
    use h5readwrite
    use writetec
    !
    integer,intent(in) :: nsta,nend
    character(len=*),intent(in) :: database
    ! 
    integer :: n,n2,i,j,k,k1,j1,j2,jmm
    integer :: im2,jm2,km2,nmax
    character(len=5) :: fname
    character(len=255) :: filename
    !
    integer :: dims(3)
    real(8),allocatable,dimension(:,:) :: x,y,z
    real(8),allocatable,dimension(:,:)  :: ro,u1,u2,u3,t,p
    real(8),allocatable,dimension(:) :: rom,u1m,u2m,tm,rom2,u1m2,u2m2,tm2
    real(8),allocatable,dimension(:,:) :: royz2,u1yz2,u2yz2,u3yz2,tyz2,pyz2
    real(8),allocatable,dimension(:,:) :: x2,y2,z2
    real(8),allocatable,dimension(:) :: y1
    real(8),allocatable,dimension(:,:,:) :: vary
    real(8) :: time,time0,timemax,scale,hstep,yhight
    real(8) :: delta,dstar,theter,utaw,ymax,zmax,var1
    real(8) :: l_ref,u_ref,ro_ref,t_ref
    real(8) :: time_ref
    !
    nmax=7513
    !
    l_ref=0.001466d0/1.265213d0
    u_ref=950.981d0
    ro_ref=0.174395d0
    t_ref=1000.d0
    time_ref=l_ref/u_ref
    !
    scale=0.001466d0/1.265213d0
    hstep=3.05d0/1000.d0
    yhight=14.66d0/1000.d0
    ! scale=1.d0
    ! hstep=0.d0
    ! yhight=10.d0
    !
    dims=h5_getdimensio('x',database//'datin/grid.h5')
    jm=(dims(1)-1)*2+1
    km=dims(2)-1
    print*,' ** dimension of the input data: ',jm,'x',km
    !
    jmm=(jm-1)/2
    !
    print*,' ** reading mean profile data... '
    allocate(rom(0:jmm),u1m(0:jmm),u2m(0:jmm),tm(0:jmm))
    call H5ReadArray(rom(0:jmm),jmm,'ro',database//'Results/meanprofile.h5')
    call H5ReadArray(u1m(0:jmm),jmm,'u1',database//'Results/meanprofile.h5')
    call H5ReadArray(u2m(0:jmm),jmm,'u2',database//'Results/meanprofile.h5')
    call H5ReadArray( tm(0:jmm),jmm,'t',database//'Results/meanprofile.h5')
    ! !
    filename=database//'islice/islice00000.h5'
    call H5ReadArray(time0,'time',trim(filename))
    print*,' ** time0 : ',time0
    !
    write(fname,'(i5.5)')nmax
    filename=database//'islice/islice'//fname//'.h5'
    call H5ReadArray(timemax,'time',trim(filename))
    print*,' ** timemax : ',timemax
    !
    print*,' ** reading input grid ... '
    !
    allocate(y(0:jm,0:km),z(0:jm,0:km))
    ! call h5_read2dfrom3d(x2,im,jm,km,'x','datin/grid.h5',islice=0)
    call H5ReadArray(y(0:jmm,:),jmm,km,'y',database//'datin/grid.h5')
    call H5ReadArray(z(0:jmm,:),jmm,km,'z',database//'datin/grid.h5')
    !
    y=y*l_ref
    z=z*l_ref
    !
    do j=jm,jmm+1,-1
      ! y(j,:)=14.66d0-y(jm-j,:)
      y(j,:)=yhight-y(jm-j,:)
      z(j,:)=z(jm-j,:)
    enddo
    !
    ymax=y(jm,km)
    zmax=z(0,km)
    print*,' ** zmax=',zmax
    print*,' ** ymax=',ymax
    print*,' **            y0        Δy1         Δymid-1        ymid         Δymid+1          Δymax   ymax'
    print*,y(0,0),y(1,0)-y(0,0),y(jmm,0)-y(jmm-1,0),y(jmm,0),y(jmm+1,0)-y(jmm,0),y(jm,0)-y(jm-1,0),y(jm,0)
    !
    y=y+hstep
    !
    open(12,file=database//'Results/meanprofile.dat')
    read(12,*)
    read(12,*)delta,dstar,theter,utaw
    close(12)
    print*,' << ',database,'meanprofile.dat'
    !
    delta=delta*l_ref
    !
    do j=0,jm
      if(y(j,0)>2.5d0*delta+hstep) then
        j2=j
        exit
      endif
    enddo
    print*,' ** j location of 2 bl thickness:',j2
    !
    dims=h5_getdimensio('x','datin/grid.h5')
    print*,dims-1
    !
    im2=dims(1)-1
    jm2=dims(2)-1
    km2=dims(3)-1
    print*,' ** dimension of grid: ',jm2,'x',km2
    !
    ! print*,' ** reading  grid ... '
    allocate(y2(0:jm2,0:km2),z2(0:jm2,0:km2))
    ! ! call h5_read2dfrom3d(x,im,jm,km,'x','datin/grid.h5',islice=0)
    call H5ReadSubset(y2,im2,jm2,km2,'y','datin/grid.h5',islice=0)
    call H5ReadSubset(z2,im2,jm2,km2,'z','datin/grid.h5',islice=0)
    !
    allocate(ro(0:jm,0:km),u1(0:jm,0:km),u2(0:jm,0:km),u3(0:jm,0:km),t(0:jm,0:km))
    allocate(vary(0:jmm,0:km,1:5))
    allocate(royz2(0:jm2,0:km2),u1yz2(0:jm2,0:km2),u2yz2(0:jm2,0:km2), &
             u3yz2(0:jm2,0:km2), tyz2(0:jm2,0:km2)  )
    !
    do n=nsta,nend
      !
      ! read the first flow field 
      if(n>nmax) then
        write(fname,'(i5.5)')n-nmax
      else
        write(fname,'(i5.5)')n
      endif
      !
      filename=database//'islice/islice'//fname//'.h5'
      !
      call H5ReadArray(time,'time',  trim(filename),explicit=.false.)
      call H5ReadArray(vary(:,:,1),jmm,km,'ro',trim(filename),explicit=.false.)
      call H5ReadArray(vary(:,:,2),jmm,km,'u1',trim(filename),explicit=.false.)
      call H5ReadArray(vary(:,:,3),jmm,km,'u2',trim(filename),explicit=.false.)
      call H5ReadArray(vary(:,:,4),jmm,km,'u3',trim(filename),explicit=.false.)
      call H5ReadArray(vary(:,:,5),jmm,km, 't',trim(filename),explicit=.false.)
      !
      if(n>nmax) then
        time=time+timemax-time0
      endif
      !
      ! store the first flow field to the lower boundary layer
      do j=0,jmm
        !
        if(j>=j2) then
          var1=(y(j,0)-y(j2,0))**2
          var1=exp(-1.d0*var1)
        else
          var1=1.d0
        endif
        !
        ro(j,:)   = (vary(j,:,1)-rom(j))*var1
        u1(j,:)   = (vary(j,:,2)-u1m(j))*var1
        u2(j,:)   = (vary(j,:,3)-u2m(j))*var1
        u3(j,:)   = (vary(j,:,4)       )*var1
         t(j,:)   = (vary(j,:,5)- tm(j))*var1
        !
      enddo
      !
      ! read the second flow field 
      n2=n+3000
      !
      if(n2>nmax) then
        write(fname,'(i5.5)')n2-nmax
      else
        write(fname,'(i5.5)')n2
      endif
      !
      filename=database//'islice/islice'//fname//'.h5'
      !
      call H5ReadArray(vary(:,:,1),jmm,km,'ro',trim(filename),explicit=.false.)
      call H5ReadArray(vary(:,:,2),jmm,km,'u1',trim(filename),explicit=.false.)
      call H5ReadArray(vary(:,:,3),jmm,km,'u2',trim(filename),explicit=.false.)
      call H5ReadArray(vary(:,:,4),jmm,km,'u3',trim(filename),explicit=.false.)
      call H5ReadArray(vary(:,:,5),jmm,km, 't',trim(filename),explicit=.false.)
      !
      ! store the first flow field to the upper boundary layer
      do j=0,jmm
        !
        if(j>=j2) then
          var1=(y(j,0)-y(j2,0))**2
          var1=exp(-1.d0*var1)
        else
          var1=1.d0
        endif
        !
        ro(jm-j,:)   = (vary(j,:,1)-rom(j))*var1
        u1(jm-j,:)   = (vary(j,:,2)-u1m(j))*var1
        u2(jm-j,:)   =-(vary(j,:,3)-u2m(j))*var1
        u3(jm-j,:)   = (vary(j,:,4)       )*var1
         t(jm-j,:)   = (vary(j,:,5)- tm(j))*var1
        !
      enddo
      !
      ! interpolation
      call regularlinearinterp(y,z,ro,u1,u2,u3,t,jm,km,                      &
                               y2,z2,royz2,u1yz2,u2yz2,u3yz2,tyz2,jm2,km2,&
                               homoy=.true.,progress=.false.)
      !
      ! write the data
      write(fname,'(i5.5)')n
      filename='inflow/islice'//fname//'.h5'
      !
      royz2=royz2*ro_ref
      u1yz2=u1yz2*u_ref
      u2yz2=u2yz2*u_ref
      u3yz2=u3yz2*u_ref
       tyz2= tyz2*t_ref
      !
      time=(time-time0)*time_ref
      !
      call H5WriteArray(time,'time',trim(filename),explicit=.false.)
      call H5WriteArray(royz2,jm2,km2,'ro',trim(filename),explicit=.false.)
      call H5WriteArray(u1yz2,jm2,km2,'u1',trim(filename),explicit=.false.)
      call H5WriteArray(u2yz2,jm2,km2,'u2',trim(filename),explicit=.false.)
      call H5WriteArray(u3yz2,jm2,km2,'u3',trim(filename),explicit=.false.)
      call H5WriteArray( tyz2,jm2,km2,'t',trim(filename),explicit=.false.)
      !
      write(*,'(1A1,3(A),I0,A,$)')char(13),'  << ',trim(filename),   &
                                   ' ... ',(n-nsta)*100/(nend-nsta),' %'
      !
    end do
    !
    deallocate(rom,u1m,u2m,tm)
    !
    allocate(rom(0:jm2),u1m(0:jm2),u2m(0:jm2),tm(0:jm2),pyz2(0:jm2,0:km2))
    !
    open(12,file='datin/inlet.prof')
    read(12,*)
    read(12,*)
    read(12,*)
    read(12,*)
    do j=0,jm2
      read(12,*)rom(j),u1m(j),u2m(j),tm(j)
    enddo
    close(12)
    print*,' >> datin/inlet.prof'
    !
    do j=0,jm2
      royz2(j,:)=royz2(j,:)+rom(j)
      u1yz2(j,:)=u1yz2(j,:)+u1m(j)
      u2yz2(j,:)=u2yz2(j,:)+u2m(j)
       tyz2(j,:)= tyz2(j,:)+ tm(j)
      !
    enddo
    !
    pyz2=royz2*tyz2/const2
    !
    call writetecbin('tecro2.plt',y2,'y',z2,'z',royz2,'ro',u1yz2,'u',  &
                          u2yz2,'v',u3yz2,'w',tyz2,'t',pyz2,'p',jm2,km2)
      
  end subroutine isliceinterp_chan
  !+-------------------------------------------------------------------+
  !| The end of the subroutine isliceinterp_chan.                      |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine reads monitor.dat and plot the mon points.        |
  !+-------------------------------------------------------------------+
  subroutine plotmonpoint
    !
    use h5readwrite
    use commvardefine, only: pi,gridfile
    !
    integer :: nmonitor,n,ios
    integer,allocatable :: imon(:),jmon(:)
    real(8),allocatable,dimension(:,:) :: x,y
    real(8),allocatable :: xmon(:),ymon(:),alfa(:)
    integer :: im,jm,km,dim(3)
    !
    dim=h5_getdimensio('x',gridfile)
    print*,dim-1
    !
    im=dim(1)-1
    jm=dim(2)-1
    km=dim(3)-1
    !
    allocate(x(0:im,0:jm),y(0:im,0:jm))
    call h5_read2dfrom3d(x,im,jm,km,'x',gridfile,kslice=0)
    call h5_read2dfrom3d(y,im,jm,km,'y',gridfile,kslice=0)
    !
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
    allocate(imon(nmonitor),jmon(nmonitor))
    open(12,file='datin/monitor.dat')
    read(12,'()')
    do n=1,nmonitor
      read(12,*)imon(n),jmon(n)
    enddo
    close(12)
    print*,' >> datin/monitor.dat'
    !
    ! allocate(alfa(nmonitor))
    ! do n=1,nmonitor
    !   alfa(n)=asin(y(imon(n),jmon(n))/sqrt((x(imon(n),jmon(n))-1.d0)**2+y(imon(n),jmon(n))**2))
    !   alfa(n)=alfa(n)/pi*180.d0
    ! enddo
    ! open(18,file='alfa_monitor.dat')
    ! do n=1,nmonitor
    !   print*,n,alfa(n)
    !   write(18,*)n,alfa(n)
    ! enddo
    ! close(18)
    ! print*,' << alfa_monitor.dat'
    ! !
    open(18,file='xyzmonitor.dat')
    write(18,'(A)')'TITLE     = "mointor points"'
    write(18,'(A)')'VARIABLES = "x" "y"'
    write(18,'(A)')'ZONE T="ZONE 001"'
    write(18,'(A)')'STRANDID=0, SOLUTIONTIME=0'
    write(18,'(A,I0,A)')'I=',nmonitor,', ZONETYPE=Ordered'
    write(18,'(A)')'DATAPACKING=POINT'
    write(18,'(A)')'DT=(SINGLE SINGLE )'
    do n=1,nmonitor
      write(18,*)x(imon(n),jmon(n)),y(imon(n),jmon(n))
    enddo
    close(18)
    print*,' << xyzmonitor.dat'
    !
  end subroutine plotmonpoint
  !+-------------------------------------------------------------------+
  !| The end of the subroutine isliceinterp.                           |
  !+-------------------------------------------------------------------+
  !
  subroutine pwallmon
    !
    use h5readwrite
    use commvardefine, only: pi,pinf
    !
    integer :: nmonitor,n,n1,i,icutoff
    integer,allocatable :: imon(:),jmon(:),nmon(:)
    real(8),allocatable,dimension(:,:) :: x,y
    integer :: im,jm,km,dim(3),nmax,nstep,nvar,istart
    real(8) :: rvar,lref,tref,var1,dfeq,feq,utaw,miu,row,tawx,ro,miuw,tintv,psdint,uref
    real(8),allocatable :: pwall(:,:),time(:),pmean(:),pf(:,:),prms(:)
    real(8),allocatable :: psd(:),wpsd(:)
    character(len=4) :: mfname
    character(len=64) :: char
    !
    dim=h5_getdimensio('x','datin/grid.2d')
    print*,dim-1
    !
    im=dim(1)-1
    jm=dim(2)-1
    km=dim(3)-1
    !
    ! allocate(x(0:im,0:jm),y(0:im,0:jm))
    ! call h5_read2dfrom3d(x,im,jm,km,'x','datin/grid.h5',kslice=0)
    ! call h5_read2dfrom3d(y,im,jm,km,'y','datin/grid.h5',kslice=0)
    !
    nmonitor=6
    nmax=679100
    !
    uref=1831.678407d0
    tref=1.d0/uref
    !
    allocate(time(1:nmax),pwall(1:nmonitor,1:nmax),imon(1:nmonitor),nmon(1:nmonitor))
    allocate(pmean(1:nmonitor),prms(1:nmonitor),pf(1:nmonitor,1:nmax))
    !
    nmon(1)=1
    nmon(2)=4
    nmon(3)=8
    nmon(4)=12
    nmon(5)=17
    nmon(6)=22
    ! nmon(6)=19
    !
    do n=1,nmonitor
      !
      write(mfname,'(I4.4)')nmon(n)
      !
      open(16,file='monitor.'//mfname//'.dat')
      read(16,*)
      i=1
      do while(i<=nmax)
        read(16,*)n1,time(i),var1,var1,pwall(n,i)
        i=i+1
      enddo
      ! read(16,'(A)')char
      ! read(16,'(A)')char
      ! !
      ! nstep=0
      ! do while(nstep<639159)
      !   !
      !   read(16,*)nstep,time(nstep),rvar,rvar,rvar,pwall(n,nstep)
      !   !
      ! enddo
      ! close(16)
      print*,' >> monitor.',mfname,'.dat'
      !
    enddo
    !
    ! time=time*tref
    ! pwall=pwall*pinf
    !
    ! open(18,file='Results/pwall_monitor.dat')
    ! do i=1,nmax
    !   write(18,"(7(1X,E15.7E3))")time(i),(pwall(n,i),n=1,nmonitor)
    ! enddo
    ! close(18)
    ! print*,' << Results/pwall_monitor.dat'
    !
    ! do i=1,nmax
    !   !
    !   if(time(i)>=600.d0) then
    !     istart=i
    !     exit
    !   endif
    !   !
    ! enddo
    istart=500000
    print*,' ** istart=',istart
    !
    do n=1,nmonitor
      !
      write(mfname,'(I4.4)')nmon(n)
      !
      pmean(n)=0.d0
      do i=istart,nmax
        pmean(n)=pmean(n)+pwall(n,i)
      enddo
      pmean(n)=pmean(n)/dble(nmax-istart+1)
      !
      print*,' ** pmean=',pmean(n)
      !
      prms(n)=0.d0
      do i=istart,nmax
        !
        pf(n,i)=pwall(n,i)-pmean(n)
        prms(n)=prms(n)+pf(n,i)**2
      enddo
      prms(n)=sqrt(prms(n)/dble(nmax-istart+1))
      !
      print*,' ** prms=',prms(n)
      !
      pf(n,:)=pf(n,:)/prms(n)
      !
      ! psd=freqspetra_overlap(pf(n,istart:nmax),8)
      ! !
      ! psdint=0.d0
      ! do i=1,size(psd)
      !   psdint=psdint+psd(i)
      ! enddo
      !
      psd=freqspetra(pf(n,istart:nmax))
      !
      print*,' ** size of pds',size(psd)
      !
      tintv=(time(2)-time(1))*dble(size(psd)*2+1)
      print*,'the base T periodcity:',tintv,'~ frequency:',1.d0/tintv
      !
      print*,' ** psd calculated .'
      !
      icutoff=1000
      !
      var1=0.d0
      dfeq=2.d0*pi/tintv
      do i=2,icutoff
        var1=var1+0.5d0*(psd(i)+psd(i+1))*dfeq
      enddo
      allocate(wpsd(1:icutoff))
      do i=1,icutoff
        ! feq=(i-1)/(time(nmax)-time(istart))
        feq=2.d0*pi*(i-1)/tintv
        wpsd(i)=feq*psd(i)/var1
      enddo
      !
      ! do i=7300,10000
      !   ! var1=tanh(dble(i-7300)*0.001d0)*2.d0
      !   var1=dble(i-7300)/5000.d0+1.d0
      !   psd(i)=psd(i)*1.d0/var1**1.01
      !   print*,i,1.d0/var1**1.01
      ! enddo
      ! ** tawx=   1.0696974238552227E-003
      ! utaw=   5.0404101841396015E-002
      ! miu=   1.d0/66420.00 !
      ! miuw=3.2955970957377128E-005
      ! ro=  1.d0 !0.42104565416025225
      ! row=0.42104565416025225
      ! tawx=   1.0696974238552227E-003
      lref=0.01048d0
      uref=1.d0
      utaw=5.564215d-2
      tawx=2.2056597514574835d-3
      !
      open(18,file='spec.pressure.'//mfname//'.dat')
      write(18,"(4(1X,A15))")'f','st','psd','wpsd'
      do i=1,icutoff
        write(18,"(4(1X,E15.7E3))")dble(i-1)/tintv,dble(i-1)/tintv*lref/uref,psd(i),wpsd(i)
        ! write(18,"(2(1X,E15.7E3))")2.d0*pi*(i-1)/(time(nmax)-time(istart))*lref,psd(i)/(0.5*lref)
        ! write(18,"(3(1X,E15.7E3))")2.d0*pi*(i-1)/tintv*lref/utaw, psd(i)*utaw/tawx/tawx*lref, wpsd(i)
        ! write(18,"(3(1X,E15.7E3))")2.d0*pi*(i-1)/tintv*lref,psd(i)/0.25d0*lref,wpsd(i)
      end do
      close(18)
      print*,' << spec.pressure.',mfname,'.dat'
      !
      deallocate(psd,wpsd)
      !
    enddo
    !
      ! !
      ! open(18,file='spec.pressure.'//mfname//'.dat')
      ! write(18,"(2(1X,A15))")'frequency','psd'
      ! do i=1,size(spec)
      !   write(18,"(2(1X,E15.7E3))")2.d0*pi*i/(timemon(nline)-timemon(nsline))*lref, &
      !                             spec(i)/(0.5d0*lref)
      ! end do
      ! close(18)
      ! print*,' << spec.pressure.',mfname,'.dat'
  end subroutine pwallmon
  !
  subroutine psdtest
    !
    use commvardefine, only: pi,pinf
    !
    integer,parameter :: dim=1000
    real(8) :: pf(1:dim),time(1:dim)
    real(8),allocatable :: psd(:)
    integer :: i
    !
    do i=1,dim
      time(i)=2.d0*pi/dim*i
      pf(i)=sin(time(i))+0.5d0*sin(2.d0*time(i)+0.5d0*pi)
    enddo
    psd=freqspetra(pf)
    !
    open(18,file='test.dat')
    do i=1,dim
      write(18,"(2(1X,E15.7E3))")time(i),pf(i)
    end do
    close(18)
    print*,' << test.dat'
    !
    open(18,file='spec_test.dat')
    do i=1,size(psd)
      write(18,"(1X,I4,1X,E15.7E3)")i,psd(i)
    end do
    close(18)
    print*,' << spec_test.dat'
    !
  end subroutine psdtest
  !
  function freqspetra(p) result(spc)
    !
    use singleton
    !
    real(8),intent(in) :: p(:)
    real(8),allocatable :: spc(:)
    !
    integer :: dim,i
    complex(8),allocatable,dimension(:) :: pfc
    !
    dim=size(p)
    !
    allocate(spc(1:dim/2))
    allocate(pfc(dim))
    !
    pfc=p*(1.d0,0.d0)
    !
    pfc=FFT(pfc,inv=.false.)/sqrt(dble(dim))*2.d0
    !
    ! pfc=FFT(pcc,inv=.true.)
    !
    ! print*,pcc(1)
    ! print*,'----------------------------------------'
    !
    do i=1,dim/2
      !
      spc(i)=real(pfc(i))**2+aimag(pfc(i))**2
      !
      ! print*,i,ZABS(pcc(i)),ZABS(pcc(dim-i+2))
      ! ! print*,i,aimag(pcc(i)),aimag(pcc(dim-i+2))
      ! print*,'----------------------------------------'
    enddo
    !
    return
    !
  end function freqspetra
  !
  function freqspetra_overlap(p,nintval) result(spc)
    !
    use singleton
    !
    real(8),intent(in) :: p(:)
    real(8),allocatable :: spc(:)
    integer,intent(in) :: nintval
    !
    integer :: dim_t,dim,i,n,ns
    complex(8),allocatable,dimension(:) :: pfc
    !
    dim_t=size(p)
    !
    print*,' ** total size of the p array:',dim_t
    !
    dim=(dim_t/(nintval+1))*2
    !
    print*,' ** size of each interval:',dim
    !
    allocate(spc(1:dim/2))
    allocate(pfc(dim))
    !
    spc=0.d0
    !
    do n=1,nintval
      !
      ns=dim/2*(n-1)
      !
      print*,n,':',ns+1,'~',ns+dim
      !
      pfc(1:dim)=p(ns+1:ns+dim)*(1.d0,0.d0)
      !
      pfc=FFT(pfc,inv=.false.)/sqrt(dble(dim))*2.d0
      !
      do i=1,dim/2
        spc(i)=spc(i)+real(pfc(i))**2+aimag(pfc(i))**2
      enddo
      !
    enddo
    !
    spc=spc/dble(nintval)
    !
    ! pfc=p*(1.d0,0.d0)
    ! !
    ! pfc=FFT(pfc,inv=.false.)/sqrt(dble(dim))*2.d0
    ! !
    ! ! pfc=FFT(pcc,inv=.true.)
    ! !
    ! ! print*,pcc(1)
    ! ! print*,'----------------------------------------'
    ! !
    ! do i=1,dim/2
    !   !
    !   spc(i)=real(pfc(i))**2+aimag(pfc(i))**2
    !   !
    !   ! print*,i,ZABS(pcc(i)),ZABS(pcc(dim-i+2))
    !   ! ! print*,i,aimag(pcc(i)),aimag(pcc(dim-i+2))
    !   ! print*,'----------------------------------------'
    ! enddo
    !
    return
    !
  end function freqspetra_overlap
  !
  subroutine spec_homo(nsta,nend)
    !
    use commvardefine,only: im,jm,km,u1,u2,u3,x,y,z,Reynolds,outfolder
    use singleton
    use h5readwrite
    !
    integer,intent(in) :: nsta,nend
    !
    complex(8), allocatable, dimension(:,:,:) :: u1ct,u2ct,u3ct,u1c,u2c,u3c
    real(8),allocatable :: spce(:),wavn(:)
    !
    character(len=4) :: fname
    character(len=64) :: filename
    real(8) :: u1m,u2m,u3m,time,engs
    integer :: i,j,k,n,kw,imm,jmm,kmm
    !
    imm=im/2
    jmm=jm/2
    kmm=km/2
    !
    allocate(u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),                    &
             u3(0:im,0:jm,0:km))
    allocate(u1c(1:im,1:jm,1:km),u2c(1:im,1:jm,1:km),                  &  
             u3c(1:im,1:jm,1:km) )
    allocate(u1ct(1:im,1:jm,1:km),u2ct(1:im,1:jm,1:km),              &  
             u3ct(1:im,1:jm,1:km) )
    allocate(spce(0:kmm**3),wavn(0:kmm**3))
    !
    do n=nsta,nend
      !
      write(fname,'(i4.4)') n
      filename=outfolder//'flowfield'//fname//'.h5'
      !
      call H5ReadArray(time,'time',trim(filename))
      print*,' ** time=',time
      call H5ReadArray(u1,im,jm,km,'u1',trim(filename))
      call H5ReadArray(u2,im,jm,km,'u2',trim(filename))
      call H5ReadArray(u3,im,jm,km,'u3',trim(filename))
      !
      u1m=0.d0
      u2m=0.d0
      u3m=0.d0
      do k=1,km
      do j=1,jm
      do i=1,im
        u1m=u1m+u1(i,j,k)
        u2m=u2m+u2(i,j,k)
        u3m=u3m+u3(i,j,k)
      enddo
      enddo
      enddo
      u1m=u1m/dble(im*jm*km)
      u2m=u2m/dble(im*jm*km)
      u3m=u3m/dble(im*jm*km)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        u1c(i,j,k)=cmplx(u1(i,j,k)-u1m,0.d0,kind=8)
        u2c(i,j,k)=cmplx(u2(i,j,k)-u2m,0.d0,kind=8)
        u3c(i,j,k)=cmplx(u3(i,j,k)-u3m,0.d0,kind=8)
      enddo
      enddo
      enddo
      !
      spce=0.d0
      !
      do k=1,km
      do j=1,jm
        u1ct(:,j,k)=FFT(u1c(:,j,k),inv=.false.)/sqrt(dble(im))
        u2ct(:,j,k)=FFT(u2c(:,j,k),inv=.false.)/sqrt(dble(im))
        u3ct(:,j,k)=FFT(u3c(:,j,k),inv=.false.)/sqrt(dble(im))
      enddo
      enddo
      !
      do i=0,imm
        do k=1,km
        do j=1,jm
          spce(i)=spce(i)+0.5d0*( (real(u1ct(i+1,j,k),8)**2+aimag(u1ct(i+1,j,k))**2) + &
                                  (real(u2ct(i+1,j,k),8)**2+aimag(u2ct(i+1,j,k))**2) + &
                                  (real(u3ct(i+1,j,k),8)**2+aimag(u3ct(i+1,j,k))**2) )
        enddo
        enddo
      enddo
      !
      do k=1,km
      do i=1,im
        u1ct(i,:,k)=FFT(u1c(i,:,k),inv=.false.)/sqrt(dble(jm))
        u2ct(i,:,k)=FFT(u2c(i,:,k),inv=.false.)/sqrt(dble(jm))
        u3ct(i,:,k)=FFT(u3c(i,:,k),inv=.false.)/sqrt(dble(jm))
      enddo
      enddo
      !
      do j=0,jmm
        do k=1,km
        do i=1,im
          spce(j)=spce(j)+0.5d0*( (real(u1ct(i,j+1,k),8)**2+aimag(u1ct(i,j+1,k))**2) + &
                                  (real(u2ct(i,j+1,k),8)**2+aimag(u2ct(i,j+1,k))**2) + &
                                  (real(u3ct(i,j+1,k),8)**2+aimag(u3ct(i,j+1,k))**2) )
        enddo
        enddo
      enddo
      !
      do j=1,jm
      do i=1,im
        u1ct(i,j,:)=FFT(u1c(i,j,:),inv=.false.)/sqrt(dble(km))
        u2ct(i,j,:)=FFT(u2c(i,j,:),inv=.false.)/sqrt(dble(km))
        u3ct(i,j,:)=FFT(u3c(i,j,:),inv=.false.)/sqrt(dble(km))
      enddo
      enddo
      !
      do k=0,kmm
        do j=1,jm
        do i=1,im
          spce(k)=spce(k)+0.5d0*( (real(u1ct(i,j,k+1),8)**2+aimag(u1ct(i,j,k+1))**2) + &
                                  (real(u2ct(i,j,k+1),8)**2+aimag(u2ct(i,j,k+1))**2) + &
                                  (real(u3ct(i,j,k+1),8)**2+aimag(u3ct(i,j,k+1))**2) )
        enddo
        enddo
      enddo
      !
      spce=spce/dble(im*jm+im*km+jm*km)
      !
      engs=0.d0
      do k=0,kmm
        engs=engs+spce(k)
      enddo
      engs=engs/dble(kmm+1)

      ! kw=0
      ! wavn=0.d0
      ! spce=0.d0
      ! do k=0,kmm-1
      ! do j=0,jmm-1
      ! do i=0,imm-1
      !   kw=kw+1
      !   wavn(kw)=sqrt(dble(i)**2+dble(j)**2+dble(k)**2)
      !   spce(kw)=0.5d0*( (real(u1ct(i+1,j+1,k+1),8)**2+aimag(u1ct(i+1,j+1,k+1))**2) + &
      !                    (real(u2ct(i+1,j+1,k+1),8)**2+aimag(u2ct(i+1,j+1,k+1))**2) + &
      !                    (real(u3ct(i+1,j+1,k+1),8)**2+aimag(u3ct(i+1,j+1,k+1))**2) )

      ! enddo
      ! enddo
      ! enddo
      !
      open(18,file='spectra'//fname//'.dat')
      do kw=0,km 
        write(18,*)dble(kw),spce(kw),spce(kw)/engs
      enddo
      close(18)
      print*,' << spectra',fname,'.dat'
      !
    enddo
    !
  end subroutine spec_homo
  !
  subroutine spec_homo_onefile(fileinput,fileoutput)
    !
    use commvardefine,only: im,jm,km,u1,u2,u3,x,y,z,Reynolds,outfolder
    use singleton
    use h5readwrite
    !
    character(len=*),intent(in) :: fileinput,fileoutput
    !
    complex(8), allocatable, dimension(:,:,:) :: u1ct,u2ct,u3ct,u1c,u2c,u3c
    real(8),allocatable :: spce(:),wavn(:)
    !
    real(8) :: u1m,u2m,u3m,time,engs
    integer :: i,j,k,kw,imm,jmm,kmm
    !
    imm=im/2
    jmm=jm/2
    kmm=km/2
    !
    allocate(u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),                    &
             u3(0:im,0:jm,0:km))
    allocate(u1c(1:im,1:jm,1:km),u2c(1:im,1:jm,1:km),                  &  
             u3c(1:im,1:jm,1:km) )
    allocate(u1ct(1:im,1:jm,1:km),u2ct(1:im,1:jm,1:km),              &  
             u3ct(1:im,1:jm,1:km) )
    allocate(spce(0:kmm**3),wavn(0:kmm**3))
    !
    call H5ReadArray(time,'time',trim(fileinput))
    print*,' ** time=',time
    call H5ReadArray(u1,im,jm,km,'u1',trim(fileinput))
    call H5ReadArray(u2,im,jm,km,'u2',trim(fileinput))
    call H5ReadArray(u3,im,jm,km,'u3',trim(fileinput))
    !
    u1m=0.d0
    u2m=0.d0
    u3m=0.d0
    do k=1,km
    do j=1,jm
    do i=1,im
      u1m=u1m+u1(i,j,k)
      u2m=u2m+u2(i,j,k)
      u3m=u3m+u3(i,j,k)
    enddo
    enddo
    enddo
    u1m=u1m/dble(im*jm*km)
    u2m=u2m/dble(im*jm*km)
    u3m=u3m/dble(im*jm*km)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      u1c(i,j,k)=cmplx(u1(i,j,k)-u1m,0.d0,kind=8)
      u2c(i,j,k)=cmplx(u2(i,j,k)-u2m,0.d0,kind=8)
      u3c(i,j,k)=cmplx(u3(i,j,k)-u3m,0.d0,kind=8)
    enddo
    enddo
    enddo
    !
    spce=0.d0
    !
    do k=1,km
    do j=1,jm
      u1ct(:,j,k)=FFT(u1c(:,j,k),inv=.false.)/sqrt(dble(im))
      u2ct(:,j,k)=FFT(u2c(:,j,k),inv=.false.)/sqrt(dble(im))
      u3ct(:,j,k)=FFT(u3c(:,j,k),inv=.false.)/sqrt(dble(im))
    enddo
    enddo
    !
    do i=0,imm
      do k=1,km
      do j=1,jm
        spce(i)=spce(i)+0.5d0*( (real(u1ct(i+1,j,k),8)**2+aimag(u1ct(i+1,j,k))**2) + &
                                (real(u2ct(i+1,j,k),8)**2+aimag(u2ct(i+1,j,k))**2) + &
                                (real(u3ct(i+1,j,k),8)**2+aimag(u3ct(i+1,j,k))**2) )
      enddo
      enddo
    enddo
    !
    do k=1,km
    do i=1,im
      u1ct(i,:,k)=FFT(u1c(i,:,k),inv=.false.)/sqrt(dble(jm))
      u2ct(i,:,k)=FFT(u2c(i,:,k),inv=.false.)/sqrt(dble(jm))
      u3ct(i,:,k)=FFT(u3c(i,:,k),inv=.false.)/sqrt(dble(jm))
    enddo
    enddo
    !
    do j=0,jmm
      do k=1,km
      do i=1,im
        spce(j)=spce(j)+0.5d0*( (real(u1ct(i,j+1,k),8)**2+aimag(u1ct(i,j+1,k))**2) + &
                                (real(u2ct(i,j+1,k),8)**2+aimag(u2ct(i,j+1,k))**2) + &
                                (real(u3ct(i,j+1,k),8)**2+aimag(u3ct(i,j+1,k))**2) )
      enddo
      enddo
    enddo
    !
    do j=1,jm
    do i=1,im
      u1ct(i,j,:)=FFT(u1c(i,j,:),inv=.false.)/sqrt(dble(km))
      u2ct(i,j,:)=FFT(u2c(i,j,:),inv=.false.)/sqrt(dble(km))
      u3ct(i,j,:)=FFT(u3c(i,j,:),inv=.false.)/sqrt(dble(km))
    enddo
    enddo
    !
    do k=0,kmm
      do j=1,jm
      do i=1,im
        spce(k)=spce(k)+0.5d0*( (real(u1ct(i,j,k+1),8)**2+aimag(u1ct(i,j,k+1))**2) + &
                                (real(u2ct(i,j,k+1),8)**2+aimag(u2ct(i,j,k+1))**2) + &
                                (real(u3ct(i,j,k+1),8)**2+aimag(u3ct(i,j,k+1))**2) )
      enddo
      enddo
    enddo
    !
    spce=spce/dble(im*jm+im*km+jm*km)
    !
    engs=0.d0
    do k=0,kmm
      engs=engs+spce(k)
    enddo
    engs=engs/dble(kmm+1)

    ! kw=0
    ! wavn=0.d0
    ! spce=0.d0
    ! do k=0,kmm-1
    ! do j=0,jmm-1
    ! do i=0,imm-1
    !   kw=kw+1
    !   wavn(kw)=sqrt(dble(i)**2+dble(j)**2+dble(k)**2)
    !   spce(kw)=0.5d0*( (real(u1ct(i+1,j+1,k+1),8)**2+aimag(u1ct(i+1,j+1,k+1))**2) + &
    !                    (real(u2ct(i+1,j+1,k+1),8)**2+aimag(u2ct(i+1,j+1,k+1))**2) + &
    !                    (real(u3ct(i+1,j+1,k+1),8)**2+aimag(u3ct(i+1,j+1,k+1))**2) )

    ! enddo
    ! enddo
    ! enddo
    !
    open(18,file=trim(fileoutput))
    do kw=0,km 
      write(18,*)dble(kw),spce(kw),spce(kw)/engs
    enddo
    close(18)
    print*,' << ',trim(fileoutput)
    !
    !
  end subroutine spec_homo_onefile
  !
  subroutine pwcorrelation(nsta,nend)
    !
    use commvardefine, only: im,jm,km,outfolder,gridfile
    use gradsolver
    use h5readwrite
    !
    integer :: nsta,nend,n,k,i,j
    !
    real(8),allocatable,dimension(:,:,:) :: u3,p
    real(8),allocatable,dimension(:,:) :: pzm,c_dpdz_w,ww,dpz2
    real(8),allocatable,dimension(:) :: z1d,dpdz
    real(8) :: time
    character(len=4) sname
    character(len=23) fname
    !
    allocate(z1d(0:km),dpdz(0:km))
    call H5ReadSubset(z1d,im,jm,km,'z',trim(gridfile),islice=0,jslice=0)
    !
    allocate(u3(0:im,0:jm,0:km),p(0:im,0:jm,0:km)   )
    allocate(c_dpdz_w(0:im,0:jm),pzm(0:im,0:jm),ww(0:im,0:jm),dpz2(0:im,0:jm))
    !
    call H5ReadArray(pzm, im,jm, 'p','Results/mean.fav.zm.h5')
    !
    c_dpdz_w=0.d0
    ww=0.d0
    dpz2=0.d0
    !
    do n=nsta,nend
      !
      write(sname,'(i4.4)')n
      fname=outfolder//'flowfield'//sname//'.h5'
      call H5ReadArray(time,'time',fname)
      write(*,'(A,F10.5)')'  ** << flow field  at t=',time
      !
      call H5ReadArray(u3,im,jm,km,'u3',fname)
      call H5ReadArray( p,im,jm,km, 'p',fname)
      !
      do k=0,km
        p(:,:,k)=p(:,:,k)-pzm
      enddo
      !
      do j=0,jm
      do i=0,im
        !
        dpdz=dfdk(p(i,j,:))/z1d(1)
        !
        do k=1,km
          c_dpdz_w(i,j)=c_dpdz_w(i,j)+dpdz(k)*u3(i,j,k)
          ww(i,j)=ww(i,j)+u3(i,j,k)*u3(i,j,k)
          dpz2(i,j)=dpz2(i,j)+dpdz(k)*dpdz(k)
        enddo
        !
      enddo
      enddo
      !
    enddo
    !
    c_dpdz_w=c_dpdz_w/dble(km*(nend-nsta+1))
    ww=ww/dble(km*(nend-nsta+1))
    dpz2=dpz2/dble(km*(nend-nsta+1))
    !
    call H5WriteArray(c_dpdz_w,im,jm,'c_dpdz_w','Results/c_dpdz_w.h5')
    call H5WriteArray(      ww,im,jm,      'ww','Results/c_dpdz_w.h5')
    call H5WriteArray(    dpz2,im,jm,    'dpz2','Results/c_dpdz_w.h5')
      
    !
  end subroutine pwcorrelation
  !
  subroutine flowrecons
    !
    use commvardefine,only: im,jm,km,Reynolds
    use h5readwrite
    use WriteTec
    !
    real(8),allocatable,dimension(:,:,:) :: ro,u1,u2,u3,p,t,x,y,z
    integer :: i,j,k
    !
    allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
    allocate(ro(0:im,0:jm,0:km),u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km), &
             u3(0:im,0:jm,0:km), p(0:im,0:jm,0:km), t(0:im,0:jm,0:km))
    !
    call H5ReadArray(x,im,jm,km,'x','datin/grid.h5')
    call H5ReadArray(y,im,jm,km,'y','datin/grid.h5')
    call H5ReadArray(z,im,jm,km,'z','datin/grid.h5')
    !
    call H5ReadArray(ro,im,jm,km,'ro','bakup/flowfield.h5')
    call H5ReadArray(u1,im,jm,km,'u1','bakup/flowfield.h5')
    call H5ReadArray(u2,im,jm,km,'u2','bakup/flowfield.h5')
    call H5ReadArray(u3,im,jm,km,'u3','bakup/flowfield.h5')
    call H5ReadArray( p,im,jm,km, 'p','bakup/flowfield.h5')
    call H5ReadArray( t,im,jm,km, 't','bakup/flowfield.h5')
    !
    do k=0,km
    do j=0,jm
      !
      do i=0,im
        !
        if(x(i,j,k)>=60.d0) then
          ro(i,j,k)=ro(i-1,j,k)
          u1(i,j,k)=u1(i-1,j,k)
          u2(i,j,k)=u2(i-1,j,k)
          u3(i,j,k)=u3(i-1,j,k)
           p(i,j,k)= p(i-1,j,k)
           t(i,j,k)= t(i-1,j,k)
        endif
        !
      enddo
      !
    enddo
    enddo
    !
    call H5WriteArray(ro,im,jm,km,'ro','flowini3d.h5')
    call H5WriteArray(u1,im,jm,km,'u1','flowini3d.h5')
    call H5WriteArray(u2,im,jm,km,'u2','flowini3d.h5')
    call H5WriteArray(u3,im,jm,km,'u3','flowini3d.h5')
    call H5WriteArray( p,im,jm,km, 'p','flowini3d.h5')
    call H5WriteArray( t,im,jm,km, 't','flowini3d.h5')
    !
    call writetecbin('tecinixy.plt',x(:,:,km/2),'x',y(:,:,km/2),'y', &
                      ro(:,:,km/2),'ro',u1(:,:,km/2),'u',            &
                      u2(:,:,km/2),'v',u3(:,:,km/2),'w',             &
                      p(:,:,km/2),'p',t(:,:,km/2),'t',im,jm)
  end subroutine flowrecons
  !
end module dataprocess
!+---------------------------------------------------------------------+
!| The end of the module statistic.                                    |
!+---------------------------------------------------------------------+