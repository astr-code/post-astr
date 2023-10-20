!+---------------------------------------------------------------------+
!| This module contains subroutines to calculate statistics.           |
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2020-04-16                                     |
!+---------------------------------------------------------------------+
module statistic
  !
  use commvardefine, only : outfolder
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| The subroutine is used to read meanflow and do spatial averaging. |
  !+-------------------------------------------------------------------+
  !| Writen by Jian Fang, 2020-04-18                                   |
  !+-------------------------------------------------------------------+
  subroutine meanflow_fav(nsta,nend)
    !
    use commvardefine, only: im,jm,km,ro_m,u1_m,u2_m,u3_m,p_m,t_m,     &
                             ro_zm,u1_zm,u2_zm,u3_zm,p_zm,t_zm,        &
                             ro_xzm,u1_xzm,u2_xzm,u3_xzm,p_xzm,t_xzm,  &
                             lihomo,ljhomo,lkhomo,x,y
    !
    use h5readwrite
    use writetec
    use basicfunction, only: integration
    !
    integer :: nsta,nend,n
    !
    integer :: nsamples
    integer :: i,j,k
    real(8) :: time
    logical :: filexists=.false.
    character(len=4) sname
    character(len=23) fname
    !
    real(8),allocatable,dimension(:,:,:) :: ro,u1,u2,u3,p,t
    !
    allocate(ro(0:im,0:jm,0:km),u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km), &
             u3(0:im,0:jm,0:km),p(0:im,0:jm,0:km), t(0:im,0:jm,0:km)   )
    allocate(ro_m(0:im,0:jm,0:km),u1_m(0:im,0:jm,0:km),                &
             u2_m(0:im,0:jm,0:km),u3_m(0:im,0:jm,0:km),                &
              p_m(0:im,0:jm,0:km), t_m(0:im,0:jm,0:km)                 )
    !
    if(nsta>=0 .and. nsta<=99999) then
      !
      ro_m=0.d0
      u1_m=0.d0
      u2_m=0.d0
      u3_m=0.d0
       p_m=0.d0
       t_m=0.d0
      do n=nsta,nend
        write(sname,'(i4.4)')n
        fname=outfolder//'flowfield'//sname//'.h5'
        call H5ReadArray(time,'time',fname)
        write(*,'(A,F10.5)')'  ** << flow field  at t=',time
        !
        call H5ReadArray(ro,im,jm,km,'ro',fname)
        call H5ReadArray(u1,im,jm,km,'u1',fname)
        call H5ReadArray(u2,im,jm,km,'u2',fname)
        call H5ReadArray(u3,im,jm,km,'u3',fname)
        call H5ReadArray( p,im,jm,km, 'p',fname)
        call H5ReadArray( t,im,jm,km, 't',fname)
        !
        ro_m=ro_m+ro      
        u1_m=u1_m+u1*ro
        u2_m=u2_m+u2*ro
        u3_m=u3_m+u3*ro
         p_m= p_m+ p
         t_m= t_m+ t*ro
        !
      enddo

      !
      nsamples=nend-nsta+1
      !
    else
      !
      call H5ReadArray(nsamples,'nsamples',outfolder//'meanflow.h5')
      !
      write(*,'(A,I5)')'  ** Total number of temporal samples:',nsamples
      !
      call H5ReadArray(ro_m,im,jm,km,'rom',outfolder//'meanflow.h5')
      call H5ReadArray(u1_m,im,jm,km,'u1m',outfolder//'meanflow.h5')
      call H5ReadArray(u2_m,im,jm,km,'u2m',outfolder//'meanflow.h5')
      call H5ReadArray(u3_m,im,jm,km,'u3m',outfolder//'meanflow.h5')
      call H5ReadArray( p_m,im,jm,km, 'pm',outfolder//'meanflow.h5')
      call H5ReadArray( t_m,im,jm,km, 'tm',outfolder//'meanflow.h5')
    endif
    !
    if(lkhomo .and. lihomo) then
      allocate( ro_xzm(0:jm),u1_xzm(0:jm),u2_xzm(0:jm),u3_xzm(0:jm),   &
                 p_xzm(0:jm), t_xzm(0:jm) )
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j)
      !
      !$OMP DO
      do j=0,jm
        ro_xzm(j)=integration(ro_m(:,j,:))
        u1_xzm(j)=integration(u1_m(:,j,:))
        u2_xzm(j)=integration(u2_m(:,j,:))
        u3_xzm(j)=integration(u3_m(:,j,:))
         p_xzm(j)=integration( p_m(:,j,:))
         t_xzm(j)=integration( t_m(:,j,:))
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      print*,' ** averaged in i and k direction'
      !
      u1_xzm=u1_xzm/ro_xzm
      u2_xzm=u2_xzm/ro_xzm
      u3_xzm=u3_xzm/ro_xzm
      !
       t_xzm= t_xzm/ro_xzm
      !
      ro_xzm=ro_xzm/real(nsamples*im*km,8)
       p_xzm= p_xzm/real(nsamples*im*km,8)
      !
      inquire(file='Results/mean.fav.xzm.h5',exist=filexists)
      if(filexists) then
        call system('mv -v Results/mean.fav.xzm.h5 Results/mean.fav.zm.bak')
      endif
      !
      call H5WriteArray(ro_xzm,jm,'ro','Results/mean.fav.xzm.h5')
      call H5WriteArray(u1_xzm,jm,'u1','Results/mean.fav.xzm.h5')
      call H5WriteArray(u2_xzm,jm,'u2','Results/mean.fav.xzm.h5')
      call H5WriteArray(u3_xzm,jm,'u3','Results/mean.fav.xzm.h5')
      call H5WriteArray( p_xzm,jm, 'p','Results/mean.fav.xzm.h5')
      call H5WriteArray( t_xzm,jm, 't','Results/mean.fav.xzm.h5')
      !
      open(18,file='Results/profile.dat')
      write(18,"(6(1X,A15))")'y','ro','u','v','T','P'
      write(18,"(6(1X,E15.7E3))")(y(0,j,0),ro_xzm(j),u1_xzm(j),        &
                                  u2_xzm(j),t_xzm(j),p_xzm(j),j=0,jm)
      close(18)
      print*,' << profile.dat ... done. '
      !
    endif
    !
    if(lkhomo) then
      ! average in z direction
      !
      allocate(ro_zm(0:im,0:jm),u1_zm(0:im,0:jm),u2_zm(0:im,0:jm),     &
               u3_zm(0:im,0:jm), p_zm(0:im,0:jm), t_zm(0:im,0:jm))
      
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
      !
      !$OMP DO
      do j=0,jm
      do i=0,im
        ro_zm(i,j)=integration(ro_m(i,j,:))
        u1_zm(i,j)=integration(u1_m(i,j,:))
        u2_zm(i,j)=integration(u2_m(i,j,:))
        u3_zm(i,j)=integration(u3_m(i,j,:))
         p_zm(i,j)=integration( p_m(i,j,:))
         t_zm(i,j)=integration( t_m(i,j,:))
      enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      u1_zm=u1_zm/ro_zm
      u2_zm=u2_zm/ro_zm
      u3_zm=u3_zm/ro_zm
      !
       t_zm= t_zm/ro_zm
      !
      ro_zm=ro_zm/dble(nsamples)
       p_zm= p_zm/dble(nsamples)
      !
      if(km>0) then
        ro_zm=ro_zm/dble(km)
        p_zm= p_zm/dble(km)
      endif
      !
      inquire(file='Results/mean.fav.zm.h5',exist=filexists)
      if(filexists) then
        call system('mv -v Results/mean.fav.zm.h5 Results/mean.fav.zm.bak')
      endif
      !
      call H5WriteArray(ro_zm,im,jm,'ro','Results/mean.fav.zm.h5')
      call H5WriteArray(u1_zm,im,jm,'u1','Results/mean.fav.zm.h5')
      call H5WriteArray(u2_zm,im,jm,'u2','Results/mean.fav.zm.h5')
      call H5WriteArray(u3_zm,im,jm,'u3','Results/mean.fav.zm.h5')
      call H5WriteArray( p_zm,im,jm, 'p','Results/mean.fav.zm.h5')
      call H5WriteArray( t_zm,im,jm, 't','Results/mean.fav.zm.h5')
      !
      call writetecbin('Results/teczm.plt',x(:,:,0),'x',y(:,:,0),'y',  &
                                       u1_zm,'u',u2_zm,'v',u3_zm,'w',  &
                                   ro_zm,'ro', p_zm,'p', t_zm,'t',im,jm)
      !
    endif
    !
    u1_m=u1_m/ro_m
    u2_m=u2_m/ro_m
    u3_m=u3_m/ro_m
     t_m= t_m/ro_m
    !
     p_m= p_m/real(nsamples,8)
    ro_m=ro_m/real(nsamples,8)
    ! !
    inquire(file='Results/mean.fav.h5',exist=filexists)
    if(filexists) then
      call system('mv -v Results/mean.fav.h5 Results/mean.fav.bak')
    endif
    !
    call H5WriteArray(ro_m,im,jm,km,'ro','Results/mean.fav.h5')
    call H5WriteArray(u1_m,im,jm,km,'u1','Results/mean.fav.h5')
    call H5WriteArray(u2_m,im,jm,km,'u2','Results/mean.fav.h5')
    call H5WriteArray(u3_m,im,jm,km,'u3','Results/mean.fav.h5')
    call H5WriteArray( p_m,im,jm,km, 'p','Results/mean.fav.h5')
    call H5WriteArray( t_m,im,jm,km, 't','Results/mean.fav.h5')
    !
  end subroutine meanflow_fav
  subroutine meanflow_rey(nsta,nend)
    !
    use commvardefine, only: im,jm,km,ro_m,u1_m,u2_m,u3_m,p_m,t_m,     &
                             ro_zm,u1_zm,u2_zm,u3_zm,p_zm,t_zm,        &
                             ro_xzm,u1_xzm,u2_xzm,u3_xzm,p_xzm,t_xzm,  &
                             lihomo,ljhomo,lkhomo,x,y
    !
    use h5readwrite
    use writetec
    use basicfunction, only: integration
    !
    integer :: nsta,nend,n
    !
    integer :: nsamples
    integer :: i,j,k
    real(8) :: time
    logical :: filexists=.false.
    character(len=4) sname
    character(len=23) fname
    !
    real(8),allocatable,dimension(:,:,:) :: ro,u1,u2,u3,p,t
    !
    allocate(ro(0:im,0:jm,0:km),u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km), &
             u3(0:im,0:jm,0:km),p(0:im,0:jm,0:km), t(0:im,0:jm,0:km)   )
    allocate(ro_m(0:im,0:jm,0:km),u1_m(0:im,0:jm,0:km),                &
             u2_m(0:im,0:jm,0:km),u3_m(0:im,0:jm,0:km),                &
              p_m(0:im,0:jm,0:km), t_m(0:im,0:jm,0:km)                 )
    !
    if(nsta>=0 .and. nsta<=99999) then
      !
      ro_m=0.d0
      u1_m=0.d0
      u2_m=0.d0
      u3_m=0.d0
       p_m=0.d0
       t_m=0.d0
      do n=nsta,nend
        !
        write(sname,'(i4.4)')n
        fname=outfolder//'flowfield'//sname//'.h5'
        call H5ReadArray(time,'time',fname)
        write(*,'(A,F10.5)')'  ** << flow field  at t=',time
        !
        ! call H5ReadArray(ro,im,jm,km,'ro',fname)
        call H5ReadArray(u1,im,jm,km,'u1',fname)
        call H5ReadArray(u2,im,jm,km,'u2',fname)
        call H5ReadArray(u3,im,jm,km,'u3',fname)
        ! call H5ReadArray( p,im,jm,km, 'p',fname)S
        !
        ! ro_m=ro_m+ro      
        u1_m=u1_m+u1
        u2_m=u2_m+u2
        u3_m=u3_m+u3
         ! p_m= p_m+ p
         ! t_m= t_m+ t
        !
      enddo

      !
      nsamples=nend-nsta+1
      !
    else
      !
      call H5ReadArray(nsamples,'nsamples',outfolder//'meanflow.h5')
      write(*,'(A,I5)')'  ** Total number of temporal samples:',nsamples
      !
      call H5ReadArray(ro_m,im,jm,km,'rom',outfolder//'meanflow.h5')
      call H5ReadArray(u1_m,im,jm,km,'u1m',outfolder//'meanflow.h5')
      call H5ReadArray(u2_m,im,jm,km,'u2m',outfolder//'meanflow.h5')
      call H5ReadArray(u3_m,im,jm,km,'u3m',outfolder//'meanflow.h5')
      call H5ReadArray( p_m,im,jm,km, 'pm',outfolder//'meanflow.h5')
      call H5ReadArray( t_m,im,jm,km, 'tm',outfolder//'meanflow.h5')
      !
    endif
    !
    if(lkhomo .and. lihomo) then
      allocate( ro_xzm(0:jm),u1_xzm(0:jm),u2_xzm(0:jm),u3_xzm(0:jm),   &
                 p_xzm(0:jm), t_xzm(0:jm) )
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j)
      !
      !$OMP DO
      do j=0,jm
        ro_xzm(j)=integration(ro_m(:,j,:))
        u1_xzm(j)=integration(u1_m(:,j,:))
        u2_xzm(j)=integration(u2_m(:,j,:))
        u3_xzm(j)=integration(u3_m(:,j,:))
         p_xzm(j)=integration( p_m(:,j,:))
         t_xzm(j)=integration( t_m(:,j,:))
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      print*,' ** averaged in i and k direction'
      !
      u1_xzm=u1_xzm/real(nsamples*im*km,8)
      u2_xzm=u2_xzm/real(nsamples*im*km,8)
      u3_xzm=u3_xzm/real(nsamples*im*km,8)
      !
       t_xzm= t_xzm/real(nsamples*im*km,8)
      !
      ro_xzm=ro_xzm/real(nsamples*im*km,8)
       p_xzm= p_xzm/real(nsamples*im*km,8)
      !
      inquire(file='Results/mean.rey.xzm.h5',exist=filexists)
      if(filexists) then
        call system('mv -v Results/mean.rey.xzm.h5 Results/mean.rey.zm.bak')
      endif
      !
      call H5WriteArray(ro_xzm,jm,'ro','Results/mean.rey.xzm.h5')
      call H5WriteArray(u1_xzm,jm,'u1','Results/mean.rey.xzm.h5')
      call H5WriteArray(u2_xzm,jm,'u2','Results/mean.rey.xzm.h5')
      call H5WriteArray(u3_xzm,jm,'u3','Results/mean.rey.xzm.h5')
      call H5WriteArray( p_xzm,jm, 'p','Results/mean.rey.xzm.h5')
      call H5WriteArray( t_xzm,jm, 't','Results/mean.rey.xzm.h5')
      !
      open(18,file='Results/profile.dat')
      write(18,"(6(1X,A15))")'y','ro','u','v','T','P'
      write(18,"(6(1X,E15.7E3))")(y(0,j,0),ro_xzm(j),u1_xzm(j),        &
                                  u2_xzm(j),t_xzm(j),p_xzm(j),j=0,jm)
      close(18)
      print*,' << profile.dat ... done. '
      !
    endif
    !
    if(lkhomo) then
      ! average in z direction
      !
      allocate(ro_zm(0:im,0:jm),u1_zm(0:im,0:jm),u2_zm(0:im,0:jm),     &
               u3_zm(0:im,0:jm), p_zm(0:im,0:jm), t_zm(0:im,0:jm))
      
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
      !
      !$OMP DO
      do j=0,jm
      do i=0,im
        ro_zm(i,j)=integration(ro_m(i,j,:))
        u1_zm(i,j)=integration(u1_m(i,j,:))
        u2_zm(i,j)=integration(u2_m(i,j,:))
        u3_zm(i,j)=integration(u3_m(i,j,:))
         p_zm(i,j)=integration( p_m(i,j,:))
         t_zm(i,j)=integration( t_m(i,j,:))
      enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      if(km==0) then
      u1_zm=u1_zm/real(nsamples,8)
      u2_zm=u2_zm/real(nsamples,8)
      u3_zm=u3_zm/real(nsamples,8)
      !
       t_zm= t_zm/real(nsamples,8)
      !
      ro_zm=ro_zm/real(nsamples,8)
       p_zm= p_zm/real(nsamples,8)
      else
      u1_zm=u1_zm/real(nsamples*km,8)
      u2_zm=u2_zm/real(nsamples*km,8)
      u3_zm=u3_zm/real(nsamples*km,8)
      !
       t_zm= t_zm/real(nsamples*km,8)
      !
      ro_zm=ro_zm/real(nsamples*km,8)
       p_zm= p_zm/real(nsamples*km,8)
      endif
      !
      inquire(file='Results/mean.rey.zm.h5',exist=filexists)
      if(filexists) then
        call system('mv -v Results/mean.rey.zm.h5 Results/mean.rey.zm.bak')
      endif
      !
      call H5WriteArray(ro_zm,im,jm,'ro','Results/mean.rey.zm.h5')
      call H5WriteArray(u1_zm,im,jm,'u1','Results/mean.rey.zm.h5')
      call H5WriteArray(u2_zm,im,jm,'u2','Results/mean.rey.zm.h5')
      call H5WriteArray(u3_zm,im,jm,'u3','Results/mean.rey.zm.h5')
      call H5WriteArray( p_zm,im,jm, 'p','Results/mean.rey.zm.h5')
      call H5WriteArray( t_zm,im,jm, 't','Results/mean.rey.zm.h5')
      !
      call writetecbin('Results/teczm.plt',x(:,:,0),'x',y(:,:,0),'y',  &
                                       u1_zm,'u',u2_zm,'v',u3_zm,'w',  &
                                   ro_zm,'ro', p_zm,'p', t_zm,'t',im,jm)
      !
    endif
    !
    u1_m=u1_m/real(nsamples,8)
    u2_m=u2_m/real(nsamples,8)
    u3_m=u3_m/real(nsamples,8)
     t_m= t_m/real(nsamples,8)
    !
     p_m= p_m/real(nsamples,8)
    ro_m=ro_m/real(nsamples,8)
    ! !
    inquire(file='Results/mean.rey.h5',exist=filexists)
    if(filexists) then
      call system('mv -v Results/mean.rey.h5 Results/mean.rey.bak')
    endif
    !
    call H5WriteArray(ro_m,im,jm,km,'ro','Results/mean.rey.h5')
    call H5WriteArray(u1_m,im,jm,km,'u1','Results/mean.rey.h5')
    call H5WriteArray(u2_m,im,jm,km,'u2','Results/mean.rey.h5')
    call H5WriteArray(u3_m,im,jm,km,'u3','Results/mean.rey.h5')
    call H5WriteArray( p_m,im,jm,km, 'p','Results/mean.rey.h5')
    call H5WriteArray( t_m,im,jm,km, 't','Results/mean.rey.h5')
    !
  end subroutine meanflow_rey
  !+-------------------------------------------------------------------+
  !| The end of the subroutine meanflow.                               |
  !+-------------------------------------------------------------------+
  !
  subroutine highorder_rey(nsta,nend)
    !
    use commvardefine, only: im,jm,km,ro_m,u1_m,u2_m,u3_m,p_m,t_m,     &
                             ro_zm,u1_zm,u2_zm,u3_zm,p_zm,t_zm,        &
                             ro_xzm,u1_xzm,u2_xzm,u3_xzm,p_xzm,t_xzm,  &
                             lihomo,ljhomo,lkhomo,x,y
    !
    use h5readwrite
    use writetec
    use basicfunction, only: integration
    !
    integer :: nsta,nend,n
    !
    integer :: nsamples
    integer :: i,j,k
    real(8) :: time
    logical :: filexists=.false.
    character(len=4) sname
    character(len=23) fname
    !
    real(8),allocatable,dimension(:,:,:) :: ro,u1,u2,u3,p,t
    !
    real(8),allocatable,dimension(:,:) :: u11_zm,u12_zm,u13_zm,u22_zm, &
                u23_zm,u33_zm,pu1_zm,pu2_zm,pu3_zm,pp_zm
    real(8) :: uf,vf,wf,pf
    !
    if(nsta>=0 .and. nsta<=99999) then
      !
      allocate(ro_zm(0:im,0:jm),u1_zm(0:im,0:jm),u2_zm(0:im,0:jm),  &
               u3_zm(0:im,0:jm),p_zm(0:im,0:jm),t_zm(0:im,0:jm))
      !
      call H5ReadArray(u1_zm,im,jm,'u1','Results/mean.rey.zm.h5')
      call H5ReadArray(u2_zm,im,jm,'u2','Results/mean.rey.zm.h5')
      call H5ReadArray(u3_zm,im,jm,'u3','Results/mean.rey.zm.h5')
      call H5ReadArray( p_zm,im,jm, 'p','Results/mean.rey.zm.h5')
      !
      allocate(u11_zm(0:im,0:jm),u12_zm(0:im,0:jm),                    &
               u22_zm(0:im,0:jm),u13_zm(0:im,0:jm),                    &
               u23_zm(0:im,0:jm),u33_zm(0:im,0:jm),                    &
               pu1_zm(0:im,0:jm),pu2_zm(0:im,0:jm),                    &
               pu3_zm(0:im,0:jm), pp_zm(0:im,0:jm))
      u11_zm=0.d0
      u12_zm=0.d0
      u13_zm=0.d0
      u22_zm=0.d0
      u23_zm=0.d0
      u33_zm=0.d0
      pu1_zm=0.d0
      pu2_zm=0.d0
      pu3_zm=0.d0
      pp_zm=0.d0
      !
      allocate(u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),u3(0:im,0:jm,0:km), &
               p(0:im,0:jm,0:km))
      do n=nsta,nend
        !
        write(sname,'(i4.4)')n
        fname=outfolder//'flowfield'//sname//'.h5'
        call H5ReadArray(time,'time',fname)
        write(*,'(A,F10.5)')'  ** << flow field  at t=',time
        !
        call H5ReadArray(u1,im,jm,km,'u1',fname)
        call H5ReadArray(u2,im,jm,km,'u2',fname)
        call H5ReadArray(u3,im,jm,km,'u3',fname)
        call H5ReadArray( p,im,jm,km, 'p',fname)
        !
        do j=0,jm
        do i=0,im
          !
          do k=1,km
            uf=u1(i,j,k)-u1_zm(i,j)
            vf=u2(i,j,k)-u2_zm(i,j)
            wf=u3(i,j,k)-u3_zm(i,j)
            pf= p(i,j,k)- p_zm(i,j)
            !
            u11_zm(i,j)=u11_zm(i,j)+uf*uf
            u12_zm(i,j)=u12_zm(i,j)+uf*vf
            u13_zm(i,j)=u13_zm(i,j)+uf*wf
            u22_zm(i,j)=u22_zm(i,j)+vf*vf
            u23_zm(i,j)=u23_zm(i,j)+vf*wf
            u33_zm(i,j)=u33_zm(i,j)+wf*wf
            pu1_zm(i,j)=pu1_zm(i,j)+pf*uf
            pu2_zm(i,j)=pu2_zm(i,j)+pf*vf
            pu3_zm(i,j)=pu3_zm(i,j)+pf*wf
             pp_zm(i,j)= pp_zm(i,j)+pf*pf
          enddo
          !
        enddo
        enddo
        !
      enddo
      !
      nsamples=nend-nsta+1
      !
    endif
    !
    u11_zm=u11_zm/real(nsamples*km,8)
    u12_zm=u12_zm/real(nsamples*km,8)
    u13_zm=u13_zm/real(nsamples*km,8)
    u22_zm=u22_zm/real(nsamples*km,8)
    u23_zm=u23_zm/real(nsamples*km,8)
    u33_zm=u33_zm/real(nsamples*km,8)
    pu1_zm=pu1_zm/real(nsamples*km,8)
    pu2_zm=pu2_zm/real(nsamples*km,8)
    pu3_zm=pu3_zm/real(nsamples*km,8)
     pp_zm= pp_zm/real(nsamples*km,8)
    !
    call H5WriteArray(u11_zm,im,jm,'u11','Results/2order.rey.zm.h5')
    call H5WriteArray(u22_zm,im,jm,'u22','Results/2order.rey.zm.h5')
    call H5WriteArray(u33_zm,im,jm,'u33','Results/2order.rey.zm.h5')
    call H5WriteArray(u12_zm,im,jm,'u12','Results/2order.rey.zm.h5')
    call H5WriteArray(u13_zm,im,jm,'u13','Results/2order.rey.zm.h5')
    call H5WriteArray(u23_zm,im,jm,'u23','Results/2order.rey.zm.h5')
    call H5WriteArray( pp_zm,im,jm, 'pp','Results/2order.rey.zm.h5')
    call H5WriteArray(pu1_zm,im,jm,'pu1','Results/2order.rey.zm.h5')
    call H5WriteArray(pu2_zm,im,jm,'pu2','Results/2order.rey.zm.h5')
    call H5WriteArray(pu3_zm,im,jm,'pu3','Results/2order.rey.zm.h5')
    !
  end subroutine highorder_rey
  !
  subroutine betchovsta(nsta,nend)
    !
    use commvardefine, only: im,jm,km,ro_m,u1_m,u2_m,u3_m,p_m,t_m,     &
                             ro_zm,u1_zm,u2_zm,u3_zm,p_zm,t_zm,        &
                             ro_xzm,u1_xzm,u2_xzm,u3_xzm,p_xzm,t_xzm,  &
                             lihomo,ljhomo,lkhomo,x,y,z
    !
    use h5readwrite
    use writetec
    use gradsolver
    use basicfunction, only: integration
    !
    integer :: nsta,nend,n
    !
    integer :: nsamples
    integer :: i,j,k
    real(8) :: time
    logical :: filexists=.false.
    character(len=4) sname
    character(len=23) fname
    !
    real(8),allocatable,dimension(:,:,:) :: ro,u1,u2,u3,p,t
    real(8),allocatable,dimension(:,:,:,:) :: du1,du2,du3
    !
    real(8),allocatable,dimension(:,:) :: a1111,a2222,a3333, &
                                          a111111,a222222,a333333, &
                                          a121211,a212122,a131311, &
                                          a313133,a323233,a232322, &
                                          a111122,a111133,a222211, &
                                          a222233,a333311,a333322, &
                                          a113232,a112323,a223131, &
                                          a221313,a331212,a332121
    real(8) :: uf,vf,wf,pf,du11,du12,du13,du21,du22,du23,du31,du32,du33
    !
    if(nsta>=0 .and. nsta<=99999) then
      !
      allocate(ro_zm(0:im,0:jm),u1_zm(0:im,0:jm),u2_zm(0:im,0:jm),  &
               u3_zm(0:im,0:jm),p_zm(0:im,0:jm),t_zm(0:im,0:jm))
      !
      call H5ReadArray(u1_zm,im,jm,'u1','Results/mean.rey.zm.h5')
      call H5ReadArray(u2_zm,im,jm,'u2','Results/mean.rey.zm.h5')
      call H5ReadArray(u3_zm,im,jm,'u3','Results/mean.rey.zm.h5')
      !
      allocate(a1111(0:im,0:jm),a2222(0:im,0:jm),a3333(0:im,0:jm), &
               a111111(0:im,0:jm),a222222(0:im,0:jm),a333333(0:im,0:jm), &
               a121211(0:im,0:jm),a212122(0:im,0:jm),a131311(0:im,0:jm), &
               a313133(0:im,0:jm),a323233(0:im,0:jm),a232322(0:im,0:jm), &
               a111122(0:im,0:jm),a111133(0:im,0:jm),a222211(0:im,0:jm), &
               a222233(0:im,0:jm),a333311(0:im,0:jm),a333322(0:im,0:jm), &
               a113232(0:im,0:jm),a112323(0:im,0:jm),a223131(0:im,0:jm), &
               a221313(0:im,0:jm),a331212(0:im,0:jm),a332121(0:im,0:jm) )
      !
      allocate(du1(3,0:im,0:jm,0:km),du2(3,0:im,0:jm,0:km),du3(3,0:im,0:jm,0:km))
      
      allocate(u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),u3(0:im,0:jm,0:km) )
      !
      a1111   = 0.d0
      a2222   = 0.d0
      a3333   = 0.d0
      a111111 = 0.d0
      a222222 = 0.d0
      a333333 = 0.d0
      a121211 = 0.d0
      a212122 = 0.d0
      a131311 = 0.d0
      a313133 = 0.d0
      a323233 = 0.d0
      a232322 = 0.d0
      a111122 = 0.d0
      a111133 = 0.d0 
      a222211 = 0.d0 
      a222233 = 0.d0 
      a333311 = 0.d0
      a333322 = 0.d0
      a113232 = 0.d0 
      a112323 = 0.d0 
      a223131 = 0.d0 
      a221313 = 0.d0 
      a331212 = 0.d0
      a332121 = 0.d0
      !
      do n=nsta,nend
        !
        write(sname,'(i4.4)')n
        fname=outfolder//'flowfield'//sname//'.h5'
        call H5ReadArray(time,'time',fname)
        write(*,'(A,F10.5)')'  ** << flow field  at t=',time
        !
        call H5ReadArray(u1,im,jm,km,'u1',fname)
        call H5ReadArray(u2,im,jm,km,'u2',fname)
        call H5ReadArray(u3,im,jm,km,'u3',fname)
        !
        do k=0,km
          u1(:,:,k)=u1(:,:,k)-u1_zm(:,:)
          u2(:,:,k)=u2(:,:,k)-u2_zm(:,:)
          u3(:,:,k)=u3(:,:,k)-u3_zm(:,:)
        enddo
        !
        if(n==nsta) then
          du1=grad(u1,x,y,z)
        else
          du1=grad(u1)
        endif
        du2=grad(u2)
        du3=grad(u3)
        !
        do j=0,jm
        do i=0,im
          !
          do k=1,km
            du11=du1(1,i,j,k); du21=du2(1,i,j,k); du31=du3(1,i,j,k)
            du12=du1(2,i,j,k); du22=du2(2,i,j,k); du32=du3(2,i,j,k)
            du13=du1(3,i,j,k); du23=du2(3,i,j,k); du33=du3(3,i,j,k)
            !
            a1111(i,j)   = a1111(i,j)   + du11*du11
            a2222(i,j)   = a2222(i,j)   + du22*du22
            a3333(i,j)   = a3333(i,j)   + du33*du33
            !
            a111111(i,j) = a111111(i,j) + du11*du11*du11
            a222222(i,j) = a222222(i,j) + du22*du22*du22
            a333333(i,j) = a333333(i,j) + du33*du33*du33
            !
            a121211(i,j) = a121211(i,j) + du12*du12*du11 
            a212122(i,j) = a121211(i,j) + du21*du21*du22
            a131311(i,j) = a121211(i,j) + du13*du13*du11
            a313133(i,j) = a121211(i,j) + du31*du31*du33
            a323233(i,j) = a121211(i,j) + du32*du32*du33
            a232322(i,j) = a121211(i,j) + du23*du23*du22
            !
            a111122(i,j) = a111122(i,j) + du11*du11*du22
            a111133(i,j) = a111133(i,j) + du11*du11*du33 
            a222211(i,j) = a222211(i,j) + du22*du22*du11 
            a222233(i,j) = a222233(i,j) + du22*du22*du33 
            a333311(i,j) = a333311(i,j) + du33*du33*du11
            a333322(i,j) = a333322(i,j) + du33*du33*du22
            !
            a113232(i,j) = a113232(i,j) + du11*du32*du32 
            a112323(i,j) = a112323(i,j) + du11*du23*du23 
            a223131(i,j) = a223131(i,j) + du22*du31*du31 
            a221313(i,j) = a221313(i,j) + du22*du13*du13 
            a331212(i,j) = a331212(i,j) + du33*du12*du12
            a332121(i,j) = a332121(i,j) + du33*du21*du21
          enddo
          !
        enddo
        enddo
        !
        nsamples=nend-nsta+1
        !
      enddo
      !
    endif
    !
    a1111   = a1111 /real(nsamples*km,8)
    a2222   = a2222 /real(nsamples*km,8)
    a3333   = a3333 /real(nsamples*km,8)
    a111111 = a111111 /real(nsamples*km,8)
    a222222 = a222222 /real(nsamples*km,8)
    a333333 = a333333 /real(nsamples*km,8)
    a121211 = a121211 /real(nsamples*km,8)
    a212122 = a121211 /real(nsamples*km,8)
    a131311 = a121211 /real(nsamples*km,8)
    a313133 = a121211 /real(nsamples*km,8)
    a323233 = a121211 /real(nsamples*km,8)
    a232322 = a121211 /real(nsamples*km,8)
    a111122 = a111122 /real(nsamples*km,8)
    a111133 = a111133 /real(nsamples*km,8) 
    a222211 = a222211 /real(nsamples*km,8) 
    a222233 = a222233 /real(nsamples*km,8) 
    a333311 = a333311 /real(nsamples*km,8)
    a333322 = a333322 /real(nsamples*km,8)
    a113232 = a113232 /real(nsamples*km,8) 
    a112323 = a112323 /real(nsamples*km,8) 
    a223131 = a223131 /real(nsamples*km,8) 
    a221313 = a221313 /real(nsamples*km,8) 
    a331212 = a331212 /real(nsamples*km,8)
    a332121 = a332121 /real(nsamples*km,8)
    !
    call H5WriteArray(a1111  ,im,jm,'A1111',  'Results/betchov_a.zm.h5')
    call H5WriteArray(a2222  ,im,jm,'A2222',  'Results/betchov_a.zm.h5')
    call H5WriteArray(a3333  ,im,jm,'A3333',  'Results/betchov_a.zm.h5')
    call H5WriteArray(a111111,im,jm,'A111111','Results/betchov_a.zm.h5')
    call H5WriteArray(a222222,im,jm,'A222222','Results/betchov_a.zm.h5')
    call H5WriteArray(a333333,im,jm,'A333333','Results/betchov_a.zm.h5')
    call H5WriteArray(a121211,im,jm,'A121211','Results/betchov_a.zm.h5')
    call H5WriteArray(a212122,im,jm,'A212122','Results/betchov_a.zm.h5')
    call H5WriteArray(a131311,im,jm,'A131311','Results/betchov_a.zm.h5')
    call H5WriteArray(a313133,im,jm,'A313133','Results/betchov_a.zm.h5')
    call H5WriteArray(a323233,im,jm,'A323233','Results/betchov_a.zm.h5')
    call H5WriteArray(a232322,im,jm,'A232322','Results/betchov_a.zm.h5')
    call H5WriteArray(a111122,im,jm,'A111122','Results/betchov_a.zm.h5')
    call H5WriteArray(a111133,im,jm,'A111133','Results/betchov_a.zm.h5')
    call H5WriteArray(a222211,im,jm,'A222211','Results/betchov_a.zm.h5')
    call H5WriteArray(a222233,im,jm,'A222233','Results/betchov_a.zm.h5')
    call H5WriteArray(a333311,im,jm,'A333311','Results/betchov_a.zm.h5')
    call H5WriteArray(a333322,im,jm,'A333322','Results/betchov_a.zm.h5')
    call H5WriteArray(a113232,im,jm,'A113232','Results/betchov_a.zm.h5')
    call H5WriteArray(a112323,im,jm,'A112323','Results/betchov_a.zm.h5')
    call H5WriteArray(a223131,im,jm,'A223131','Results/betchov_a.zm.h5')
    call H5WriteArray(a221313,im,jm,'A221313','Results/betchov_a.zm.h5')
    call H5WriteArray(a331212,im,jm,'A331212','Results/betchov_a.zm.h5')
    call H5WriteArray(a332121,im,jm,'A332121','Results/betchov_a.zm.h5')
    !
  end subroutine betchovsta
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to process budget data.                   |
  !+-------------------------------------------------------------------+
  subroutine budgetnorm
    !
    use commvardefine, only: im,jm,km,outfolder,lihomo,ljhomo,lkhomo
    use basicfunction, only: integration
    use h5readwrite
    !
    real(8) :: rsamples
    integer :: nsamples,i,j
    logical :: filexists
    real(8),allocatable,dimension(:,:,:) :: pu1,pu2,pu3,u1re,u2re,u3re,&
                                            disspa,predil,visdif1,     &
                                            visdif2,visdif3,           &
                                            sgmam11m,sgmam22m,sgmam33m,&
                                            sgmam12m,sgmam13m,sgmam23m
    real(8),allocatable,dimension(:,:) :: u1re_zm,u2re_zm,u3re_zm,     &
                                          pu1_zm,pu2_zm,pu3_zm,        &
                                          disspa_zm,predil_zm,         &
                                     visdif1_zm,visdif2_zm,visdif3_zm, &
                                        sgmam11_zm,sgmam22_zm,         &
                                          sgmam33_zm,sgmam12_zm,       &
                                          sgmam13_zm,sgmam23_zm
    real(8),allocatable,dimension(:) :: u1re_xzm,u2re_xzm,u3re_xzm,    &
                                          pu1_xzm,pu2_xzm,pu3_xzm,     &
                                          disspa_xzm,predil_xzm,       &
                                          visdif1_xzm,visdif2_xzm,     &
                                          visdif3_xzm,                 &
                                          sgmam11_xzm,sgmam22_xzm,     &
                                          sgmam33_xzm,sgmam12_xzm,     &
                                          sgmam13_xzm,sgmam23_xzm
    !
    allocate(pu1(0:im,0:jm,0:km),pu2(0:im,0:jm,0:km),                  &
             pu3(0:im,0:jm,0:km),u1re(0:im,0:jm,0:km),                 &
             u2re(0:im,0:jm,0:km),u3re(0:im,0:jm,0:km),                &
             disspa(0:im,0:jm,0:km),predil(0:im,0:jm,0:km)) 
    allocate(visdif1(0:im,0:jm,0:km),visdif2(0:im,0:jm,0:km),          &
             visdif3(0:im,0:jm,0:km))
    !
    call H5ReadArray(nsamples,'nsamples',outfolder//'budget.h5')
    !
    call H5ReadArray(disspa,im,jm,km,'disspa',outfolder//'budget.h5')
    call H5ReadArray(predil,im,jm,km,'predil',outfolder//'budget.h5')
    !
    call H5ReadArray(pu1,im,jm,km,'pu1',outfolder//'budget.h5')
    call H5ReadArray(pu2,im,jm,km,'pu2',outfolder//'budget.h5')
    call H5ReadArray(pu3,im,jm,km,'pu3',outfolder//'budget.h5')
    call H5ReadArray(u1re,im,jm,km,'u1rem',outfolder//'budget.h5')
    call H5ReadArray(u2re,im,jm,km,'u2rem',outfolder//'budget.h5')
    call H5ReadArray(u3re,im,jm,km,'u3rem',outfolder//'budget.h5')
    call H5ReadArray(visdif1,im,jm,km,'visdif1',outfolder//'budget.h5')
    call H5ReadArray(visdif2,im,jm,km,'visdif2',outfolder//'budget.h5')
    call H5ReadArray(visdif3,im,jm,km,'visdif3',outfolder//'budget.h5')
    !
    allocate(sgmam11m(0:im,0:jm,0:km),sgmam22m(0:im,0:jm,0:km),        &
             sgmam33m(0:im,0:jm,0:km),sgmam12m(0:im,0:jm,0:km),        &
             sgmam13m(0:im,0:jm,0:km),sgmam23m(0:im,0:jm,0:km)         )
    call H5ReadArray(sgmam11m,im,jm,km,'sgmam11',outfolder//'budget.h5')
    call H5ReadArray(sgmam22m,im,jm,km,'sgmam22',outfolder//'budget.h5')
    call H5ReadArray(sgmam33m,im,jm,km,'sgmam33',outfolder//'budget.h5')
    call H5ReadArray(sgmam12m,im,jm,km,'sgmam12',outfolder//'budget.h5')
    call H5ReadArray(sgmam13m,im,jm,km,'sgmam13',outfolder//'budget.h5')
    call H5ReadArray(sgmam23m,im,jm,km,'sgmam23',outfolder//'budget.h5')
    !
    rsamples=real(nsamples,8)
    print*,' ** number of samples: ',rsamples
    !
    u1re=u1re/rsamples
    u2re=u2re/rsamples
    u3re=u3re/rsamples
    !
    pu1= pu1/rsamples
    pu2= pu2/rsamples
    pu3= pu3/rsamples
    !
    disspa= disspa/rsamples
    predil= predil/rsamples
    !
    visdif1=visdif1/rsamples
    visdif2=visdif2/rsamples
    visdif3=visdif3/rsamples
    !
    sgmam11m=sgmam11m/rsamples
    sgmam22m=sgmam22m/rsamples
    sgmam33m=sgmam33m/rsamples
    sgmam12m=sgmam12m/rsamples
    sgmam13m=sgmam13m/rsamples
    sgmam23m=sgmam23m/rsamples
    !
    inquire(file='Results/budget.h5',exist=filexists)
    if(filexists) call system('mv -v Results/budget.h5 Results/budget.bak')
    !
    call H5WriteArray(disspa,im,jm,km,'disspa','Results/budget.h5')
    call H5WriteArray(predil,im,jm,km,'predil','Results/budget.h5')
    call H5WriteArray(pu1,im,jm,km,'pu1','Results/budget.h5')
    call H5WriteArray(pu2,im,jm,km,'pu2','Results/budget.h5')
    call H5WriteArray(pu3,im,jm,km,'pu3','Results/budget.h5')
    call H5WriteArray(u1re,im,jm,km,'u1rem','Results/budget.h5')
    call H5WriteArray(u2re,im,jm,km,'u2rem','Results/budget.h5')
    call H5WriteArray(u3re,im,jm,km,'u3rem','Results/budget.h5')
    call H5WriteArray(visdif1,im,jm,km,'visdif1','Results/budget.h5')
    call H5WriteArray(visdif2,im,jm,km,'visdif2','Results/budget.h5')
    call H5WriteArray(visdif3,im,jm,km,'visdif3','Results/budget.h5')
    !
    call H5WriteArray(sgmam11m,im,jm,km,'sgmam11','Results/budget.h5')
    call H5WriteArray(sgmam22m,im,jm,km,'sgmam22','Results/budget.h5')
    call H5WriteArray(sgmam33m,im,jm,km,'sgmam33','Results/budget.h5')
    call H5WriteArray(sgmam12m,im,jm,km,'sgmam12','Results/budget.h5')
    call H5WriteArray(sgmam13m,im,jm,km,'sgmam13','Results/budget.h5')
    call H5WriteArray(sgmam23m,im,jm,km,'sgmam23','Results/budget.h5')
    !
    if(lkhomo .and. lihomo) then
      !
      allocate(u1re_xzm(0:jm),u2re_xzm(0:jm),u3re_xzm(0:jm), &
               pu1_xzm(0:jm),pu2_xzm(0:jm),pu3_xzm(0:jm),    &
               disspa_xzm(0:jm),predil_xzm(0:jm),            &
               visdif1_xzm(0:jm),visdif2_xzm(0:jm),          &
               visdif3_xzm(0:jm))
      !
      allocate(sgmam11_xzm(0:jm),sgmam22_xzm(0:jm),        &
               sgmam33_xzm(0:jm),sgmam12_xzm(0:jm),        &
               sgmam13_xzm(0:jm),sgmam23_xzm(0:jm)         )
      !
      do j=0,jm
        !
        disspa_xzm(j)=integration(disspa(:,j,:))/real(im*km,8)
        predil_xzm(j)=integration(predil(:,j,:))/real(im*km,8)
        pu1_xzm(j)=integration(pu1(:,j,:))/real(im*km,8)
        pu2_xzm(j)=integration(pu2(:,j,:))/real(im*km,8)
        pu3_xzm(j)=integration(pu3(:,j,:))/real(im*km,8)
        u1re_xzm(j)=integration(u1re(:,j,:))/real(im*km,8)
        u2re_xzm(j)=integration(u2re(:,j,:))/real(im*km,8)
        u3re_xzm(j)=integration(u3re(:,j,:))/real(im*km,8)
        visdif1_xzm(j)=integration(visdif1(:,j,:))/real(im*km,8)
        visdif2_xzm(j)=integration(visdif2(:,j,:))/real(im*km,8)
        visdif3_xzm(j)=integration(visdif3(:,j,:))/real(im*km,8)
        !
        sgmam11_xzm(j)=integration(sgmam11m(:,j,:))/real(im*km,8)
        sgmam22_xzm(j)=integration(sgmam22m(:,j,:))/real(im*km,8)
        sgmam33_xzm(j)=integration(sgmam33m(:,j,:))/real(im*km,8)
        sgmam12_xzm(j)=integration(sgmam12m(:,j,:))/real(im*km,8)
        sgmam13_xzm(j)=integration(sgmam13m(:,j,:))/real(im*km,8)
        sgmam23_xzm(j)=integration(sgmam23m(:,j,:))/real(im*km,8)
      enddo
      !
      inquire(file='Results/budget.xzm.h5',exist=filexists)
      if(filexists) call system('mv -v Results/budget.xzm.h5 Results/budget.xzm.bak')
      !
      call H5WriteArray(disspa_xzm,jm,'disspa','Results/budget.xzm.h5')
      call H5WriteArray(predil_xzm,jm,'predil','Results/budget.xzm.h5')
      call H5WriteArray(pu1_xzm,jm,'pu1','Results/budget.xzm.h5')
      call H5WriteArray(pu2_xzm,jm,'pu2','Results/budget.xzm.h5')
      call H5WriteArray(pu3_xzm,jm,'pu3','Results/budget.xzm.h5')
      call H5WriteArray(u1re_xzm,jm,'u1rem','Results/budget.xzm.h5')
      call H5WriteArray(u2re_xzm,jm,'u2rem','Results/budget.xzm.h5')
      call H5WriteArray(u3re_xzm,jm,'u3rem','Results/budget.xzm.h5')
      call H5WriteArray(visdif1_xzm,jm,'visdif1','Results/budget.xzm.h5')
      call H5WriteArray(visdif2_xzm,jm,'visdif2','Results/budget.xzm.h5')
      call H5WriteArray(visdif3_xzm,jm,'visdif3','Results/budget.xzm.h5')
      !
      call H5WriteArray(sgmam11_xzm,jm,'sgmam11','Results/budget.xzm.h5')
      call H5WriteArray(sgmam22_xzm,jm,'sgmam22','Results/budget.xzm.h5')
      call H5WriteArray(sgmam33_xzm,jm,'sgmam33','Results/budget.xzm.h5')
      call H5WriteArray(sgmam12_xzm,jm,'sgmam12','Results/budget.xzm.h5')
      call H5WriteArray(sgmam13_xzm,jm,'sgmam13','Results/budget.xzm.h5')
      call H5WriteArray(sgmam23_xzm,jm,'sgmam23','Results/budget.xzm.h5')
      !
      print*,' ** normlising budets in k direction'
      !
      deallocate(u1re_xzm,u2re_xzm,u3re_xzm,pu1_xzm,pu2_xzm,pu3_xzm,           &
                 disspa_xzm,predil_xzm,visdif1_xzm,visdif2_xzm,visdif3_xzm,   &
                 sgmam11_xzm,sgmam22_xzm,sgmam33_xzm,sgmam12_xzm,sgmam13_xzm, &
                 sgmam23_xzm)
      !
    endif
    !
    if(lkhomo) then
      allocate(u1re_zm(0:im,0:jm),u2re_zm(0:im,0:jm),u3re_zm(0:im,0:jm), &
               pu1_zm(0:im,0:jm),pu2_zm(0:im,0:jm),pu3_zm(0:im,0:jm),    &
               disspa_zm(0:im,0:jm),predil_zm(0:im,0:jm),                &
               visdif1_zm(0:im,0:jm),visdif2_zm(0:im,0:jm),visdif3_zm(0:im,0:jm))
      !
      allocate(sgmam11_zm(0:im,0:jm),sgmam22_zm(0:im,0:jm),        &
               sgmam33_zm(0:im,0:jm),sgmam12_zm(0:im,0:jm),        &
               sgmam13_zm(0:im,0:jm),sgmam23_zm(0:im,0:jm)         )
      !
      do j=0,jm
      do i=0,im
        !
        disspa_zm(i,j)=integration(disspa(i,j,:))/dble(km)
        predil_zm(i,j)=integration(predil(i,j,:))/dble(km)
        pu1_zm(i,j)=integration(pu1(i,j,:))/dble(km)
        pu2_zm(i,j)=integration(pu2(i,j,:))/dble(km)
        pu3_zm(i,j)=integration(pu3(i,j,:))/dble(km)
        u1re_zm(i,j)=integration(u1re(i,j,:))/dble(km)
        u2re_zm(i,j)=integration(u2re(i,j,:))/dble(km)
        u3re_zm(i,j)=integration(u3re(i,j,:))/dble(km)
        visdif1_zm(i,j)=integration(visdif1(i,j,:))/dble(km)
        visdif2_zm(i,j)=integration(visdif2(i,j,:))/dble(km)
        visdif3_zm(i,j)=integration(visdif3(i,j,:))/dble(km)
        !
        sgmam11_zm(i,j)=integration(sgmam11m(i,j,:))/dble(km)
        sgmam22_zm(i,j)=integration(sgmam22m(i,j,:))/dble(km)
        sgmam33_zm(i,j)=integration(sgmam33m(i,j,:))/dble(km)
        sgmam12_zm(i,j)=integration(sgmam12m(i,j,:))/dble(km)
        sgmam13_zm(i,j)=integration(sgmam13m(i,j,:))/dble(km)
        sgmam23_zm(i,j)=integration(sgmam23m(i,j,:))/dble(km)
      enddo
      enddo
      !
      inquire(file='Results/budget.zm.h5',exist=filexists)
      if(filexists) call system('mv -v Results/budget.zm.h5 Results/budget.zm.bak')
      !
      call H5WriteArray(disspa_zm,im,jm,'disspa','Results/budget.zm.h5')
      call H5WriteArray(predil_zm,im,jm,'predil','Results/budget.zm.h5')
      call H5WriteArray(pu1_zm,im,jm,'pu1','Results/budget.zm.h5')
      call H5WriteArray(pu2_zm,im,jm,'pu2','Results/budget.zm.h5')
      call H5WriteArray(pu3_zm,im,jm,'pu3','Results/budget.zm.h5')
      call H5WriteArray(u1re_zm,im,jm,'u1rem','Results/budget.zm.h5')
      call H5WriteArray(u2re_zm,im,jm,'u2rem','Results/budget.zm.h5')
      call H5WriteArray(u3re_zm,im,jm,'u3rem','Results/budget.zm.h5')
      call H5WriteArray(visdif1_zm,im,jm,'visdif1','Results/budget.zm.h5')
      call H5WriteArray(visdif2_zm,im,jm,'visdif2','Results/budget.zm.h5')
      call H5WriteArray(visdif3_zm,im,jm,'visdif3','Results/budget.zm.h5')
      !
      call H5WriteArray(sgmam11_zm,im,jm,'sgmam11','Results/budget.zm.h5')
      call H5WriteArray(sgmam22_zm,im,jm,'sgmam22','Results/budget.zm.h5')
      call H5WriteArray(sgmam33_zm,im,jm,'sgmam33','Results/budget.zm.h5')
      call H5WriteArray(sgmam12_zm,im,jm,'sgmam12','Results/budget.zm.h5')
      call H5WriteArray(sgmam13_zm,im,jm,'sgmam13','Results/budget.zm.h5')
      call H5WriteArray(sgmam23_zm,im,jm,'sgmam23','Results/budget.zm.h5')
      !
      print*,' ** normlising budets in k direction'
      !
      deallocate(u1re_zm,u2re_zm,u3re_zm,pu1_zm,pu2_zm,pu3_zm,           &
                 disspa_zm,predil_zm,visdif1_zm,visdif2_zm,visdif3_zm,   &
                 sgmam11_zm,sgmam22_zm,sgmam33_zm,sgmam12_zm,sgmam13_zm, &
                 sgmam23_zm)
    endif

    deallocate(pu1,pu2,pu3,u1re,u2re,u3re,disspa,predil,visdif1,       &
               visdif2,visdif3,sgmam11m,sgmam22m,sgmam33m,             &
               sgmam12m,sgmam13m,sgmam23m)
    !
  end subroutine budgetnorm
  !+-------------------------------------------------------------------+
  !| The end of the subroutine budgetnorm.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate statistics in budget.h5      |
  !+-------------------------------------------------------------------+
  subroutine TKEBudgetSta
    !
    use commvardefine, only: im,jm,km,lihomo,ljhomo,lkhomo,x,y,z,Reynolds
    use h5readwrite
    use readwrite, only: H5_ReadGrid
    use basicfunction, only: integration,miucal
    use gradsolver
    use writetec
    !
    integer :: nsamples,i,j,k
    !
    real(8),allocatable,dimension(:,:,:) :: pu1,pu2,pu3,u1re,u2re,u3re,&
                                            sigma11,sigma22,sigma33,   &
                                            sigma12,sigma13,sigma23,   &
                                            disspa,predil,visdif1,     &
                                            visdif2,visdif3,u1,u2,u3,  &
                                            p,t,ro,                    &
                                            usigma1,usigma2,usigma3,   &
                                            u11,u12,u13,u22,u23,u33,   &
                                            u111,u112,u113,u122,u222,  &
                                            u223,u133,u233,u333,tke
    real(8),allocatable,dimension(:,:,:) :: sgmam11m,sgmam22m,sgmam33m,&
                                            sgmam12m,sgmam13m,sgmam23m
    real(8),allocatable,dimension(:,:,:,:) :: du1re,du2re,du3re,       &
                                              du1,du2,du3,gengrad
    ! TKE budget terms
    real(8),allocatable,dimension(:,:,:) :: p_tran,dissip,pres_dila,   &
                                            visc_acce,pres_acce,       &
                                            visc_diff,convec,produc,   &
                                            turb_tran
    real(8),allocatable,dimension(:,:) :: p_tran_zm,dissip_zm,         &
                                          pres_dila_zm,visc_acce_zm,   &
                                          pres_acce_zm,visc_diff_zm,   &
                                          convec_zm,produc_zm,         &
                                          turb_tran_zm
    !
    real(8),allocatable,dimension(:,:,:) :: term1,term2,term3
    !
    real(8) :: miu,var1,var2,var3
    real(8) :: f11,f12,f13,f22,f23,f33,d11,d12,d13,d21,d22,d23,        &
               d31,d32,d33,dkk,g11,g12,g13,g22,g23,g33,e11,e12,e13,    &
               e21,e22,e23,e31,e32,e33,ekk,usg1,usg2,usg3
    logical :: filexists
    character(len=255) :: string
    !
    !+--------------------------+
    !| programme starts         |
    !+--------------------------+
    !
    allocate(pu1(0:im,0:jm,0:km),pu2(0:im,0:jm,0:km),                  &
             pu3(0:im,0:jm,0:km),u1re(0:im,0:jm,0:km),                 &
             u2re(0:im,0:jm,0:km),u3re(0:im,0:jm,0:km),                &
             disspa(0:im,0:jm,0:km),predil(0:im,0:jm,0:km)) 
    allocate(visdif1(0:im,0:jm,0:km),visdif2(0:im,0:jm,0:km),          &
             visdif3(0:im,0:jm,0:km))
    !
    allocate(sgmam11m(0:im,0:jm,0:km),sgmam22m(0:im,0:jm,0:km),        &
             sgmam33m(0:im,0:jm,0:km),sgmam12m(0:im,0:jm,0:km),        &
             sgmam13m(0:im,0:jm,0:km),sgmam23m(0:im,0:jm,0:km)         )
    !
    call H5ReadArray(disspa,im,jm,km,'disspa','Results/budget.h5')
    call H5ReadArray(predil,im,jm,km,'predil','Results/budget.h5')
    call H5ReadArray(pu1,im,jm,km,'pu1','Results/budget.h5')
    call H5ReadArray(pu2,im,jm,km,'pu2','Results/budget.h5')
    call H5ReadArray(pu3,im,jm,km,'pu3','Results/budget.h5')
    call H5ReadArray(u1re,im,jm,km,'u1rem','Results/budget.h5')
    call H5ReadArray(u2re,im,jm,km,'u2rem','Results/budget.h5')
    call H5ReadArray(u3re,im,jm,km,'u3rem','Results/budget.h5')
    call H5ReadArray(visdif1,im,jm,km,'visdif1','Results/budget.h5')
    call H5ReadArray(visdif2,im,jm,km,'visdif2','Results/budget.h5')
    call H5ReadArray(visdif3,im,jm,km,'visdif3','Results/budget.h5')
    !
    call H5ReadArray(sgmam11m,im,jm,km,'sgmam11','Results/budget.h5')
    call H5ReadArray(sgmam22m,im,jm,km,'sgmam22','Results/budget.h5')
    call H5ReadArray(sgmam33m,im,jm,km,'sgmam33','Results/budget.h5')
    call H5ReadArray(sgmam12m,im,jm,km,'sgmam12','Results/budget.h5')
    call H5ReadArray(sgmam13m,im,jm,km,'sgmam13','Results/budget.h5')
    call H5ReadArray(sgmam23m,im,jm,km,'sgmam23','Results/budget.h5')
    !
    allocate(gengrad(1:3,0:im,0:jm,0:km))
    allocate(term1(0:im,0:jm,0:km),term2(0:im,0:jm,0:km),              &
             term3(0:im,0:jm,0:km))
    !
    ! +-----------------------------------------+
    ! | start to calculate the TKE budget terms |
    ! +-----------------------------------------+
    !
    call H5_ReadGrid
    !
    allocate(u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),                    &
             u3(0:im,0:jm,0:km),p(0:im,0:jm,0:km),t(0:im,0:jm,0:km))
    !
    call H5ReadArray(u1,im,jm,km,'u1','Results/mean.fav.h5')
    call H5ReadArray(u2,im,jm,km,'u2','Results/mean.fav.h5')
    call H5ReadArray(u3,im,jm,km,'u3','Results/mean.fav.h5')
    !
    u1re=u1re-u1
    u2re=u2re-u2
    u3re=u3re-u3
    !
    allocate(p_tran(0:im,0:jm,0:km),dissip(0:im,0:jm,0:km),            &
             pres_dila(0:im,0:jm,0:km),visc_acce(0:im,0:jm,0:km),      &
             pres_acce(0:im,0:jm,0:km),visc_diff(0:im,0:jm,0:km),      &
             convec(0:im,0:jm,0:km),produc(0:im,0:jm,0:km),            &
             turb_tran(0:im,0:jm,0:km))
    !
    ! pressure transport term
    call H5ReadArray(p,im,jm,km,'p','Results/mean.fav.h5')
    call H5ReadArray(t,im,jm,km,'t','Results/mean.fav.h5')
    !
    pu1=pu1-p*(u1+u1re)
    pu2=pu2-p*(u2+u2re)
    pu3=pu3-p*(u3+u3re)
    !
    gengrad=grad(pu1,x,y,z)
    p_tran(:,:,:)=gengrad(1,:,:,:)
    !
    gengrad=grad(pu2)
    p_tran(:,:,:)=p_tran(:,:,:)+gengrad(2,:,:,:)
    !
    gengrad=grad(pu3)
    p_tran(:,:,:)=p_tran(:,:,:)+gengrad(3,:,:,:)
    !
    p_tran=-p_tran
    print*,' ** pressure transport term calculated !'
    !
    ! dissipation terms
    allocate(du1(3,0:im,0:jm,0:km),du2(3,0:im,0:jm,0:km),              &
             du3(3,0:im,0:jm,0:km),du1re(1:3,0:im,0:jm,0:km),          &
             du2re(1:3,0:im,0:jm,0:km),du3re(1:3,0:im,0:jm,0:km))
    !
    call H5ReadArray(du1(1,:,:,:),im,jm,km,'du1dx','Results/mean.fav.h5')
    call H5ReadArray(du1(2,:,:,:),im,jm,km,'du1dy','Results/mean.fav.h5')
    call H5ReadArray(du1(3,:,:,:),im,jm,km,'du1dz','Results/mean.fav.h5')
    call H5ReadArray(du2(1,:,:,:),im,jm,km,'du2dx','Results/mean.fav.h5')
    call H5ReadArray(du2(2,:,:,:),im,jm,km,'du2dy','Results/mean.fav.h5')
    call H5ReadArray(du2(3,:,:,:),im,jm,km,'du2dz','Results/mean.fav.h5')
    call H5ReadArray(du3(1,:,:,:),im,jm,km,'du3dx','Results/mean.fav.h5')
    call H5ReadArray(du3(2,:,:,:),im,jm,km,'du3dy','Results/mean.fav.h5')
    call H5ReadArray(du3(3,:,:,:),im,jm,km,'du3dz','Results/mean.fav.h5')
    !
    du1re=grad(u1re)
    du2re=grad(u2re)
    du3re=grad(u3re)
    !
    allocate(sigma11(0:im,0:jm,0:km),sigma12(0:im,0:jm,0:km),          &
             sigma13(0:im,0:jm,0:km),sigma22(0:im,0:jm,0:km),          &
             sigma23(0:im,0:jm,0:km),sigma33(0:im,0:jm,0:km),          &
             usigma1(0:im,0:jm,0:km),usigma2(0:im,0:jm,0:km),          &
             usigma3(0:im,0:jm,0:km) )
    !
    string='  ** calculating dissip, pres_dila, visc_diff terms ...'
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      miu=miucal(t(i,j,k))/Reynolds
      !
      d11=du1(1,i,j,k); d12=du1(2,i,j,k); d13=du1(3,i,j,k)
      d21=du2(1,i,j,k); d22=du2(2,i,j,k); d23=du2(3,i,j,k)
      d31=du3(1,i,j,k); d32=du3(2,i,j,k); d33=du3(3,i,j,k)
      dkk=(d11+d22+d33)/3.d0
      !
      f11=2.d0*miu*(d11-dkk)
      f22=2.d0*miu*(d22-dkk)
      f33=2.d0*miu*(d33-dkk)
      f12=miu*(d12+d21)
      f13=miu*(d13+d31)
      f23=miu*(d23+d32)
      !
      e11=du1re(1,i,j,k); e12=du1re(2,i,j,k); e13=du1re(3,i,j,k)
      e21=du2re(1,i,j,k); e22=du2re(2,i,j,k); e23=du2re(3,i,j,k)
      e31=du3re(1,i,j,k); e32=du3re(2,i,j,k); e33=du3re(3,i,j,k)
      ekk=(e11+e22+e33)/3.d0
      !
      g11=2.d0*miu*(e11-ekk)
      g22=2.d0*miu*(e22-ekk)
      g33=2.d0*miu*(e33-ekk)
      g12=miu*(e12+e21)
      g13=miu*(e13+e31)
      g23=miu*(e23+e32)
      !
      ! sigma11(i,j,k)=f11; sigma12(i,j,k)=f12; sigma13(i,j,k)=f13
      ! sigma22(i,j,k)=f22; sigma23(i,j,k)=f23
      ! sigma33(i,j,k)=f33
      !
      sigma11(i,j,k)=sgmam11m(i,j,k)
      sigma22(i,j,k)=sgmam22m(i,j,k)
      sigma33(i,j,k)=sgmam33m(i,j,k)
      sigma12(i,j,k)=sgmam12m(i,j,k)
      sigma13(i,j,k)=sgmam13m(i,j,k)
      sigma23(i,j,k)=sgmam23m(i,j,k)
      !
      ! usg1=u1(i,j,k)*f11+u1(i,j,k)*g11+u1re(i,j,k)*f11+                &
      !      u2(i,j,k)*f12+u2(i,j,k)*g12+u2re(i,j,k)*f12+                &
      !      u3(i,j,k)*f13+u3(i,j,k)*g13+u3re(i,j,k)*f13
      ! usg2=u1(i,j,k)*f12+u1(i,j,k)*g12+u1re(i,j,k)*f12+                &
      !      u2(i,j,k)*f22+u2(i,j,k)*g22+u2re(i,j,k)*f22+                &
      !      u3(i,j,k)*f23+u3(i,j,k)*g23+u3re(i,j,k)*f23
      ! usg3=u1(i,j,k)*f13+u1(i,j,k)*g13+u1re(i,j,k)*f13+                &
      !      u2(i,j,k)*f23+u2(i,j,k)*g23+u2re(i,j,k)*f23+                &
      !      u3(i,j,k)*f33+u3(i,j,k)*g33+u3re(i,j,k)*f33
      usg1=(u1(i,j,k)+u1re(i,j,k))*sgmam11m(i,j,k)+                    &
           (u2(i,j,k)+u2re(i,j,k))*sgmam12m(i,j,k)+                    &
           (u3(i,j,k)+u3re(i,j,k))*sgmam13m(i,j,k) 
      usg2=(u1(i,j,k)+u1re(i,j,k))*sgmam12m(i,j,k)+                    &
           (u2(i,j,k)+u2re(i,j,k))*sgmam22m(i,j,k)+                    &
           (u3(i,j,k)+u3re(i,j,k))*sgmam23m(i,j,k) 
      usg3=(u1(i,j,k)+u1re(i,j,k))*sgmam13m(i,j,k)+                    &
           (u2(i,j,k)+u2re(i,j,k))*sgmam23m(i,j,k)+                    &
           (u3(i,j,k)+u3re(i,j,k))*sgmam33m(i,j,k) 
      !
      usigma1(i,j,k)=visdif1(i,j,k)-usg1
      usigma2(i,j,k)=visdif2(i,j,k)-usg2
      usigma3(i,j,k)=visdif3(i,j,k)-usg3
      !
      ! dissipation term
      ! var1=f11*d11+f22*d22+f33*d33+f12*(d12+d21)+f13*(d13+d31)+f23*(d23+d32)
      ! var2=f11*e11+f22*e22+f33*e33+f12*(e12+e21)+f13*(e13+e31)+f23*(e23+e32)
      ! var3=g11*d11+g22*d22+g33*d33+g12*(d12+d21)+g13*(d13+d31)+g23*(d23+d32)
      ! !
      ! dissip(i,j,k)=disspa(i,j,k)-var1-var2-var3
      var1=sgmam11m(i,j,k)*(d11+e11)+sgmam22m(i,j,k)*(d22+e22)+        &
           sgmam33m(i,j,k)*(d33+e33)+sgmam12m(i,j,k)*(d12+d21+e12+e21)+&
           sgmam13m(i,j,k)*(d13+d31+e13+e31)+                          &
           sgmam23m(i,j,k)*(d23+d32+e23+e32) 
      dissip(i,j,k)=disspa(i,j,k)-var1 
      !
      ! if(i==800 .and. j==0 .and. k==1) then
      !   print*,disspa(i,j,k),var1
      ! endif
      !
      ! Pressure-dilatation term
      pres_dila(i,j,k)=predil(i,j,k)-p(i,j,k)*(d11+d22+d33+e11+e22+e33)
      !
    enddo
    enddo
      write(*,'(1A1,A,I4,A2,$)')char(13),trim(string),100*k/km,' %'
    enddo
    write(*,*)
    !
    print*,' ** pressure dilatation term calculated !'
    !
    ! print*,' ** dissipation at wall:',disspa(800,0,1),dissip(800,0,1)
    !
    dissip=-dissip
    print*,' ** dissipation term calculated !'
    !
    ! viscous diffusion term
    gengrad=grad(usigma1)
    visc_diff=gengrad(1,:,:,:)
    gengrad=grad(usigma2)
    visc_diff=visc_diff+gengrad(2,:,:,:)
    gengrad=grad(usigma3)
    visc_diff=visc_diff+gengrad(3,:,:,:)
    print*,' ** viscous diffusion term calculated !'
    !
    !
    ! viscous accelaration and pressure accelaration terms
    gengrad=grad(p)
    pres_acce=gengrad(1,:,:,:)*u1re+gengrad(2,:,:,:)*u2re+             &
              gengrad(3,:,:,:)*u3re
    pres_acce=-pres_acce
    print*,' ** pressure accelaration term calculated !'
    !
    ! viscous accelaration term
    write(*,'(A)',advance='no')'  ** calculating viscous accelaration term ... '
    !
    gengrad=grad(sigma11)
    visc_acce=gengrad(1,:,:,:)*u1re
    write(*,'(A)',advance='no')'16%.'
    !
    gengrad=grad(sigma12)
    visc_acce=visc_acce+gengrad(1,:,:,:)*u2re+gengrad(2,:,:,:)*u1re
    write(*,'(A)',advance='no')'33%.'
    !
    gengrad=grad(sigma13)
    visc_acce=visc_acce+gengrad(1,:,:,:)*u3re+gengrad(3,:,:,:)*u1re
    write(*,'(A)',advance='no')'50%.'
    !
    gengrad=grad(sigma22)
    visc_acce=visc_acce+gengrad(2,:,:,:)*u2re
    write(*,'(A)',advance='no')'66%.'
    !
    gengrad=grad(sigma23)
    visc_acce=visc_acce+gengrad(2,:,:,:)*u3re+gengrad(3,:,:,:)*u2re
    write(*,'(A)',advance='no')'83%.'
    !
    gengrad=grad(sigma33)
    visc_acce=visc_acce+gengrad(3,:,:,:)*u3re
    write(*,'(A)',advance='yes')'100% '
    !

    print*,' ** viscous accelaration term calculated !'
    !
    ! release memory for the production and turbulent transport terms
    deallocate(pu1,pu2,pu3,u1re,u2re,u3re,disspa,predil,visdif1,       &
               visdif2,visdif3,t,p,du1re,du2re,du3re)
    deallocate(sigma11,sigma12,sigma13,sigma22,sigma23,sigma33,        &
               usigma1,usigma2,usigma3,sgmam11m,sgmam22m,sgmam33m,     &
               sgmam12m,sgmam13m,sgmam23m )
    print*,' ** large memory released.'
    !
    allocate(u11(0:im,0:jm,0:km),u12(0:im,0:jm,0:km),                  &
             u13(0:im,0:jm,0:km),u22(0:im,0:jm,0:km),                  &
             u23(0:im,0:jm,0:km),u33(0:im,0:jm,0:km),                  &
             u111(0:im,0:jm,0:km),u112(0:im,0:jm,0:km),                &
             u113(0:im,0:jm,0:km),u122(0:im,0:jm,0:km),                &
             u222(0:im,0:jm,0:km),u223(0:im,0:jm,0:km),                &
             u133(0:im,0:jm,0:km),u233(0:im,0:jm,0:km),                &
             u333(0:im,0:jm,0:km),ro(0:im,0:jm,0:km),                  &
             tke(0:im,0:jm,0:km))
    !
    call H5ReadArray(ro,im,jm,km,'ro','Results/mean.fav.h5')
    !
    call H5ReadArray(u11,im,jm,km,'u11','Results/2order.fav.h5')
    call H5ReadArray(u22,im,jm,km,'u22','Results/2order.fav.h5')
    call H5ReadArray(u33,im,jm,km,'u33','Results/2order.fav.h5')
    call H5ReadArray(u12,im,jm,km,'u12','Results/2order.fav.h5')
    call H5ReadArray(u13,im,jm,km,'u13','Results/2order.fav.h5')
    call H5ReadArray(u23,im,jm,km,'u23','Results/2order.fav.h5')
    !
    call H5ReadArray(u111,im,jm,km,'u111','Results/3order.fav.h5')
    call H5ReadArray(u112,im,jm,km,'u112','Results/3order.fav.h5')
    call H5ReadArray(u113,im,jm,km,'u113','Results/3order.fav.h5')
    call H5ReadArray(u122,im,jm,km,'u122','Results/3order.fav.h5')
    call H5ReadArray(u222,im,jm,km,'u222','Results/3order.fav.h5')
    call H5ReadArray(u223,im,jm,km,'u223','Results/3order.fav.h5')
    call H5ReadArray(u133,im,jm,km,'u133','Results/3order.fav.h5')
    call H5ReadArray(u233,im,jm,km,'u233','Results/3order.fav.h5')
    call H5ReadArray(u333,im,jm,km,'u333','Results/3order.fav.h5')
    !
    ! convection term
    tke=0.5d0*(u11+u22+u33)
    term1=ro*u1*tke
    term2=ro*u2*tke
    term3=ro*u3*tke
    gengrad=grad(term1)
    convec=gengrad(1,:,:,:)
    gengrad=grad(term2)
    convec=convec+gengrad(2,:,:,:)
    gengrad=grad(term3)
    convec=convec+gengrad(3,:,:,:)
    !
    convec=-convec
    print*,' ** convection term calculated !'
    !
    ! turbulent transport terms
    term1=0.5d0*ro*(u111+u122+u133)
    term2=0.5d0*ro*(u112+u222+u233)
    term3=0.5d0*ro*(u113+u223+u333)
    gengrad=grad(term1)
    turb_tran=gengrad(1,:,:,:)
    gengrad=grad(term2)
    turb_tran=turb_tran+gengrad(2,:,:,:)
    gengrad=grad(term3)
    turb_tran=turb_tran+gengrad(3,:,:,:)
    !
    turb_tran=-turb_tran
    print*,' ** turbulent transport term calculated !'
    !
    ! production term
    produc=ro*(u11*du1(1,:,:,:)+u22*du2(2,:,:,:)+u33*du3(3,:,:,:)+     &
               u12*(du1(2,:,:,:)+du2(1,:,:,:))+                        &
               u13*(du1(3,:,:,:)+du3(1,:,:,:))+                        &
               u23*(du2(3,:,:,:)+du3(2,:,:,:)) )
    !
    produc=-produc
    print*,' ** production term calculated !'
    !
    deallocate(du1,du2,du3,term1,term2,term3,gengrad,ro,u1,u2,u3,      &
               u11,u12,u13,u22,u23,u33,u111,u112,u113,u122,u222,u223,  &
               u133,u233,u333)
    print*,' ** large memory released.'
    !
    ! write 3-D budget data
    inquire(file='Results/tke_budget.h5',exist=filexists)
    if(filexists) call system('mv -v Results/tke_budget.h5 Results/tke_budget.bak')
    !
    call H5WriteArray(p_tran,im,jm,km,'pressure_transport',            &
                                                'Results/tke_budget.h5') 
    call H5WriteArray(pres_dila,im,jm,km,'pressure_dilatation',        &
                                                'Results/tke_budget.h5')
    call H5WriteArray(dissip,im,jm,km,'dissipation',                   &
                                                'Results/tke_budget.h5')
    call H5WriteArray(produc,im,jm,km,'production',                    &
                                                'Results/tke_budget.h5')
    call H5WriteArray(turb_tran,im,jm,km,'turbulent_transport',        &
                                                'Results/tke_budget.h5')
    call H5WriteArray(convec,im,jm,km,'convection',                    &
                                                'Results/tke_budget.h5')
    call H5WriteArray(visc_acce,im,jm,km,'viscous_accelaration',       &
                                                'Results/tke_budget.h5')
    call H5WriteArray(pres_acce,im,jm,km,'pressure_accelaration',      &
                                                'Results/tke_budget.h5')
    call H5WriteArray(visc_diff,im,jm,km,'viscous_diffusion',          &
                                                'Results/tke_budget.h5')
    !
    if(lkhomo) then
      ! 
      print*,' ** averaging in the spanwise direction'
      !
      allocate(p_tran_zm(0:im,0:jm),dissip_zm(0:im,0:jm),              &
               pres_dila_zm(0:im,0:jm),visc_acce_zm(0:im,0:jm),        &
               pres_acce_zm(0:im,0:jm),visc_diff_zm(0:im,0:jm),        &
               convec_zm(0:im,0:jm),produc_zm(0:im,0:jm),              &
               turb_tran_zm(0:im,0:jm))
      !
      do j=0,jm
      do i=0,im
        !
        p_tran_zm(i,j)=integration(p_tran(i,j,:))/dble(km)
        dissip_zm(i,j)=integration(dissip(i,j,:))/dble(km)
        convec_zm(i,j)=integration(convec(i,j,:))/dble(km)
        produc_zm(i,j)=integration(produc(i,j,:))/dble(km)
        !
        pres_dila_zm(i,j)=integration(pres_dila(i,j,:))/dble(km)
        visc_acce_zm(i,j)=integration(visc_acce(i,j,:))/dble(km)
        pres_acce_zm(i,j)=integration(pres_acce(i,j,:))/dble(km)
        visc_diff_zm(i,j)=integration(visc_diff(i,j,:))/dble(km)
        turb_tran_zm(i,j)=integration(turb_tran(i,j,:))/dble(km)
        !
      enddo
      enddo
      !
      inquire(file='Results/tke_budget.zm.h5',exist=filexists)
      if(filexists) call system('mv -v Results/tke_budget.zm.h5 Results/tke_budget.zm.bak')
      !
      call H5WriteArray(convec_zm,im,jm,'convection',                  &
                                             'Results/tke_budget.zm.h5')
      call H5WriteArray(turb_tran_zm,im,jm,'turbulent_transport',      &
                                             'Results/tke_budget.zm.h5')
      call H5WriteArray(p_tran_zm,im,jm,'pressure_transport',          &
                                             'Results/tke_budget.zm.h5')
      call H5WriteArray(visc_diff_zm,im,jm,'viscous_diffusion',        &
                                             'Results/tke_budget.zm.h5')
      call H5WriteArray(produc_zm,im,jm,'production',                  &
                                             'Results/tke_budget.zm.h5')
      call H5WriteArray(pres_dila_zm,im,jm,'pressure_dilatation',      &
                                             'Results/tke_budget.zm.h5')
      call H5WriteArray(pres_acce_zm,im,jm,'pressure_accelaration',    &
                                             'Results/tke_budget.zm.h5')
      call H5WriteArray(visc_acce_zm,im,jm,'viscous_accelaration',     &
                                             'Results/tke_budget.zm.h5')
      call H5WriteArray(dissip_zm,im,jm,'dissipation',                 &
                                             'Results/tke_budget.zm.h5')
      !
      call writetecbin('Results/tecbudget_zm.plt',                     &
                               x(:,:,0),'x',y(:,:,0),'y',              &
                                  convec_zm,'convection',              &
                              turb_tran_zm,'turbulent transport',      &
                              p_tran_zm,'pressure transport',          &
                              visc_diff_zm,'viscous diffusion',        &
                              produc_zm,'production',                  &
                              pres_dila_zm,'pressure dilatation',      &
                              pres_acce_zm,'pressure accelaration',    &
                              visc_acce_zm,'viscous accelaration',     &
                              dissip_zm,'dissipation',im,jm)
    endif
    !
  end subroutine TKEBudgetSta
  !!
  subroutine TKEBudget_zm
    !
    use commvardefine, only: im,jm,km,lihomo,ljhomo,lkhomo,Reynolds,gridfile
    use h5readwrite
    use basicfunction, only: integration,miucal
    use gradsolver
    use writetec
    !
    integer :: nsamples,i,j
    !
    real(8),allocatable,dimension(:,:) :: x,y,u1,u2,p,t,ro,            &
                                          pu1,pu2,u1re,u2re,           &
                                          disspa,predil,visdif1,       &
                                          visdif2,                     &
                                          sigma11,sigma12,sigma22,     &
                                          usigma1,usigma2,             &
                                          u11,u12,u22,u33,tke,         &
                                          u111,u112,u122,u222,u133,u233
    real(8),allocatable,dimension(:,:,:) :: du1re,du2re,du1,du2,gengrad
    ! TKE budget terms
    real(8),allocatable,dimension(:,:) :: p_tran,dissip,pres_dila,   &
                                          visc_acce,pres_acce,       &
                                          visc_diff,convec,produc,   &
                                          turb_tran
    real(8),allocatable,dimension(:) :: convec_xzm,turb_tran_xzm,    &
                             p_tran_xzm,visc_diff_xzm,produc_xzm,    &
                             pres_dila_xzm,pres_acce_xzm,visc_acce_xzm, &
                             dissip_xzm
    !
    real(8),allocatable,dimension(:,:) :: term1,term2
    !
    real(8) :: rsamples,miu,var1,var2,var3
    real(8) :: f11,f12,f22,d11,d12,d21,d22,dkk,g11,g12,g22,e11,e12,    &
               e21,e22,ekk,usg1,usg2
    logical :: filexists
    character(len=255) :: string
    !
    !+--------------------------+
    !| programme starts         |
    !+--------------------------+
    !
    allocate(pu1(0:im,0:jm),pu2(0:im,0:jm),u1re(0:im,0:jm),            &
             u2re(0:im,0:jm),disspa(0:im,0:jm),predil(0:im,0:jm)) 
    allocate(visdif1(0:im,0:jm),visdif2(0:im,0:jm))
    !
    allocate(sigma11(0:im,0:jm),sigma12(0:im,0:jm),sigma22(0:im,0:jm), &
             usigma1(0:im,0:jm),usigma2(0:im,0:jm) )
    !
    call H5ReadArray(disspa,im,jm,'disspa','Results/budget.zm.h5')
    call H5ReadArray(predil,im,jm,'predil','Results/budget.zm.h5')
    call H5ReadArray(pu1,im,jm,'pu1','Results/budget.zm.h5')
    call H5ReadArray(pu2,im,jm,'pu2','Results/budget.zm.h5')
    call H5ReadArray(u1re,im,jm,'u1rem','Results/budget.zm.h5')
    call H5ReadArray(u2re,im,jm,'u2rem','Results/budget.zm.h5')
    call H5ReadArray(visdif1,im,jm,'visdif1','Results/budget.zm.h5')
    call H5ReadArray(visdif2,im,jm,'visdif2','Results/budget.zm.h5')
    !
    call H5ReadArray(sigma11,im,jm,'sgmam11','Results/budget.zm.h5')
    call H5ReadArray(sigma22,im,jm,'sgmam22','Results/budget.zm.h5')
    call H5ReadArray(sigma12,im,jm,'sgmam12','Results/budget.zm.h5')
    !
    allocate(x(0:im,0:jm),y(0:im,0:jm),u1(0:im,0:jm),u2(0:im,0:jm),    &
             p(0:im,0:jm),t(0:im,0:jm))
    !
    call h5_read2dfrom3d(x,im,jm,km,'x',trim(gridfile),kslice=0)
    call h5_read2dfrom3d(y,im,jm,km,'y',trim(gridfile),kslice=0)
    call H5ReadArray(u1,im,jm,'u1','Results/mean.fav.zm.h5')
    call H5ReadArray(u2,im,jm,'u2','Results/mean.fav.zm.h5')
    call H5ReadArray(p,im,jm,'p','Results/mean.fav.zm.h5')
    call H5ReadArray(t,im,jm,'t','Results/mean.fav.zm.h5')
    !
    u1re=u1re-u1
    u2re=u2re-u2
    !
    allocate(gengrad(1:2,0:im,0:jm),term1(0:im,0:jm),term2(0:im,0:jm))
    !
    ! +-----------------------------------------+
    ! | start to calculate the TKE budget terms |
    ! +-----------------------------------------+
    !
    allocate(p_tran(0:im,0:jm),dissip(0:im,0:jm),pres_dila(0:im,0:jm), &
             visc_acce(0:im,0:jm),pres_acce(0:im,0:jm),                &
             visc_diff(0:im,0:jm),convec(0:im,0:jm),produc(0:im,0:jm), &
             turb_tran(0:im,0:jm))
    !
    ! pressure transport term
    pu1=pu1-p*(u1+u1re)
    pu2=pu2-p*(u2+u2re)
    !
    gengrad=grad(pu1,x,y)
    p_tran=gengrad(1,:,:)
    !
    gengrad=grad(pu2)
    p_tran=p_tran+gengrad(2,:,:)
    !
    p_tran=-p_tran
    print*,' ** pressure transport term calculated !'
    !
    !
    ! dissipation terms
    allocate(du1(2,0:im,0:jm),du2(2,0:im,0:jm),                        &
             du1re(2,0:im,0:jm),du2re(2,0:im,0:jm))
    ! !
    call H5ReadArray(du1(1,:,:),im,jm,'du1dx','Results/mean.fav.zm.h5')
    call H5ReadArray(du1(2,:,:),im,jm,'du1dy','Results/mean.fav.zm.h5')
    call H5ReadArray(du2(1,:,:),im,jm,'du2dx','Results/mean.fav.zm.h5')
    call H5ReadArray(du2(2,:,:),im,jm,'du2dy','Results/mean.fav.zm.h5')
    !
    du1re=grad(u1re)
    du2re=grad(u2re)
    ! !    !
    string='  ** calculating dissipation, pressure-dilatation, viscous diffusion terms ...'
    ! !
    do j=0,jm
    do i=0,im
      !
      miu=miucal(t(i,j))/Reynolds
      !
      d11=du1(1,i,j); d12=du1(2,i,j)
      d21=du2(1,i,j); d22=du2(2,i,j)
      !
      dkk=(d11+d22)/3.d0
      !
      f11=sigma11(i,j)
      f22=sigma22(i,j)
      f12=sigma12(i,j)
      !
      e11=du1re(1,i,j); e12=du1re(2,i,j)
      e21=du2re(1,i,j); e22=du2re(2,i,j)
      ekk=(e11+e22)/3.d0
      !
      g11=2.d0*miu*(e11-ekk)
      g22=2.d0*miu*(e22-ekk)
      g12=miu*(e12+e21)
      !
      usg1=(u1(i,j)+u1re(i,j))*f11+(u2(i,j)+u2re(i,j))*f12
      usg2=(u1(i,j)+u1re(i,j))*f12+(u2(i,j)+u2re(i,j))*f22
      !
      usigma1(i,j)=visdif1(i,j)-usg1
      usigma2(i,j)=visdif2(i,j)-usg2
      !
      ! dissipation term
      var1=f11*(d11+e11)+f22*(d22+e22)+f12*(d12+d21+e12+e21)
      !
      dissip(i,j)=disspa(i,j)-var1
      !
      ! Pressure-dilatation term
      pres_dila(i,j)=predil(i,j)-p(i,j)*(d11+d22+e11+e22)
      !
    enddo
      write(*,'(1A1,A,I4,A2,$)')char(13),trim(string),100*j/jm,' %'
    enddo
    write(*,*)
    !
    print*,' ** pressure dilatation term calculated !'
    !
    dissip=-dissip
    print*,' ** dissipation term calculated !'
    !
    ! viscous diffusion term
    gengrad=grad(usigma1)
    visc_diff=gengrad(1,:,:)
    gengrad=grad(usigma2)
    visc_diff=visc_diff+gengrad(2,:,:)
    print*,' ** viscous diffusion term calculated !'
    ! !
    ! !
    ! viscous accelaration and pressure accelaration terms
    gengrad=grad(p)
    pres_acce=gengrad(1,:,:)*u1re+gengrad(2,:,:)*u2re
    pres_acce=-pres_acce
    print*,' ** pressure accelaration term calculated !'
    !
    ! viscous accelaration term
    write(*,'(A)',advance='no')'  ** calculating viscous accelaration term ... '
    !
    gengrad=grad(sigma11)
    visc_acce=gengrad(1,:,:)*u1re
    write(*,'(A)',advance='no')'33%.'
    !
    gengrad=grad(sigma12)
    visc_acce=visc_acce+gengrad(1,:,:)*u2re+gengrad(2,:,:)*u1re
    write(*,'(A)',advance='no')'67%.'
    !
    gengrad=grad(sigma22)
    visc_acce=visc_acce+gengrad(2,:,:)*u2re
    write(*,'(A)',advance='yes')'100%.'
    !
    print*,' ** viscous accelaration term calculated !'
    ! !
    allocate(u11(0:im,0:jm),u12(0:im,0:jm),u22(0:im,0:jm),             &
             u33(0:im,0:jm),u111(0:im,0:jm),u112(0:im,0:jm),           &
             u122(0:im,0:jm),u222(0:im,0:jm),u133(0:im,0:jm),          &
             u233(0:im,0:jm),tke(0:im,0:jm),ro(0:im,0:jm))
    !
    call H5ReadArray(ro,im,jm,'ro','Results/mean.fav.zm.h5')
    !
    call H5ReadArray(u11,im,jm,'u11','Results/2order.fav.zm.h5')
    call H5ReadArray(u22,im,jm,'u22','Results/2order.fav.zm.h5')
    call H5ReadArray(u33,im,jm,'u33','Results/2order.fav.zm.h5')
    call H5ReadArray(u12,im,jm,'u12','Results/2order.fav.zm.h5')
    !
    call H5ReadArray(u111,im,jm,'u111','Results/3order.fav.zm.h5')
    call H5ReadArray(u112,im,jm,'u112','Results/3order.fav.zm.h5')
    call H5ReadArray(u122,im,jm,'u122','Results/3order.fav.zm.h5')
    call H5ReadArray(u222,im,jm,'u222','Results/3order.fav.zm.h5')
    call H5ReadArray(u133,im,jm,'u133','Results/3order.fav.zm.h5')
    call H5ReadArray(u233,im,jm,'u233','Results/3order.fav.zm.h5')
    ! !
    ! convection term
    tke=0.5d0*(u11+u22+u33)
    term1=ro*u1*tke
    term2=ro*u2*tke
    gengrad=grad(term1)
    convec=gengrad(1,:,:)
    gengrad=grad(term2)
    convec=convec+gengrad(2,:,:)
    !
    convec=-convec
    print*,' ** convection term calculated !'
    ! !
    ! turbulent transport terms
    term1=0.5d0*ro*(u111+u122+u133)
    term2=0.5d0*ro*(u112+u222+u233)
    gengrad=grad(term1)
    turb_tran=gengrad(1,:,:)
    gengrad=grad(term2)
    turb_tran=turb_tran+gengrad(2,:,:)
    !
    turb_tran=-turb_tran
    print*,' ** turbulent transport term calculated !'
    ! !
    ! production term
    produc=ro*(u11*du1(1,:,:)+u22*du2(2,:,:)+u12*(du1(2,:,:)+du2(1,:,:)))
    !
    produc=-produc
    print*,' ** production term calculated !'
    !
    inquire(file='Results/tke_budget.zm.h5',exist=filexists)
    !
    if(filexists) call system('mv -v Results/tke_budget.zm.h5 Results/tke_budget.zm.bak')
    !
    call H5WriteArray(convec,im,jm,'convection',                       &
                                             'Results/tke_budget.zm.h5')
    call H5WriteArray(turb_tran,im,jm,'turbulent_transport',           &
                                             'Results/tke_budget.zm.h5')
    call H5WriteArray(p_tran,im,jm,'pressure_transport',               &
                                             'Results/tke_budget.zm.h5')
    call H5WriteArray(visc_diff,im,jm,'viscous_diffusion',             &
                                             'Results/tke_budget.zm.h5')
    call H5WriteArray(produc,im,jm,'production',                       &
                                             'Results/tke_budget.zm.h5')
    call H5WriteArray(pres_dila,im,jm,'pressure_dilatation',           &
                                             'Results/tke_budget.zm.h5')
    call H5WriteArray(pres_acce,im,jm,'pressure_accelaration',         &
                                             'Results/tke_budget.zm.h5')
    call H5WriteArray(visc_acce,im,jm,'viscous_accelaration',          &
                                             'Results/tke_budget.zm.h5')
    call H5WriteArray(dissip,im,jm,'dissipation',                      &
                                             'Results/tke_budget.zm.h5')
    !
    call writetecbin('Results/tecbudget_zm.plt', x,'x',y,'y',          &
                                     convec,'convection',              &
                                 turb_tran,'turbulent transport',      &
                                 p_tran,'pressure transport',          &
                                 visc_diff,'viscous diffusion',        &
                                 produc,'production',                  &
                                 pres_dila,'pressure dilatation',      &
                                 pres_acce,'pressure accelaration',    &
                                 visc_acce,'viscous accelaration',     &
                                   dissip,'dissipation',im,jm)
    !
    if(lihomo .and. lkhomo) then
      !
      allocate(convec_xzm(0:jm),turb_tran_xzm(0:jm),p_tran_xzm(0:jm),   &
               visc_diff_xzm(0:jm),produc_xzm(0:jm),pres_dila_xzm(0:jm),&
               pres_acce_xzm(0:jm),visc_acce_xzm(0:jm),dissip_xzm(0:jm))
      !
      do j=0,jm
        convec_xzm(j)   =integration(   convec(:,j))/dble(im)
        turb_tran_xzm(j)=integration(turb_tran(:,j))/dble(im)
        p_tran_xzm(j)   =integration(   p_tran(:,j))/dble(im)
        visc_diff_xzm(j)=integration(visc_diff(:,j))/dble(im)
        produc_xzm(j)   =integration(   produc(:,j))/dble(im)
        pres_dila_xzm(j)=integration(pres_dila(:,j))/dble(im)
        pres_acce_xzm(j)=integration(pres_acce(:,j))/dble(im)
        visc_acce_xzm(j)=integration(visc_acce(:,j))/dble(im)
        dissip_xzm(j)   =integration(   dissip(:,j))/dble(im)
      enddo
      !
      if(filexists) call system('mv -v Results/tke_budget.xzm.h5 Results/tke_budget.xzm.bak')
      !
      call H5WriteArray(convec_xzm,jm,'convection','Results/tke_budget.xzm.h5')
      call H5WriteArray(turb_tran_xzm,jm,'turbulent_transport','Results/tke_budget.xzm.h5')
      call H5WriteArray(p_tran_xzm,jm,'pressure_transport','Results/tke_budget.xzm.h5')
      call H5WriteArray(visc_diff_xzm,jm,'viscous_diffusion','Results/tke_budget.xzm.h5')
      call H5WriteArray(produc_xzm,jm,'production','Results/tke_budget.xzm.h5')
      call H5WriteArray(pres_dila_xzm,jm,'pressure_dilatation','Results/tke_budget.xzm.h5')
      call H5WriteArray(pres_acce_xzm,jm,'pressure_accelaration','Results/tke_budget.xzm.h5')
      call H5WriteArray(visc_acce_xzm,jm,'viscous_accelaration','Results/tke_budget.xzm.h5')
      call H5WriteArray(dissip_xzm,jm,'dissipation','Results/tke_budget.xzm.h5')
      !
    endif
    !
  end subroutine TKEBudget_zm
  !
  !+-------------------------------------------------------------------+
  !| The subroutine is used to process 2nd and 3rd-order statistics.   |
  !+-------------------------------------------------------------------+
  !| Writen by Jian Fang, 2020-04-18                                   |
  !+-------------------------------------------------------------------+
  subroutine highordermom
    !
    use commvardefine, only: im,jm,km,ro_m,u1_m,u2_m,u3_m,p_m,t_m,     &
                             ro_xzm,u1_xzm,u2_xzm,u3_xzm,p_xzm,t_xzm,  &
                             ro_zm,u1_zm,u2_zm,u3_zm,p_zm,t_zm,        &
                             lihomo,ljhomo,lkhomo,x,y
    use h5readwrite
    use writetec
    use basicfunction, only: integration
    !
    real(8),allocatable,dimension(:,:,:) :: u11,u12,u13,u22,u23,u33,   &
                                            u111,u112,u113,u122,u123,  &
                                            u133,u222,u223,u233,u333,  &
                                            tt,pp,tu1,tu2,tu3
    real(8),allocatable,dimension(:) :: u11_xzm,u12_xzm,u13_xzm,       &
               u22_xzm,u23_xzm,u33_xzm,u111_xzm,u112_xzm,u113_xzm,     &
               u122_xzm,u123_xzm,u133_xzm,u222_xzm,u223_xzm,u233_xzm,  &
               u333_xzm,tt_xzm,pp_xzm,tu1_xzm,tu2_xzm,tu3_xzm
    real(8),allocatable,dimension(:,:) :: u11_zm,u12_zm,u13_zm,u22_zm, &
                u23_zm,u33_zm,u111_zm,u112_zm,u113_zm,u122_zm,u123_zm, &
           u133_zm,u222_zm,u223_zm,u233_zm,u333_zm,tt_zm,pp_zm,tu1_zm, &
          tu2_zm,tu3_zm
    !
    integer :: nsamples,i,j,k
    real(8) :: rsamples
    logical :: filexists=.false.
    character(len=64) :: file2order,file3order
    !
    allocate(u11(0:im,0:jm,0:km),u12(0:im,0:jm,0:km),                  &
             u22(0:im,0:jm,0:km),u13(0:im,0:jm,0:km),                  &
             u23(0:im,0:jm,0:km),u33(0:im,0:jm,0:km),                  &
             u111(0:im,0:jm,0:km),u112(0:im,0:jm,0:km),                &
             u113(0:im,0:jm,0:km),u122(0:im,0:jm,0:km),                &
             u123(0:im,0:jm,0:km),u133(0:im,0:jm,0:km),                &
             u222(0:im,0:jm,0:km),u223(0:im,0:jm,0:km),                &
             u233(0:im,0:jm,0:km),u333(0:im,0:jm,0:km),                &
             tt(0:im,0:jm,0:km),pp(0:im,0:jm,0:km),                    &
             tu1(0:im,0:jm,0:km),tu2(0:im,0:jm,0:km),                  &
             tu3(0:im,0:jm,0:km))
    !
    ! file2order='meanflow.h5'
    ! file3order='meanflow.h5'
    
    file2order='2ndsta.h5'
    file3order='3rdsta.h5'
    call H5ReadArray(nsamples,'nsamples',outfolder//trim(file2order))
    !
    write(*,'(A,I5)')'  ** Total number of temporal samples:',nsamples
    !
    call H5ReadArray(u11,im,jm,km,'u11',outfolder//trim(file2order))
    call H5ReadArray(u12,im,jm,km,'u12',outfolder//trim(file2order))
    call H5ReadArray(u13,im,jm,km,'u13',outfolder//trim(file2order))
    call H5ReadArray(u22,im,jm,km,'u22',outfolder//trim(file2order))
    call H5ReadArray(u23,im,jm,km,'u23',outfolder//trim(file2order))
    call H5ReadArray(u33,im,jm,km,'u33',outfolder//trim(file2order))
    !
    call H5ReadArray( tt,im,jm,km, 'tt',outfolder//trim(file2order))
    call H5ReadArray( pp,im,jm,km, 'pp',outfolder//trim(file2order))
    call H5ReadArray(tu1,im,jm,km,'tu1',outfolder//trim(file2order))
    call H5ReadArray(tu2,im,jm,km,'tu2',outfolder//trim(file2order))
    call H5ReadArray(tu3,im,jm,km,'tu3',outfolder//trim(file2order))
    !
    call H5_ReadArray3d(u111,im,jm,km,'u111',outfolder//trim(file3order))
    call H5_ReadArray3d(u222,im,jm,km,'u222',outfolder//trim(file3order))
    call H5_ReadArray3d(u333,im,jm,km,'u333',outfolder//trim(file3order))
    call H5_ReadArray3d(u112,im,jm,km,'u112',outfolder//trim(file3order))
    call H5_ReadArray3d(u113,im,jm,km,'u113',outfolder//trim(file3order))
    call H5_ReadArray3d(u122,im,jm,km,'u122',outfolder//trim(file3order))
    call H5_ReadArray3d(u133,im,jm,km,'u133',outfolder//trim(file3order))
    call H5_ReadArray3d(u223,im,jm,km,'u223',outfolder//trim(file3order))
    call H5_ReadArray3d(u233,im,jm,km,'u233',outfolder//trim(file3order))
    call H5_ReadArray3d(u123,im,jm,km,'u123',outfolder//trim(file3order))
    !
    if(lihomo .and. lkhomo) then
      !
      allocate(u11_xzm(0:jm),u12_xzm(0:jm),                    &
               u22_xzm(0:jm),u13_xzm(0:jm),                    &
               u23_xzm(0:jm),u33_xzm(0:jm),                    &
               u111_xzm(0:jm),u112_xzm(0:jm),                  &
               u113_xzm(0:jm),u122_xzm(0:jm),                  &
               u123_xzm(0:jm),u133_xzm(0:jm),                  &
               u222_xzm(0:jm),u223_xzm(0:jm),                  &
               u233_xzm(0:jm),u333_xzm(0:jm),                  &
                tt_xzm(0:jm),pp_xzm(0:jm),                     &
               tu1_xzm(0:jm),tu2_xzm(0:jm),                    &
               tu3_xzm(0:jm))
      !
      rsamples=real(nsamples*im*km,8)
      !
      print*,' ** averaging along i & k direction ... '
      !
      do j=0,jm
        !
        u11_xzm(j) = integration(u11(:,j,:))/(ro_xzm(j)*rsamples)
        u22_xzm(j) = integration(u22(:,j,:))/(ro_xzm(j)*rsamples)
        u33_xzm(j) = integration(u33(:,j,:))/(ro_xzm(j)*rsamples)
        u12_xzm(j) = integration(u12(:,j,:))/(ro_xzm(j)*rsamples)
        u13_xzm(j) = integration(u13(:,j,:))/(ro_xzm(j)*rsamples)
        u23_xzm(j) = integration(u23(:,j,:))/(ro_xzm(j)*rsamples)
        ! 
        u111_xzm(j)=integration(u111(:,j,:))/(ro_xzm(j)*rsamples)
        u222_xzm(j)=integration(u222(:,j,:))/(ro_xzm(j)*rsamples)
        u333_xzm(j)=integration(u333(:,j,:))/(ro_xzm(j)*rsamples)
        u112_xzm(j)=integration(u112(:,j,:))/(ro_xzm(j)*rsamples)
        u113_xzm(j)=integration(u113(:,j,:))/(ro_xzm(j)*rsamples)
        u122_xzm(j)=integration(u122(:,j,:))/(ro_xzm(j)*rsamples)
        u133_xzm(j)=integration(u133(:,j,:))/(ro_xzm(j)*rsamples)
        u223_xzm(j)=integration(u223(:,j,:))/(ro_xzm(j)*rsamples)
        u233_xzm(j)=integration(u233(:,j,:))/(ro_xzm(j)*rsamples)
        u123_xzm(j)=integration(u123(:,j,:))/(ro_xzm(j)*rsamples)
        
         tt_xzm(j) = integration(tt(:,j,:))/(ro_xzm(j)*rsamples)
        tu1_xzm(j) = integration(tu1(:,j,:))/(ro_xzm(j)*rsamples)
        tu2_xzm(j) = integration(tu2(:,j,:))/(ro_xzm(j)*rsamples)
        tu3_xzm(j) = integration(tu3(:,j,:))/(ro_xzm(j)*rsamples)
        !
        pp_xzm(j) = integration(pp(:,j,:))/rsamples
        !
      end do
      !
      print*,' ** averaging along i & k direction ... done.'
      pp_xzm= pp_xzm-p_xzm**2
      !
      u11_xzm=u11_xzm-u1_xzm**2
      u22_xzm=u22_xzm-u2_xzm**2
      u33_xzm=u33_xzm-u3_xzm**2
      u12_xzm=u12_xzm-u1_xzm*u2_xzm
      u13_xzm=u13_xzm-u1_xzm*u3_xzm
      u23_xzm=u23_xzm-u2_xzm*u3_xzm
       tt_xzm= tt_xzm-t_xzm**2
      tu1_xzm=tu1_xzm-u1_xzm*t_xzm
      tu2_xzm=tu2_xzm-u2_xzm*t_xzm
      tu3_xzm=tu3_xzm-u3_xzm*t_xzm
      !
      u111_xzm=u111_xzm-u1_xzm**3-3.d0*u11_xzm*u1_xzm
      u222_xzm=u222_xzm-u2_xzm**3-3.d0*u22_xzm*u2_xzm
      u333_xzm=u333_xzm-u3_xzm**3-3.d0*u33_xzm*u3_xzm
      u112_xzm=u112_xzm-u1_xzm**2*u2_xzm-u11_xzm*u2_xzm-2.d0*u12_xzm*u1_xzm
      u113_xzm=u113_xzm-u1_xzm**2*u3_xzm-u11_xzm*u3_xzm-2.d0*u13_xzm*u1_xzm
      u122_xzm=u122_xzm-u1_xzm*u2_xzm**2-u22_xzm*u1_xzm-2.d0*u12_xzm*u2_xzm
      u133_xzm=u133_xzm-u1_xzm*u3_xzm**2-u33_xzm*u1_xzm-2.d0*u13_xzm*u3_xzm
      u223_xzm=u223_xzm-u3_xzm*u2_xzm**2-u22_xzm*u3_xzm-2.d0*u23_xzm*u2_xzm
      u233_xzm=u233_xzm-u2_xzm*u3_xzm**2-u33_xzm*u2_xzm-2.d0*u23_xzm*u3_xzm
      u123_xzm=u123_xzm-u1_xzm*u2_xzm*u3_xzm-u12_xzm*u3_xzm-u13_xzm*u2_xzm-u23_xzm*u1_xzm
      !
      inquire(file='Results/2order.fav.xzm.h5',exist=filexists)
      !
      if(filexists) then
        call system('mv -v Results/2order.fav.xzm.h5 Results/2order.fav.xzm.bak')
      endif
      !
      call H5WriteArray(u11_xzm,jm,'u11','Results/2order.fav.xzm.h5')
      call H5WriteArray(u22_xzm,jm,'u22','Results/2order.fav.xzm.h5')
      call H5WriteArray(u33_xzm,jm,'u33','Results/2order.fav.xzm.h5')
      call H5WriteArray(u12_xzm,jm,'u12','Results/2order.fav.xzm.h5')
      call H5WriteArray(u13_xzm,jm,'u13','Results/2order.fav.xzm.h5')
      call H5WriteArray(u23_xzm,jm,'u23','Results/2order.fav.xzm.h5')
      call H5WriteArray( tt_xzm,jm, 'tt','Results/2order.fav.xzm.h5')
      call H5WriteArray(tu1_xzm,jm,'tu1','Results/2order.fav.xzm.h5')
      call H5WriteArray(tu2_xzm,jm,'tu2','Results/2order.fav.xzm.h5')
      call H5WriteArray(tu3_xzm,jm,'tu3','Results/2order.fav.xzm.h5')
      call H5WriteArray( pp_xzm,jm, 'pp','Results/2order.fav.xzm.h5')
      
      inquire(file='Results/3order.fav.xzm.h5',exist=filexists)
      !
      if(filexists) then
        call system('mv -v Results/3order.fav.xzm.h5 Results/3order.fav.xzm.bak')
      endif
      !
      call H5WriteArray(u111_xzm,jm,'u111','Results/3order.fav.xzm.h5')
      call H5WriteArray(u222_xzm,jm,'u222','Results/3order.fav.xzm.h5')
      call H5WriteArray(u333_xzm,jm,'u333','Results/3order.fav.xzm.h5')
      call H5WriteArray(u112_xzm,jm,'u112','Results/3order.fav.xzm.h5')
      call H5WriteArray(u113_xzm,jm,'u113','Results/3order.fav.xzm.h5')
      call H5WriteArray(u122_xzm,jm,'u122','Results/3order.fav.xzm.h5')
      call H5WriteArray(u133_xzm,jm,'u133','Results/3order.fav.xzm.h5')
      call H5WriteArray(u223_xzm,jm,'u223','Results/3order.fav.xzm.h5')
      call H5WriteArray(u233_xzm,jm,'u233','Results/3order.fav.xzm.h5')
      call H5WriteArray(u123_xzm,jm,'u123','Results/3order.fav.xzm.h5')
      !
      deallocate(u11_xzm,u12_xzm,u13_xzm,u22_xzm,u23_xzm,u33_xzm,u111_xzm,    &
                 u112_xzm,u113_xzm,u122_xzm,u123_xzm,u133_xzm,u222_xzm,      &
                 u223_xzm,u233_xzm,u333_xzm,tt_xzm,pp_xzm,tu1_xzm,tu2_xzm,tu3_xzm)
      !
    endif
    !
    if(lkhomo) then
      !
      allocate(u11_zm(0:im,0:jm),u12_zm(0:im,0:jm),                    &
               u22_zm(0:im,0:jm),u13_zm(0:im,0:jm),                    &
               u23_zm(0:im,0:jm),u33_zm(0:im,0:jm),                    &
               u111_zm(0:im,0:jm),u112_zm(0:im,0:jm),                  &
               u113_zm(0:im,0:jm),u122_zm(0:im,0:jm),                  &
               u123_zm(0:im,0:jm),u133_zm(0:im,0:jm),                  &
               u222_zm(0:im,0:jm),u223_zm(0:im,0:jm),                  &
               u233_zm(0:im,0:jm),u333_zm(0:im,0:jm),                  &
                tt_zm(0:im,0:jm),pp_zm(0:im,0:jm),                     &
               tu1_zm(0:im,0:jm),tu2_zm(0:im,0:jm),                    &
               tu3_zm(0:im,0:jm))
      !
      if(km==0) then
        rsamples=dble(nsamples)
      else
        rsamples=real(nsamples*km,8)
      endif
      !
      print*,' ** averaging along k direction ... '
      print*,allocated(ro_zm),allocated(u11),allocated(u11_zm)
      !
      !
      do j=0,jm
      do i=0,im
        !
        u11_zm(i,j) = integration(u11(i,j,:))/(ro_zm(i,j)*rsamples)
        u22_zm(i,j) = integration(u22(i,j,:))/(ro_zm(i,j)*rsamples)
        u33_zm(i,j) = integration(u33(i,j,:))/(ro_zm(i,j)*rsamples)
        u12_zm(i,j) = integration(u12(i,j,:))/(ro_zm(i,j)*rsamples)
        u13_zm(i,j) = integration(u13(i,j,:))/(ro_zm(i,j)*rsamples)
        u23_zm(i,j) = integration(u23(i,j,:))/(ro_zm(i,j)*rsamples)
        ! 
        u111_zm(i,j)=integration(u111(i,j,:))/(ro_zm(i,j)*rsamples)
        u222_zm(i,j)=integration(u222(i,j,:))/(ro_zm(i,j)*rsamples)
        u333_zm(i,j)=integration(u333(i,j,:))/(ro_zm(i,j)*rsamples)
        u112_zm(i,j)=integration(u112(i,j,:))/(ro_zm(i,j)*rsamples)
        u113_zm(i,j)=integration(u113(i,j,:))/(ro_zm(i,j)*rsamples)
        u122_zm(i,j)=integration(u122(i,j,:))/(ro_zm(i,j)*rsamples)
        u133_zm(i,j)=integration(u133(i,j,:))/(ro_zm(i,j)*rsamples)
        u223_zm(i,j)=integration(u223(i,j,:))/(ro_zm(i,j)*rsamples)
        u233_zm(i,j)=integration(u233(i,j,:))/(ro_zm(i,j)*rsamples)
        u123_zm(i,j)=integration(u123(i,j,:))/(ro_zm(i,j)*rsamples)
        
         tt_zm(i,j) = integration(tt(i,j,:))/(ro_zm(i,j)*rsamples)
        tu1_zm(i,j) = integration(tu1(i,j,:))/(ro_zm(i,j)*rsamples)
        tu2_zm(i,j) = integration(tu2(i,j,:))/(ro_zm(i,j)*rsamples)
        tu3_zm(i,j) = integration(tu3(i,j,:))/(ro_zm(i,j)*rsamples)
        !
        pp_zm(i,j) = integration(pp(i,j,:))/rsamples
        !
      end do
      end do
      !
      print*,' ** averaging along k direction ... done.'
      pp_zm= pp_zm-p_zm**2
      !
      u11_zm=u11_zm-u1_zm**2
      u22_zm=u22_zm-u2_zm**2
      u33_zm=u33_zm-u3_zm**2
      u12_zm=u12_zm-u1_zm*u2_zm
      u13_zm=u13_zm-u1_zm*u3_zm
      u23_zm=u23_zm-u2_zm*u3_zm
       tt_zm= tt_zm-t_zm**2
      tu1_zm=tu1_zm-u1_zm*t_zm
      tu2_zm=tu2_zm-u2_zm*t_zm
      tu3_zm=tu3_zm-u3_zm*t_zm
      !
      u111_zm=u111_zm-u1_zm**3-3.d0*u11_zm*u1_zm
      u222_zm=u222_zm-u2_zm**3-3.d0*u22_zm*u2_zm
      u333_zm=u333_zm-u3_zm**3-3.d0*u33_zm*u3_zm
      u112_zm=u112_zm-u1_zm**2*u2_zm-u11_zm*u2_zm-2.d0*u12_zm*u1_zm
      u113_zm=u113_zm-u1_zm**2*u3_zm-u11_zm*u3_zm-2.d0*u13_zm*u1_zm
      u122_zm=u122_zm-u1_zm*u2_zm**2-u22_zm*u1_zm-2.d0*u12_zm*u2_zm
      u133_zm=u133_zm-u1_zm*u3_zm**2-u33_zm*u1_zm-2.d0*u13_zm*u3_zm
      u223_zm=u223_zm-u3_zm*u2_zm**2-u22_zm*u3_zm-2.d0*u23_zm*u2_zm
      u233_zm=u233_zm-u2_zm*u3_zm**2-u33_zm*u2_zm-2.d0*u23_zm*u3_zm
      u123_zm=u123_zm-u1_zm*u2_zm*u3_zm-u12_zm*u3_zm-u13_zm*u2_zm-u23_zm*u1_zm
      !
      inquire(file='Results/2order.fav.zm.h5',exist=filexists)
      !
      if(filexists) then
        call system('mv -v Results/2order.fav.zm.h5 Results/2order.fav.zm.bak')
      endif
      !
      call H5WriteArray(u11_zm,im,jm,'u11','Results/2order.fav.zm.h5')
      call H5WriteArray(u22_zm,im,jm,'u22','Results/2order.fav.zm.h5')
      call H5WriteArray(u33_zm,im,jm,'u33','Results/2order.fav.zm.h5')
      call H5WriteArray(u12_zm,im,jm,'u12','Results/2order.fav.zm.h5')
      call H5WriteArray(u13_zm,im,jm,'u13','Results/2order.fav.zm.h5')
      call H5WriteArray(u23_zm,im,jm,'u23','Results/2order.fav.zm.h5')
      call H5WriteArray( tt_zm,im,jm, 'tt','Results/2order.fav.zm.h5')
      call H5WriteArray(tu1_zm,im,jm,'tu1','Results/2order.fav.zm.h5')
      call H5WriteArray(tu2_zm,im,jm,'tu2','Results/2order.fav.zm.h5')
      call H5WriteArray(tu3_zm,im,jm,'tu3','Results/2order.fav.zm.h5')
      call H5WriteArray( pp_zm,im,jm, 'pp','Results/2order.fav.zm.h5')
      
      inquire(file='Results/3order.fav.zm.h5',exist=filexists)
      !
      if(filexists) then
        call system('mv -v Results/3order.fav.zm.h5 Results/3order.fav.zm.bak')
      endif
      !
      call H5WriteArray(u111_zm,im,jm,'u111','Results/3order.fav.zm.h5')
      call H5WriteArray(u222_zm,im,jm,'u222','Results/3order.fav.zm.h5')
      call H5WriteArray(u333_zm,im,jm,'u333','Results/3order.fav.zm.h5')
      call H5WriteArray(u112_zm,im,jm,'u112','Results/3order.fav.zm.h5')
      call H5WriteArray(u113_zm,im,jm,'u113','Results/3order.fav.zm.h5')
      call H5WriteArray(u122_zm,im,jm,'u122','Results/3order.fav.zm.h5')
      call H5WriteArray(u133_zm,im,jm,'u133','Results/3order.fav.zm.h5')
      call H5WriteArray(u223_zm,im,jm,'u223','Results/3order.fav.zm.h5')
      call H5WriteArray(u233_zm,im,jm,'u233','Results/3order.fav.zm.h5')
      call H5WriteArray(u123_zm,im,jm,'u123','Results/3order.fav.zm.h5')
      !
      call writetecbin('Results/tec2order.plt',x(:,:,0),'x',y(:,:,0),'y', &
                               u11_zm,'uu',u22_zm,'vv',u33_zm,'ww',    &
                               u12_zm,'uv',u13_zm,'uw',u23_zm,'vw',im,jm)
      !
      deallocate(u11_zm,u12_zm,u13_zm,u22_zm,u23_zm,u33_zm,u111_zm,    &
                 u112_zm,u113_zm,u122_zm,u123_zm,u133_zm,u222_zm,      &
                 u223_zm,u233_zm,u333_zm,tt_zm,pp_zm,tu1_zm,tu2_zm,tu3_zm)
      !
    endif
    !
    ! return
    !
    rsamples=real(nsamples,8)
    !
    u11=u11/ro_m/rsamples
    u12=u12/ro_m/rsamples
    u13=u13/ro_m/rsamples
    u22=u22/ro_m/rsamples
    u23=u23/ro_m/rsamples
    u33=u33/ro_m/rsamples
    !
    u111=u111/ro_m/rsamples
    u222=u222/ro_m/rsamples
    u333=u333/ro_m/rsamples
    u112=u112/ro_m/rsamples
    u113=u113/ro_m/rsamples
    u122=u122/ro_m/rsamples
    u133=u133/ro_m/rsamples
    u223=u223/ro_m/rsamples
    u233=u233/ro_m/rsamples
    u123=u123/ro_m/rsamples
    !
     tt =  tt/ro_m/rsamples
    tu1 = tu1/ro_m/rsamples
    tu2 = tu2/ro_m/rsamples
    tu3 = tu3/ro_m/rsamples
    !
    pp = pp/rsamples
    !
    pp= pp-p_m**2
    !
    u11=u11-u1_m**2
    u22=u22-u2_m**2
    u33=u33-u3_m**2
    u12=u12-u1_m*u2_m
    u13=u13-u1_m*u3_m
    u23=u23-u2_m*u3_m
     tt= tt-t_m**2
    tu1=tu1-u1_m*t_m
    tu2=tu2-u2_m*t_m
    tu3=tu3-u3_m*t_m
    !
    u111=u111-u1_m**3-3.d0*u11*u1_m
    u222=u222-u2_m**3-3.d0*u22*u2_m
    u333=u333-u3_m**3-3.d0*u33*u3_m
    u112=u112-u1_m**2*u2_m-u11*u2_m-2.d0*u12*u1_m
    u113=u113-u1_m**2*u3_m-u11*u3_m-2.d0*u13*u1_m
    u122=u122-u1_m*u2_m**2-u22*u1_m-2.d0*u12*u2_m
    u133=u133-u1_m*u3_m**2-u33*u1_m-2.d0*u13*u3_m
    u223=u223-u3_m*u2_m**2-u22*u3_m-2.d0*u23*u2_m
    u233=u233-u2_m*u3_m**2-u33*u2_m-1.d0*u23*u3_m
    u123=u123-u1_m*u2_m*u3_m-u12*u3_m-u13*u2_m-u23*u1_m
    ! 
    inquire(file='Results/2order.fav.h5',exist=filexists)
    !
    if(filexists) then
      call system('mv -v Results/2order.fav.h5 Results/2order.fav.bak')
    endif
    !
    call H5WriteArray(u11,im,jm,km,'u11','Results/2order.fav.h5')
    call H5WriteArray(u22,im,jm,km,'u22','Results/2order.fav.h5')
    call H5WriteArray(u33,im,jm,km,'u33','Results/2order.fav.h5')
    call H5WriteArray(u12,im,jm,km,'u12','Results/2order.fav.h5')
    call H5WriteArray(u13,im,jm,km,'u13','Results/2order.fav.h5')
    call H5WriteArray(u23,im,jm,km,'u23','Results/2order.fav.h5')
    call H5WriteArray( tt,im,jm,km, 'tt','Results/2order.fav.h5')
    call H5WriteArray(tu1,im,jm,km,'tu1','Results/2order.fav.h5')
    call H5WriteArray(tu2,im,jm,km,'tu2','Results/2order.fav.h5')
    call H5WriteArray(tu3,im,jm,km,'tu3','Results/2order.fav.h5')
    call H5WriteArray( pp,im,jm,km, 'pp','Results/2order.fav.h5')
    !
    inquire(file='Results/3order.fav.h5',exist=filexists)
    !
    if(filexists) then
      call system('mv -v Results/3order.fav.h5 Results/3order.fav.bak')
    endif
    !
    call H5WriteArray(u111,im,jm,km,'u111','Results/3order.fav.h5')
    call H5WriteArray(u222,im,jm,km,'u222','Results/3order.fav.h5')
    call H5WriteArray(u333,im,jm,km,'u333','Results/3order.fav.h5')
    call H5WriteArray(u112,im,jm,km,'u112','Results/3order.fav.h5')
    call H5WriteArray(u113,im,jm,km,'u113','Results/3order.fav.h5')
    call H5WriteArray(u122,im,jm,km,'u122','Results/3order.fav.h5')
    call H5WriteArray(u133,im,jm,km,'u133','Results/3order.fav.h5')
    call H5WriteArray(u223,im,jm,km,'u223','Results/3order.fav.h5')
    call H5WriteArray(u233,im,jm,km,'u233','Results/3order.fav.h5')
    call H5WriteArray(u123,im,jm,km,'u123','Results/3order.fav.h5')
    !
    deallocate(u11,u12,u13,u22,u23,u33,u111,u112,u113,u122,u123,u133,  &
               u222,u223,u233,u333,tt,pp,tu1,tu2,tu3)
    !
  end subroutine highordermom
  !+-------------------------------------------------------------------+
  !| The end of the subroutine meanflow.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| The subroutine is used to calculate the gradient of mean flow.    |
  !+-------------------------------------------------------------------+
  !| Writen by Jian Fang, 2020-04-18                                   |
  !+-------------------------------------------------------------------+
  subroutine dutmean
    !
    use commvardefine, only: im,jm,km,x,y,z,u1_m,u2_m,u3_m,t_m,u1_zm,  &
                             u2_zm,u3_zm,t_zm,du1_m,du2_m,du3_m,dt_m,  &
                             du1_zm,du2_zm,du3_zm,dt_zm,               &
                             du1_xzm,du2_xzm,du3_xzm,dt_xzm,           &
                             lihomo,ljhomo,lkhomo
    use gradsolver
    use basicfunction, only: integration
    use h5readwrite
    !
    integer :: i,j,k
    !
    if(km>0) then
      allocate(du1_m(1:3,0:im,0:jm,0:km),du2_m(1:3,0:im,0:jm,0:km),      &
               du3_m(1:3,0:im,0:jm,0:km), dt_m(1:3,0:im,0:jm,0:km))
      !
      du1_m=grad(u1_m,x,y,z) 
      du2_m=grad(u2_m)
      du3_m=grad(u3_m)
       dt_m=grad(t_m)
      !
      call H5WriteArray(du1_m(1,:,:,:),im,jm,km,'du1dx','Results/mean.fav.h5')
      call H5WriteArray(du1_m(2,:,:,:),im,jm,km,'du1dy','Results/mean.fav.h5')
      call H5WriteArray(du1_m(3,:,:,:),im,jm,km,'du1dz','Results/mean.fav.h5')
      call H5WriteArray(du2_m(1,:,:,:),im,jm,km,'du2dx','Results/mean.fav.h5')
      call H5WriteArray(du2_m(2,:,:,:),im,jm,km,'du2dy','Results/mean.fav.h5')
      call H5WriteArray(du2_m(3,:,:,:),im,jm,km,'du2dz','Results/mean.fav.h5')
      call H5WriteArray(du3_m(1,:,:,:),im,jm,km,'du3dx','Results/mean.fav.h5')
      call H5WriteArray(du3_m(2,:,:,:),im,jm,km,'du3dy','Results/mean.fav.h5')
      call H5WriteArray(du3_m(3,:,:,:),im,jm,km,'du3dz','Results/mean.fav.h5')
      call H5WriteArray( dt_m(1,:,:,:),im,jm,km, 'dtdx','Results/mean.fav.h5')
      call H5WriteArray( dt_m(2,:,:,:),im,jm,km, 'dtdy','Results/mean.fav.h5')
      call H5WriteArray( dt_m(3,:,:,:),im,jm,km, 'dtdz','Results/mean.fav.h5')
    endif
    !
    if(lihomo .and. lkhomo) then
      !
      allocate(du1_xzm(1:3,0:jm),du2_xzm(1:3,0:jm),du3_xzm(1:3,0:jm),  &
               dt_xzm(1:3,0:jm))
      !
      du1_xzm=0.d0
      du2_xzm=0.d0
      du3_xzm=0.d0
       dt_xzm=0.d0
      !
      do j=0,jm
        du1_xzm(1,j)=integration(du1_m(1,:,j,:))
        du1_xzm(2,j)=integration(du1_m(2,:,j,:))
        du1_xzm(3,j)=integration(du1_m(3,:,j,:))
        du2_xzm(1,j)=integration(du2_m(1,:,j,:))
        du2_xzm(2,j)=integration(du2_m(2,:,j,:))
        du2_xzm(3,j)=integration(du2_m(3,:,j,:))
        du3_xzm(1,j)=integration(du3_m(1,:,j,:))
        du3_xzm(2,j)=integration(du3_m(2,:,j,:))
        du3_xzm(3,j)=integration(du3_m(3,:,j,:))
         dt_xzm(1,j)=integration(dt_m(1,:,j,:))
         dt_xzm(2,j)=integration(dt_m(2,:,j,:))
         dt_xzm(3,j)=integration(dt_m(3,:,j,:))
      enddo
      !
      du1_xzm=du1_xzm/real(im*km,8)
      du2_xzm=du2_xzm/real(im*km,8)
      du3_xzm=du3_xzm/real(im*km,8)
       dt_xzm= dt_xzm/real(im*km,8)
      !
      call H5WriteArray(du1_xzm(1,:),jm,'du1dx','Results/mean.fav.xzm.h5')
      call H5WriteArray(du1_xzm(2,:),jm,'du1dy','Results/mean.fav.xzm.h5')
      call H5WriteArray(du1_xzm(3,:),jm,'du1dz','Results/mean.fav.xzm.h5')
      call H5WriteArray(du2_xzm(1,:),jm,'du2dx','Results/mean.fav.xzm.h5')
      call H5WriteArray(du2_xzm(2,:),jm,'du2dy','Results/mean.fav.xzm.h5')
      call H5WriteArray(du2_xzm(3,:),jm,'du2dz','Results/mean.fav.xzm.h5')
      call H5WriteArray(du3_xzm(1,:),jm,'du3dx','Results/mean.fav.xzm.h5')
      call H5WriteArray(du3_xzm(2,:),jm,'du3dy','Results/mean.fav.xzm.h5')
      call H5WriteArray(du3_xzm(3,:),jm,'du3dz','Results/mean.fav.xzm.h5')
      call H5WriteArray( dt_xzm(1,:),jm, 'dtdx','Results/mean.fav.xzm.h5')
      call H5WriteArray( dt_xzm(2,:),jm, 'dtdy','Results/mean.fav.xzm.h5')
      call H5WriteArray( dt_xzm(3,:),jm, 'dtdz','Results/mean.fav.xzm.h5')
      !
    endif
    if(lkhomo) then
      !
      allocate(du1_zm(1:3,0:im,0:jm),du2_zm(1:3,0:im,0:jm),              &
               du3_zm(1:3,0:im,0:jm), dt_zm(1:3,0:im,0:jm))
      !
      du1_zm=0.d0
      du2_zm=0.d0
      du3_zm=0.d0
       dt_zm=0.d0
      !
      du1_zm(1:2,:,:)=grad_xy(u1_zm,x,y) 
      du2_zm(1:2,:,:)=grad(u2_zm)
      du3_zm(1:2,:,:)=grad(u3_zm)
       dt_zm(1:2,:,:)=grad(t_zm)
      !
      ! do j=0,jm
      ! do i=0,im
      !   du1_zm(1,i,j)=integration(du1_m(1,i,j,:))
      !   du1_zm(2,i,j)=integration(du1_m(2,i,j,:))
      !   du1_zm(3,i,j)=integration(du1_m(3,i,j,:))
      !   du2_zm(1,i,j)=integration(du2_m(1,i,j,:))
      !   du2_zm(2,i,j)=integration(du2_m(2,i,j,:))
      !   du2_zm(3,i,j)=integration(du2_m(3,i,j,:))
      !   du3_zm(1,i,j)=integration(du3_m(1,i,j,:))
      !   du3_zm(2,i,j)=integration(du3_m(2,i,j,:))
      !   du3_zm(3,i,j)=integration(du3_m(3,i,j,:))
      !    dt_zm(1,i,j)=integration(dt_m(1,i,j,:))
      !    dt_zm(2,i,j)=integration(dt_m(2,i,j,:))
      !    dt_zm(3,i,j)=integration(dt_m(3,i,j,:))
      ! enddo
      ! enddo
      ! !
      ! if(km>0) then
      !   du1_zm=du1_zm/dble(km)
      !   du2_zm=du2_zm/dble(km)
      !   du3_zm=du3_zm/dble(km)
      !    dt_zm= dt_zm/dble(km)
      ! endif
      !
      call H5WriteArray(du1_zm(1,:,:),im,jm,'du1dx','Results/mean.fav.zm.h5')
      call H5WriteArray(du1_zm(2,:,:),im,jm,'du1dy','Results/mean.fav.zm.h5')
      call H5WriteArray(du1_zm(3,:,:),im,jm,'du1dz','Results/mean.fav.zm.h5')
      call H5WriteArray(du2_zm(1,:,:),im,jm,'du2dx','Results/mean.fav.zm.h5')
      call H5WriteArray(du2_zm(2,:,:),im,jm,'du2dy','Results/mean.fav.zm.h5')
      call H5WriteArray(du2_zm(3,:,:),im,jm,'du2dz','Results/mean.fav.zm.h5')
      call H5WriteArray(du3_zm(1,:,:),im,jm,'du3dx','Results/mean.fav.zm.h5')
      call H5WriteArray(du3_zm(2,:,:),im,jm,'du3dy','Results/mean.fav.zm.h5')
      call H5WriteArray(du3_zm(3,:,:),im,jm,'du3dz','Results/mean.fav.zm.h5')
      call H5WriteArray( dt_zm(1,:,:),im,jm, 'dtdx','Results/mean.fav.zm.h5')
      call H5WriteArray( dt_zm(2,:,:),im,jm, 'dtdy','Results/mean.fav.zm.h5')
      call H5WriteArray( dt_zm(3,:,:),im,jm, 'dtdz','Results/mean.fav.zm.h5')
      !
    endif
    !
  end subroutine dutmean
  !+-------------------------------------------------------------------+
  !| The end of the subroutine dutmean.                                |
  !+-------------------------------------------------------------------+
  !
  subroutine flowdataclean(nsta,nend)
    !
    use commvardefine, only : im,jm,km
    use h5readwrite
    !
    integer,intent(in) :: nsta,nend
    !
    integer :: n,nstep
    real(8) :: time
    character(len=4) sname
    character(len=23) fname
    !
    real(8),allocatable,dimension(:,:,:) :: u1,u2,u3,p,t,ro
    !
    allocate(u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),u3(0:im,0:jm,0:km), &
             p(0:im,0:jm,0:km),t(0:im,0:jm,0:km),ro(0:im,0:jm,0:km))

    do n=nsta,nend
      !
      write(sname,'(i4.4)')n
      fname=outfolder//'flowfield'//sname//'.h5'
      call H5ReadArray(time,'time',fname)
      call H5ReadArray(nstep,'nstep',fname)
      write(*,'(A,F10.5,A,I0)')'  ** << flow field  at t=',time,' nstep=',nstep
      !
      call H5ReadArray(ro,im,jm,km,'ro',fname)
      call H5ReadArray(u1,im,jm,km,'u1',fname)
      call H5ReadArray(u2,im,jm,km,'u2',fname)
      call H5ReadArray(u3,im,jm,km,'u3',fname)
      call H5ReadArray(p,im,jm,km,  'p',fname)
      call H5ReadArray(t,im,jm,km,  't',fname)
      !
      call system('rm -v '//fname)
      !
      call H5WriteArray(time,'time',fname)
      call H5WriteArray(nstep,'nstep',fname)
      call H5WriteArray(ro,im,jm,km,'ro',fname)
      call H5WriteArray(u1,im,jm,km,'u1',fname)
      call H5WriteArray(u2,im,jm,km,'u2',fname)
      call H5WriteArray(u3,im,jm,km,'u3',fname)
      call H5WriteArray(p,im,jm,km,  'p',fname)
      call H5WriteArray(t,im,jm,km,  't',fname)
      !
    enddo
    !
  end subroutine flowdataclean
  !
  subroutine dataextra(nsta,nend)
    !
    use commvardefine, only : im,jm,km
    use h5readwrite
    !
    integer,intent(in) :: nsta,nend
    !
    integer :: n,nstep
    real(8) :: time
    character(len=4) sname
    character(len=23) fname
    !
    real(8),allocatable,dimension(:,:,:) :: u1,u2,u3,p,t,ro
    !
    allocate(u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),u3(0:im,0:jm,0:km), &
             p(0:im,0:jm,0:km),t(0:im,0:jm,0:km),ro(0:im,0:jm,0:km))

    do n=nsta,nend
      !
      write(sname,'(i4.4)')n
      fname=outfolder//'flowfield'//sname//'.h5'
      call H5ReadArray(time,'time',fname)
      call H5ReadArray(nstep,'nstep',fname)
      write(*,'(A,F10.5,A,I0)')'  ** << flow field  at t=',time,' nstep=',nstep
      !
      call H5ReadArray(ro,im,jm,km,'ro',fname)
      ! call H5ReadArray(u1,im,jm,km,'u1',fname)
      ! call H5ReadArray(u2,im,jm,km,'u2',fname)
      ! call H5ReadArray(u3,im,jm,km,'u3',fname)
      ! call H5ReadArray(p,im,jm,km,  'p',fname)
      ! call H5ReadArray(t,im,jm,km,  't',fname)
      !
      ! call system('rm -v '//fname)
      !
      open(18,file='data/density'//sname//'.dat',form='unformatted')
      write(18)im,jm,km
      write(18)time
      write(18)ro(1:im,1:jm,1:km)
      close(18)
      print*,' << data/density',sname,'.dat'
      ! call H5WriteArray(time,'time',fname)
      ! call H5WriteArray(nstep,'nstep',fname)
      ! call H5WriteArray(ro,im,jm,km,'ro',fname)
      ! call H5WriteArray(u1,im,jm,km,'u1',fname)
      ! call H5WriteArray(u2,im,jm,km,'u2',fname)
      ! call H5WriteArray(u3,im,jm,km,'u3',fname)
      ! call H5WriteArray(p,im,jm,km,  'p',fname)
      ! call H5WriteArray(t,im,jm,km,  't',fname)
      !
    enddo
    !
  end subroutine dataextra
  !
  subroutine vortsta(nsta,nend)
    !
    use commvardefine, only : im,jm,km,x,y,z
    use readwrite,     only : H5_ReadGrid
    use basicfunction, only : omega
    use gradsolver,    only : grad
    use h5readwrite
    use writetec
    !
    integer,intent(in) :: nsta,nend
    !
    integer :: n,i,j,k
    real(8) :: time
    logical :: filexists=.false.
    character(len=4) sname
    character(len=23) fname
    !
    real(8),allocatable,dimension(:,:,:) :: u1,u2,u3,omegax,omegay,omegaz
    real(8),allocatable,dimension(:,:,:,:) :: du1,du2,du3
    real(8),allocatable,dimension(:,:) :: omega11,omega22,omega33,omega12, &
                                          omega13,omega23,omegaxm,omegaym, &
                                          omegazm
    !
    allocate(u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),u3(0:im,0:jm,0:km), &
             omegax(0:im,0:jm,0:km),omegay(0:im,0:jm,0:km),omegaz(0:im,0:jm,0:km))
    allocate(du1(1:3,0:im,0:jm,0:km),du2(1:3,0:im,0:jm,0:km),          &
             du3(1:3,0:im,0:jm,0:km))
    allocate(omega11(0:im,0:jm),omega22(0:im,0:jm),omega33(0:im,0:jm), &
             omega12(0:im,0:jm),omega13(0:im,0:jm),omega23(0:im,0:jm), &
             omegaxm(0:im,0:jm),omegaym(0:im,0:jm),omegazm(0:im,0:jm))
    !
    call H5_ReadGrid
    !
    omegaxm=0.d0
    omegaym=0.d0
    omegazm=0.d0
    !
    omega11=0.d0
    omega22=0.d0
    omega33=0.d0
    omega12=0.d0
    omega13=0.d0
    omega23=0.d0
    !
    do n=nsta,nend
      !
      write(sname,'(i4.4)')n
      fname=outfolder//'flowfield'//sname//'.h5'
      call H5ReadArray(time,'time',fname)
      write(*,'(A,F10.5)')'  ** << flow field  at t=',time
      !
      call H5ReadArray(u1,im,jm,km,'u1',fname)
      call H5ReadArray(u2,im,jm,km,'u2',fname)
      call H5ReadArray(u3,im,jm,km,'u3',fname)
      !
      du1=grad(u1,x,y,z)
      du2=grad(u2)
      du3=grad(u3)
      !
      omegax=omega(du1,du2,du3,'x')
      omegay=omega(du1,du2,du3,'y')
      omegaz=omega(du1,du2,du3,'z')
      !
      do k=1,km
        !
        do j=0,jm
        do i=0,im
          !
          omegaxm(i,j)=omegaxm(i,j)+omegax(i,j,k)
          omegaym(i,j)=omegaym(i,j)+omegay(i,j,k)
          omegazm(i,j)=omegazm(i,j)+omegaz(i,j,k)
          !
          omega11(i,j)=omega11(i,j)+omegax(i,j,k)*omegax(i,j,k)
          omega22(i,j)=omega22(i,j)+omegay(i,j,k)*omegay(i,j,k)
          omega33(i,j)=omega33(i,j)+omegaz(i,j,k)*omegaz(i,j,k)
          omega12(i,j)=omega12(i,j)+omegax(i,j,k)*omegay(i,j,k)
          omega13(i,j)=omega13(i,j)+omegax(i,j,k)*omegaz(i,j,k)
          omega23(i,j)=omega23(i,j)+omegay(i,j,k)*omegaz(i,j,k)
        enddo
        enddo
        !
      enddo
      !
    enddo
    !
    omegaxm=omegaxm/dble((nend-nsta+1)*km)
    omegaym=omegaym/dble((nend-nsta+1)*km)
    omegazm=omegazm/dble((nend-nsta+1)*km)
    !
    omega11=omega11/dble((nend-nsta+1)*km)
    omega22=omega22/dble((nend-nsta+1)*km)
    omega33=omega33/dble((nend-nsta+1)*km)
    omega12=omega12/dble((nend-nsta+1)*km)
    omega13=omega13/dble((nend-nsta+1)*km)
    omega23=omega23/dble((nend-nsta+1)*km)
    !
    call writetecbin('tecomega.plt',x(:,:,0),'x',y(:,:,0),'y',  &
                                                 omegaxm,'omegax',  &
                                                 omegaym,'omegay',  &
                                                 omegazm,'omegaz',  &
                                                 omega11,'omega11',  &
                                                 omega22,'omega22',  &
                                                 omega33,'omega33',  &
                                                 omega12,'omega12',  &
                                                 omega13,'omega13',  &
                                                 omega23,'omega23',im,jm)
    !
    open(16,file='omega.bin',form='unformatted')
    write(16)im,jm
    write(16)x(:,:,0),y(:,:,0)
    write(16)omegaxm,omegaym,omegazm
    write(16)omega11,omega22,omega33
    write(16)omega12,omega13,omega23
    close(16)
    print*,' << omega.bin'
    !
  end subroutine vortsta
  !
end module statistic
!+---------------------------------------------------------------------+
!| The end of the module statistic.                                    |
!+---------------------------------------------------------------------+