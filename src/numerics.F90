module numerics
  !
  use commvardefine, only : pi
  use fdm
  !
  implicit none
  !
  integer :: i
  !
  real(8),allocatable :: WaveK(:)
  complex(8),allocatable :: WaveMK(:)
  real(8) :: dissp
  complex(8),parameter :: cunit=(0.d0,1.d0)
  !
  contains
  !
  subroutine wave_propa
    !
    integer :: num_node
    real(8),allocatable,dimension(:) :: x,f,fh,df,f_exact,rhs
    !
    integer :: i,nstep,irk,nfile
    real(8) :: dx,time,deltat,cpeed,error
    character(len=8) :: scheme
    character(len=2) :: fname
    !
    scheme='cmp5ld'
    !
    if(trim(scheme)=='eu5') then
      print*,' ** reconstruction scheme: 5th-order explicit upwind'
    elseif(trim(scheme)=='cmp5') then
      print*,' ** reconstruction scheme: CMP5'
    elseif(trim(scheme)=='cmp5ld') then
      print*,' ** reconstruction scheme: CMP5-LD'
    elseif(trim(scheme)=='mp5') then
      print*,' ** reconstruction scheme: MP5'
    elseif(trim(scheme)=='mp5ld') then
      print*,' ** reconstruction scheme: MP5-LD'
    elseif(trim(scheme)=='eu5-ld') then
      print*,' ** reconstruction scheme: low-dissipative 5th-order explicit upwind'
    elseif(trim(scheme)=='cu5') then
      print*,' ** reconstruction scheme: 5th-order compact upwind'
    elseif(trim(scheme)=='cu5-ld') then
      print*,' ** reconstruction scheme: low-dissipative 5th-order compact upwind'
    else
      print*,scheme
      stop ' !! scheme not defined @ errana1d_multi_domain !!'
    endif
    !
    ! cpeed=0.5d0
    cpeed=1.d0
    !
    ! num_node=664
    num_node=100
    !
    allocate(x(0:num_node),f(-4:num_node+4),df(0:num_node), &
             rhs(0:num_node),f_exact(0:num_node))
    !
    do i=0,num_node
      x(i)=2.d0/dble(num_node)*dble(i)-1.d0
      ! x(i)=2.d0*pi/dble(num_node)*dble(i)
    enddo
    dx=x(1)-x(0)
    !
    do i=0,num_node
      f(i)=hybridwave(x(i))
      ! f(i)=wavepacket(x(i))
      ! exp(-50.d0*(x(i)-1.5)**2)*sin(0.838242d0/dx*x(i))
      ! f(i)=sin(x(i))
    enddo
    !
    open(16,file='f_init.dat')
    do i=0,num_node
      write(16,*)x(i),f(i)
    enddo
    close(16)
    print*,' << f_init.dat'
    !
    time=0.d0
    deltat=1.d-4
    !
    nstep=0
    nfile=0
    !
    open(18,file='error_wavepack_'//trim(scheme)//'.dat')
    do while(nstep<20000)
      !
      nstep=nstep+1
      time=time+deltat
      !
      do irk=1,4
        !
        ! periodic condition
        f(-4:-1)=f(num_node-4:num_node-1)
        f(num_node+1:num_node+4)=f(1:4)
        !
        df=dfcal(f,num_node,trim(scheme))/dx
        !
        rhs=-cpeed*df
        !
        call rk4(f(0:num_node),rhs,deltat,irk)
        !
      enddo
      !
      print*,nstep,time
      !
      ! if(mod(nstep,1000)==0) then
      !   !
      !   do i=0,num_node
      !     f_exact(i)=wavepacket(x(i)-time*cpeed)
      !   enddo
      !   !
      !   error=0.d0
      !   do i=1,num_node
      !     error=error+(f(i)-f_exact(i))**2
      !   enddo
      !   error=sqrt(error/dble(num_node))
      !   !
      !   write(18,*)nstep,time,error
      !   !
      ! endif
      !
      if(mod(nstep,10000)==0) then
        !
        nfile=nfile+1
        write(fname,'(i2.2)')nfile
        !
        do i=0,num_node
          ! f_exact(i)=wavepacket(x(i)-time*cpeed)
          f_exact(i)=hybridwave(x(i)-time*cpeed)
        enddo
        !
        open(16,file='fsolution_'//trim(scheme)//'_'//fname//'.dat')
        do i=0,num_node
          write(16,*)x(i),f(i),f_exact(i)
        enddo
        close(16)
        print*,' << fsolution',fname,'.dat'
        !
      endif
      !
    enddo
    close(18)
    print*,' <<'//'error_wavepack_'//trim(scheme)//'.dat'
    !
    !
    contains
    !
    function hybridwave(xin) result(fpro)
      !
      real(8),intent(in) :: xin
      real(8) ::  fpro
      !
      real(8) :: a,alfa,beter,zeter,deter
      real(8) :: Gfc1,Gfc2,Gfc3,Ffc1,Ffc2,Ffc3
      !
      a=0.5d0
      zeter=-0.7d0
      deter=0.005d0
      alfa=10.d0
      beter=dlog(2.d0)/(36.d0*deter**2)
      !
      if(xin>=-0.8d0 .and. xin<=-0.6d0) then
        Gfc1=dexp(-beter*(xin-zeter)**2)
        Gfc2=dexp(-beter*(xin-(zeter+deter))**2)
        Gfc3=dexp(-beter*(xin-(zeter-deter))**2)
        fpro=1.d0/6.0*(Gfc3+Gfc2+4.d0*Gfc1)
      elseif(xin>=-0.4d0 .and. xin<=-0.2d0) then
        fpro=1.d0
      elseif(xin>=0.d0 .and. xin<=0.2d0) then
        fpro=1.d0-dabs(10.d0*(xin-0.1d0))
      elseif(xin>=0.4d0 .and. xin<=0.6d0) then
        Ffc1=1.d0-(alfa**2)*((xin-a)**2)
        Ffc1=dsqrt(max(Ffc1,0.d0))
        Ffc2=1.d0-(alfa**2)*((xin-a-deter)**2)
        Ffc2=dsqrt(max(Ffc2,0.d0))
        Ffc3=1.d0-(alfa**2)*((xin-a+deter)**2)
        Ffc3=dsqrt(max(Ffc3,0.d0))
        fpro=1.d0/6.0*(Ffc3+Ffc2+4.d0*Ffc1)
      else
        fpro=0.d0
      end if
      !
    end function hybridwave
    !
    function wavepacket(xin) result(fpro)
      !
      real(8),intent(in) :: xin
      real(8) ::  fpro
      !
      fpro=exp(-50.d0*(xin-1.5d0)**2)*sin(0.838242d0/dx*xin)
      !
      return
      !
    end function wavepacket
    !
  end subroutine wave_propa
  !
  subroutine errana1d_multi_domain
    !
    integer :: num_node,num_rank,max_rank,tol_node
    real(8),allocatable,dimension(:,:) :: x,f,fh,df,df_exact
    !
    integer :: i,jrank,jrplus,jrmius,ncorr
    real(8) :: error
    character(len=8) :: scheme
    character(len=2) :: rkname
    !
    scheme='cu5'
    !
    if(trim(scheme)=='cu5') then
      print*,' ** reconstruction scheme: 5th-order compact upwind'
    elseif(trim(scheme)=='cu5-ld') then
      print*,' ** reconstruction scheme: low-dissipative 5th-order compact upwind'
    else
      print*,scheme
      stop ' !! scheme not defined @ errana1d_multi_domain !!'
    endif
    !
    num_rank=4
    write(rkname,'(i2.2)')num_rank
    max_rank=num_rank-1
    !
    tol_node=8
    !
    open(18,file='error_'//trim(scheme)//'_rk_'//rkname//'.dat')
    !
    do while(tol_node<=512)
      !
      ! num_rank=num_rank*2
      ! max_rank=num_rank-1
      !
      tol_node=tol_node*2
      num_node=tol_node/num_rank
      !
      allocate(x(0:max_rank,0:num_node),f(0:max_rank,-4:num_node+4),   &
               fh(0:max_rank,-1:num_node),df(0:max_rank,0:num_node),   &
               df_exact(0:max_rank,0:num_node))
      !
      do jrank=0,max_rank
        do i=0,num_node
          x(jrank,i)=2.d0*pi/dble(tol_node)*dble(i+jrank*num_node)
          !
          f(jrank,i)=sin(4.d0*x(jrank,i))
          df_exact(jrank,i)=4.d0*cos(4.d0*x(jrank,i))
        enddo
      enddo
      !
      ! data swap
      do jrank=0,max_rank
        !
        jrplus=jrank+1
        jrmius=jrank-1
        !
        if(jrplus>max_rank)  jrplus=0
        if(jrmius<0)         jrmius=max_rank
        !
        f(jrank,-4:-1)=f(jrmius,num_node-4:num_node-1)
        f(jrank,num_node+1:num_node+4)=f(jrplus,1:4)
        !
      enddo
      !
      do jrank=0,max_rank
        call recons(f(jrank,:),fh(jrank,:),num_node,scheme=trim(scheme),mode='parallel',cstep='pred')
      enddo
      !
      ! ! data swap
      ! do jrank=0,max_rank
      !   !
      !   jrplus=jrank+1
      !   jrmius=jrank-1
      !   !
      !   if(jrplus>max_rank)  jrplus=0
      !   if(jrmius<0)         jrmius=max_rank
      !   !
      !   fh(jrank,-1)=fh(jrmius,num_node-1)
      !   fh(jrank,num_node)=fh(jrplus,0)
      !   !
      ! enddo
      !
      ncorr=1
      do while(ncorr<=0)
        !
        ncorr=ncorr+1
        !
        do jrank=0,max_rank
          call recons(f(jrank,:),fh(jrank,:),num_node,scheme=trim(scheme),mode='parallel',cstep='corr')
        enddo
        !
        ! data swap
        do jrank=0,max_rank
          !
          jrplus=jrank+1
          jrmius=jrank-1
          !
          if(jrplus>max_rank)  jrplus=0
          if(jrmius<0)         jrmius=max_rank
          !
          fh(jrank,-1)=fh(jrmius,num_node-1)
          fh(jrank,num_node)=fh(jrplus,0)
          !
        enddo
        !
      enddo
      !
      !
      do jrank=0,max_rank
        do i=0,num_node
          df(jrank,i)=(fh(jrank,i)-fh(jrank,i-1))/(x(0,1)-x(0,0))
        enddo
        !
      enddo
      !
      error=0.d0
      do jrank=0,max_rank
        do i=1,num_node
          error=error+(df(jrank,i)-df_exact(jrank,i))**2
        enddo
      enddo
      error=sqrt(error/dble(tol_node))
      !
      write(18,*)num_rank,error
      !
      open(16,file='f_'//trim(scheme)//'.dat')
      do jrank=0,max_rank
      do i=0,num_node
        ! write(16,*)x(jrank,i),fh(jrank,i),f(jrank,i)
        write(16,*)x(jrank,i),df(jrank,i),df_exact(jrank,i)
        ! write(16,*)i,f(jrank,i)
      enddo
      enddo
      close(16)
      print*,' << f_',trim(scheme),'_parallel.dat'
      !
      deallocate(x,f,fh,df,df_exact)
      !
    enddo
    close(18)
    print*,' << error_'//trim(scheme)//'_rk_'//rkname//'.dat'
    !
  end subroutine errana1d_multi_domain
  !
  subroutine errana1d
    !
    integer :: n
    real(8),allocatable :: x(:),f(:),fh(:),df(:),df_exact(:)
    !
    integer :: i
    real(8) :: error
    character(len=8) :: scheme
    !
    scheme='cu5'
    !
    if(trim(scheme)=='cu5') then
      print*,' ** reconstruction scheme: 5th-order compact upwind'
    elseif(trim(scheme)=='cu5-ld') then
      print*,' ** reconstruction scheme: low-dissipative 5th-order compact upwind'
    else
      print*,scheme
      stop ' !! scheme not defined @ errana1d !!'
    endif
    !
    open(18,file='error_'//trim(scheme)//'.dat')
    n=8
    !
    do while(n<=512)
      !
      n=n*2
      !
      allocate(x(0:n),f(-4:n+4),fh(-1:n),df(0:n),df_exact(0:n))
      !
      do i=0,n
        x(i)=2.d0*pi/dble(n)*dble(i)
        f(i)=sin(4.d0*x(i))
        df_exact(i)=4.d0*cos(4.d0*x(i))
      enddo
      f(-4:-1)=f(n-4:n-1)
      f(n+1:n+4)=f(1:4)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! reconstruct interface value
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call recons(f,fh,n,scheme=trim(scheme),mode='sequential')
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! end of reconstruction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      do i=0,n
        df(i)=(fh(i)-fh(i-1))/(x(1)-x(0))
      enddo
      ! df=diffcc_recon(f,n,3,500)/(x(1)-x(0))
      !
      error=0.d0
      do i=0,n-1
        error=error+(df(i)-df_exact(i))**2
        !
        ! print*,i,df(i)-df_exact(i)
        !
      enddo
      error=sqrt(error/dble(n))
      !
      open(16,file='f_'//trim(scheme)//'.dat')
      do i=0,n
        write(16,*)x(i),df(i),df_exact(i)
      enddo
      close(16)
      print*,' << f_',trim(scheme),'.dat'
      !
      write(18,*)n,error
      !
      deallocate(x,f,fh,df,df_exact)
      !
    enddo
    !
    close(18)
    print*,' << error_',trim(scheme),'.dat'
    !
  end subroutine errana1d
  !
  subroutine tensoropt
    !
    integer :: i,j,k
    !
    do j=1,3
    do i=1,3
      write(*,'(2(A,I0,I0),A)',advance='no')'a',i,j,'*a',j,i,'+'
    enddo
    enddo
    !
    write(*,*)''
    write(*,*)'--------------------------------------------------------------'
    !
    do k=1,3
    do j=1,3
    do i=1,3
      write(*,'(3(A,I0,I0),A)',advance='no')'a',i,j,'*a',j,k,'*a',k,i,'+'
    enddo
    write(*,*)
    enddo
    enddo
    !
    write(*,*)
    !
  end subroutine tensoropt
  !
  subroutine WaveAnays
    !
    im=128
    !
    allocate(WaveK(0:im),WaveMK(0:im))
    !
    do i=0,im
      WaveK(i)=pi/im*i*1.d0
    end do
    !
    ! call RealSchemeAnay
    ! call schemederive_explicit(6)
    call schemederive_compact(6)
    ! call lowdiss_compactscheme(7,-3,3)
    ! call WA_7CC_LD
    ! call compact_scheme_analysis
    !
  end subroutine WaveAnays
  !
  !+--------------------------------------------------------------------+
  !| this subroutine is to analyse the spectral property using numerical|
  !| test                                                               | 
  !+--------------------------------------------------------------------+ 
  !| CHANGE RECORD                                                      |
  !| -------------                                                      |
  !| 10-02-2022: Created by J. Fang @ Warrington                        |
  !+--------------------------------------------------------------------+ 
  subroutine RealSchemeAnay
    !
    integer,parameter :: ima=128,nmax=64
    real(8) :: dx
    real(8) :: wavenumber(1:nmax)
    real(8) :: x(0:ima),rdf(0:ima),cdf(0:ima),ref(-hm:ima+hm),imf(-hm:ima+hm)
    complex(8) :: f(-hm:ima+hm),df(0:ima)
    complex(8) :: modified_wavenumber(1:nmax)
    !
    integer:: i,n
    !
    open(16,file='realtest.dat')
    do n=1,nmax
      !
      wavenumber(n)=dble(n)
      !
      do i=0,ima
        x(i)=2.d0*pi/ima*i
        f(i)=exp(wavenumber(n)*x(i)*cunit)
      enddo
      dx=x(1)
      !
      wavenumber(n)=wavenumber(n)
      !
      ! b.c.
      f(-4:-1)      =f(ima-4:ima-1)
      f(ima+1:ima+4)=f(1:4)
      !
      ref=real(f,8)
      imf=aimag(f)
      !
      ! rdf=diffcc_recon(ref,ima,3,5000) 
      ! cdf=diffcc_recon(imf,ima,3,5000) 
      rdf=diff5ec(ref,ima)
      cdf=diff5ec(imf,ima)
      !
      df=rdf*(1.d0,0.d0)+cdf*cunit
      df=df/dx
      !
      modified_wavenumber(n)=0.d0
      do i=1,ima
        modified_wavenumber(n)=modified_wavenumber(n)-cunit*df(i)/f(i)
      enddo
      modified_wavenumber(n)=modified_wavenumber(n)/dble(ima)
      modified_wavenumber(n)=modified_wavenumber(n)
      !
      print*,'modified wavenumber=',wavenumber(n)*dx,modified_wavenumber(n)*dx
      write(16,"(4(1X,E15.7E3))")wavenumber(n)*dx,real(modified_wavenumber(n))*dx,imag(modified_wavenumber(n))*dx, &
                                abs(real(modified_wavenumber(n))-wavenumber(n)**2)
      !
      open(12,file='test.dat')
      write(12,"(5(E15.7E3,2x))")(x(i),f(i),df(i),i=0,ima)
      close(12)
      print*,' << test.dat'
      !
    enddo
    close(16)
    print*,' << realtest.dat'
    !
  end subroutine RealSchemeAnay
  !+--------------------------------------------------------------------+
  !| the end of the subroutine   subroutine                             | 
  !+--------------------------------------------------------------------+ 
  !
  subroutine compact_scheme_analysis
    !
    integer :: nw,i,j,k
    real(8),allocatable :: weight(:,:,:),disint(:,:),todisint(:,:)
    real(8),allocatable :: bwreso(:,:,:),bwdisp(:,:,:)
    real(8) :: dissp_ref(0:im),dismax,disintmin
    real(8) :: w(2),spect_r(0:im),spect_i(0:im)
    integer :: cw1(2),cw2(2)
    !
    nw=500
    !
    allocate(weight(0:nw,0:nw,1:2))
    allocate(bwreso(0:nw,0:nw,0:im),bwdisp(0:nw,0:nw,0:im))
    allocate(disint(0:nw,0:nw),todisint(0:nw,0:nw))
    !
    do k=0,im
      dissp_ref(k)=-0.5d0*exx(10.d0*(WaveK(k)-0.9d0*pi))-0.5d0
    enddo
    open(18,file='exx.dat')
    do k=0,im
      write(18,*)WaveK(k),dissp_ref(k)
    enddo
    close(18)
    print*,' << exx.dat'
    !
    disint=1.d5
    todisint=1.d5
    disintmin=1.d10
    do j=0,nw
    do i=0,nw
      !
      weight(i,j,1)= 0.5d0+2.d0/nw*i-0.5d0*2.d0
      weight(i,j,2)= 0.5d0+2.d0/nw*j-0.5d0*2.d0
      ! print*,i,j,weight(i,j,1),weight(i,j,2)
      call schemeoptimise_compact(5,-3,2,weight(i,j,:),bwreso(i,j,:),bwdisp(i,j,:))
      !
      do k=0,im
        if(bwdisp(i,j,k)>0.d0 .or. bwdisp(i,j,36)<-1.d-3) then
        ! if(bwreso(i,j,k)>WaveK(k) .or. bwdisp(i,j,k)>0.d0) then
           goto 100
        endif
      enddo
      !
      dismax=bwdisp(i,j,im)
      !
      disint(i,j)=0.d0
      do k=1,im
        if(abs(bwdisp(i,j,k))>abs(dissp_ref(k))) then
          disint(i,j)=disint(i,j)+(dissp_ref(k)-bwdisp(i,j,k))**2
        endif
        todisint(i,j)=todisint(i,j)+dissp_ref(k)
      enddo
      !
      if(disint(i,j)<disintmin) then
        disintmin=disint(i,j)
        !
        cw1=fracnum(weight(i,j,1))
        cw2=fracnum(weight(i,j,2))
        !
        write(*,'(2(A,I0))',advance='no')' ** weight at i-3: ',cw1(1),'/',cw1(2)
        write(*,'(2(A,I0))',advance='no')', i+2: ',cw2(1),'/',cw2(2)
        print*,' dissipation: ',disintmin
        !
        ! print*,weight(i,j,1),weight(i,j,2),disintmin
        !
        w(1)=weight(i,j,1)
        w(2)=weight(i,j,2)
      endif
      !
      ! print*,i,j,dismax
      !
      100 continue
      !
    enddo
    enddo
    !
    !
    open(18,file='tecspectra.dat')
    write(18,'(A)')'TITLE = "spectra"'
    write(18,'(A)')'VARIABLES = "C-3", "C2", "dissp" ,"abs_dissp"'
    write(18,'(3(A,I0),A)')'ZONE T="BIG ZONE", I=',nw+1,', J=',nw+1,', K=',1,', DATAPACKING=POINT'
    do j=0,nw
    do i=0,nw
      write(18,'(4(1X,E15.7E3))')weight(i,j,1),weight(i,j,2),disint(i,j),todisint(i,j)
    enddo
    enddo
    close(18)
    print*,' << tecspectra.dat'
    !
    call schemeoptimise_compact(5,-3,2,w(:),spect_r,spect_i)
    !
    open(18,file='spectra.dat')
    write(18,"(3(1x,A18))")'K','MKreal','MKimag'
    write(18,"(3(1x,E18.11))")(WaveK(i),spect_r(i),spect_i(i),i=0,im)
    ! write(18,"(1x,A18)")'dissp'
    ! write(18,*)dissip
    close(18)
    print*,' << Out the spectra.dat...done!'
  end subroutine compact_scheme_analysis
  !
  subroutine WA_5UpWdE
    !
    do i=0,im
      WaveMK(i)=(0.d0,-1.d0)* (                                        &
             -2d0/60.d0*exp(-3d0*cunit*WaveK(i))              +&
             15d0/60.d0*exp(-2d0*cunit*WaveK(i))               &
                    - 1*exp(-1d0*cunit*WaveK(i))              +&
             20d0/60.d0*exp( 0d0*cunit*WaveK(i))              +&
             30d0/60.d0*exp( 1d0*cunit*WaveK(i))               &              
             -3d0/60.d0*exp( 2d0*cunit*WaveK(i))               )
    end do
    !
    open(18,file='5OUpWE.dat')
    write(18,"(3(1x,A18))")'K','MKreal','MKimag'
    write(18,"(3(1x,E18.11))")(WaveK(i),real(WaveMK(i)),               &
                                                aimag(WaveMK(i)),i=0,im)
    dissp=0.d0
    do i=1,im
      dissp=dissp+(WaveK(i)-WaveK(i-1))*aimag(WaveMK(i))
    end do
    !
    write(18,"(3(1x,A18))")'dissp'
    write(18,*)dissp
    !
    close(18)
    print*,' << Out the 5OUpWE.dat...done!'
    !
  end subroutine WA_5UpWdE
  !
  subroutine WA_5UpWdE_LD
    !
    integer :: n
    real(8) :: am3,am2,am1,a0,a1,a2,a3,dissp(0:10),a(-2:3)
    complex(8) :: WaveMK(0:10,0:im)
    !
    do n=0,10
      !
      a(3)=dble(n)/600.d0
      !
      a(-2)=-1.d0*a(3)+ 1.d0/30.d0
      a(-1)= 5.d0*a(3)-13.d0/60.d0
      a(0)=-10.d0*a(3)+47.d0/60.d0
      a(1)= 10.d0*a(3)+ 9.d0/20.d0
      a(2)= -5.d0*a(3)- 1.d0/20.d0
      !
      am3=-a(-2)
      am2=a(-2)-a(-1)
      am1=a(-1)-a(0)
      a0 = a(0)-a(1)
      a1 = a(1)-a(2)
      a2 = a(2)-a(3)
      a3 = a(3)
      !
      print*,am3,am2,am1,a0,a1,a2,a3
      !
      do i=0,im
        WaveMK(n,i)=(0.d0,-1.d0)* (                                      &
                      am3*exp(-3d0*cunit*WaveK(i))              +&
                      am2*exp(-2d0*cunit*WaveK(i))              +&
                      am1*exp(-1d0*cunit*WaveK(i))              +&
                       a0*exp( 0d0*cunit*WaveK(i))              +&
                       a1*exp( 1d0*cunit*WaveK(i))              +&              
                       a2*exp( 2d0*cunit*WaveK(i))              +&              
                       a3*exp( 3d0*cunit*WaveK(i))              )
      end do
      !
      dissp(n)=0.d0
      do i=1,im
        dissp(n)=dissp(n)+(WaveK(i)-WaveK(i-1))*aimag(WaveMK(n,i))
      end do
      !
    enddo
    !
    open(18,file='5OUpWE_LD.dat')
    do i=0,im
      write(18,"(23(1x,E18.11))")WaveK(i),(real(WaveMK(n,i)),n=0,10),(aimag(WaveMK(n,i)),n=0,10)
    enddo
    write(18,"(A)")'dissp'
    write(18,"(11(1x,E18.11))")(dissp(n),n=0,10)
    !
    close(18)
    print*,' << Out the 5OUpWE_LD.dat...done!'
    !
  end subroutine WA_5UpWdE_LD
  !
  subroutine WA_5CC_LD
    !
    integer :: n
    real(8) :: am3,am2,am1,a0,a1,a2,a3,dissp(0:10),a(-1:2),b(-2:2),alfa(-1:1)
    complex(8) :: WaveMK(0:10,0:im)
    !
    do n=0,10
      !
      a(2)=1.d0/360.d0*dble(n)
      !
      a(-1)=1.d0/18.d0-a(2)
      a(0)=19.d0/18.d0-a(2)*9.d0
      a(1)=10.d0/18.d0+a(2)*9.d0
      !
      !
      alfa(-1)=0.5d0-6.d0*a(2)
      alfa(1) =1.d0/6.d0+6.d0*a(2)
      !
      b(-2)=-a(-1)
      b(-1)= a(-1)-a(0)
      b(0) = a(0) -a(1)
      b(1) = a(1) -a(2)
      b(2) = a(2)
      !
      print*,b(:)
      !
      do i=0,im
        WaveMK(n,i)=(0.d0,-1.d0)* (                                      &
                      b(-2)*exp(-2d0*cunit*WaveK(i))              +&
                      b(-1)*exp(-1d0*cunit*WaveK(i))              +&
                       b(0)*exp( 0d0*cunit*WaveK(i))              +&
                       b(1)*exp( 1d0*cunit*WaveK(i))              +&              
                       b(2)*exp( 2d0*cunit*WaveK(i))              )
        WaveMK(n,i)=WaveMK(n,i)/(1.d0+alfa(-1)*exp(-1.d0*cunit*WaveK(i))+  &
                                      alfa(1) *exp( 1.d0*cunit*WaveK(i)))
      end do
      !
      dissp(n)=0.d0
      do i=1,im
        dissp(n)=dissp(n)+(WaveK(i)-WaveK(i-1))*aimag(WaveMK(n,i))
      end do
      !
    enddo
    !
    open(18,file='5Occ_LD.dat')
    do i=0,im
      write(18,"(23(1x,E18.11))")WaveK(i),(real(WaveMK(n,i)),n=0,10),(aimag(WaveMK(n,i)),n=0,10)
    enddo
    write(18,"(A)")'dissp'
    write(18,"(11(1x,E18.11))")(dissp(n),n=0,10)
    !
    close(18)
    print*,' << Out the 5Occ_LD.dat...done!'
    !
  end subroutine WA_5CC_LD
  !
  subroutine WA_5OCC_LD
    !
    integer :: n
    real(8) :: am3,am2,am1,a0,a1,a2,a3,dissp,a(-2:2),b(-3:2),alfa(-1:1)
    complex(8) :: WaveMK(0:im)
    !
    !
    a(-2)=-0.00444553d0
    a(-1)= 0.08101861d0
    a(0) = 0.99428149d0
    a(1) = 0.66721233d0
    a(2) = 0.01751043d0
    !
    !
    alfa(-1)=0.50163016d0
    alfa(1) =0.25394716d0
    !
    b(-3)=-a(-2)
    b(-2)= a(-2)-a(-1)
    b(-1)= a(-1)-a(0)
    b(0) = a(0) -a(1)
    b(1) = a(1) -a(2)
    b(2) = a(2)
    !
    do i=0,im
      WaveMK(i)=(0.d0,-1.d0)* (                                        &
                    b(-3)*exp(-3d0*cunit*WaveK(i))              +&
                    b(-2)*exp(-2d0*cunit*WaveK(i))              +&
                    b(-1)*exp(-1d0*cunit*WaveK(i))              +&
                     b(0)*exp( 0d0*cunit*WaveK(i))              +&
                     b(1)*exp( 1d0*cunit*WaveK(i))              +&              
                     b(2)*exp( 2d0*cunit*WaveK(i))              )
      WaveMK(i)=WaveMK(i)/(1.d0+alfa(-1)*exp(-1.d0*cunit*WaveK(i))+  &
                                alfa(1) *exp( 1.d0*cunit*WaveK(i)))
    end do
    !
    dissp=0.d0
    do i=1,im
      dissp=dissp+(WaveK(i)-WaveK(i-1))*aimag(WaveMK(i))
    end do
    !
    open(18,file='5occ.dat')
    write(18,"(3(1x,A18))")'K','MKreal','MKimag'
    write(18,"(3(1x,E18.11))")(WaveK(i),real(WaveMK(i)),               &
                                                aimag(WaveMK(i)),i=0,im)
    write(18,"(1x,A18)")'dissp'
    write(18,*)dissp
    close(18)
    print*,' << Out the 5occ.dat...done!'
    !
  end subroutine WA_5OCC_LD
  !
  subroutine WA_7CC_LD
    !
    integer :: n
    real(8) :: am3,am2,am1,a0,a1,a2,a3,dissp(0:10),a(-2:3),b(-3:3),alfa(-1:1)
    complex(8) :: WaveMK(0:10,0:im)
    !
    do n=0,10
      !
      a(3)=-1.d0/4800.d0*dble(n)
      !
      a(-2)= -1.d0/240.d0-a(3)
      a(-1)= 19.d0/240.d0+a(3)*15.d0
      a(0) =239.d0/240.d0+a(3)*80.d0
      a(1)=  53.d0/80.d0 -a(3)*80.d0
      a(2)=   1.d0/60.d0 -a(3)*15.d0
      !
      !
      alfa(-1)= 0.5d0+60.d0*a(3)
      alfa(1) =0.25d0-60.d0*a(3)
      !
      print*,n,'|',alfa(-1),alfa(1)
      print*,a
      !
      ! b(-2)=-a(-1)
      ! b(-1)= a(-1)-a(0)
      ! b(0) = a(0) -a(1)
      ! b(1) = a(1) -a(2)
      ! b(2) = a(2)
      !
      b(-3)=-a(-2)
      b(-2)= a(-2)-a(-1)
      b(-1)= a(-1)-a(0)
      b( 0)= a(0) -a(1)
      b( 1)= a(1) -a(2)
      b( 2)= a(2) -a(3)
      b( 3)= a(3)  
      !
      ! print*,b(:)
      !
      do i=0,im
        WaveMK(n,i)=(0.d0,-1.d0)* (                                      &
                      b(-3)*exp(-3d0*cunit*WaveK(i))              +&
                      b(-2)*exp(-2d0*cunit*WaveK(i))              +&
                      b(-1)*exp(-1d0*cunit*WaveK(i))              +&
                       b(0)*exp( 0d0*cunit*WaveK(i))              +&
                       b(1)*exp( 1d0*cunit*WaveK(i))              +&              
                       b(2)*exp( 2d0*cunit*WaveK(i))              +&              
                       b(3)*exp( 3d0*cunit*WaveK(i))              )
        WaveMK(n,i)=WaveMK(n,i)/(1.d0+alfa(-1)*exp(-1.d0*cunit*WaveK(i))+  &
                                      alfa(1) *exp( 1.d0*cunit*WaveK(i)))
      end do
      !
      dissp(n)=0.d0
      do i=1,im
        dissp(n)=dissp(n)+(WaveK(i)-WaveK(i-1))*aimag(WaveMK(n,i))
      end do
      !
    enddo
    !
    open(18,file='7Occ_LD.dat')
    do i=0,im
      write(18,"(23(1x,E18.11))")WaveK(i),(real(WaveMK(n,i)),n=0,10),(aimag(WaveMK(n,i)),n=0,10)
    enddo
    ! write(18,"(A)")'dissp'
    ! write(18,"(11(1x,E18.11))")(dissp(n),n=0,10)
    !
    close(18)
    print*,' << Out the 7Occ_LD.dat...done!'
    !
  end subroutine WA_7CC_LD
  !
  subroutine WA_5OCC
    !
    do i=0,im
      WaveMK(i)=(0.d0,-1.d0)* (                                        &
              -1.d0/3.d0*exp(-2.d0*cunit*WaveK(i))             &
             -18.d0/3.d0*exp(-1.d0*cunit*WaveK(i))             &
                   +3.d0*exp( 0.d0*cunit*WaveK(i))             &
             +10.d0/3.d0*exp( 1.d0*cunit*WaveK(i))         ) /(&
                    3.d0*exp(-1.d0*cunit*WaveK(i))+            & 
                    6.d0*exp(0.d0*cunit*WaveK(i))+             & 
                    1.d0*exp( 1.d0*cunit*WaveK(i))             )
    end do
    !
    open(18,file='5cc.dat')
    write(18,"(3(1x,A18))")'K','MKreal','MKimag'
    write(18,"(3(1x,E18.11))")(WaveK(i),real(WaveMK(i)),               &
                                                aimag(WaveMK(i)),i=0,im)
    close(18)
    print*,' << Out the 5cc.dat...done!'
    !
  end subroutine WA_5OCC
  !
  subroutine WA_6OCE
    !
    do i=0,im
      WaveMK(i)=(0.d0,-1.d0)* (                                        &
                 -1.d0/60.d0*exp(-3.d0*cunit*WaveK(i))         &
                 +9.d0/60.d0*exp(-2.d0*cunit*WaveK(i))         &
                -45.d0/60.d0*exp(-1.d0*cunit*WaveK(i))         &  
                +45.d0/60.d0*exp( 1.d0*cunit*WaveK(i))         &
                 -9.d0/60.d0*exp( 2.d0*cunit*WaveK(i))         &
                 +1.d0/60.d0*exp( 3.d0*cunit*WaveK(i))         )
                  
    end do
    !
    open(18,file='6ce.dat')
    write(18,"(3(1x,A18))")'K','MKreal','MKimag'
    write(18,"(3(1x,E18.11))")(WaveK(i),real(WaveMK(i)),               &
                                                aimag(WaveMK(i)),i=0,im)
    close(18)
    print*,' << Out the 6ce.dat...done!'
    !
  end subroutine WA_6OCE
  !
  subroutine WA_6OCC
    !
    do i=0,im
      WaveMK(i)=(0.d0,-1.d0)* (                                        &
                 -1d0/36d0*exp(-2.d0*cunit*WaveK(i))+          &
                 -7d0/9d0*exp(-1.d0*cunit*WaveK(i))+           &  
                  7d0/9d0*exp( 1.d0*cunit*WaveK(i))+           & 
                  1d0/36d0*exp(2.d0*cunit*WaveK(i))       ) /(&
                  1d0/3d0*exp(-1.d0*cunit*WaveK(i))+           & 
                    1.d0*exp(0.d0*cunit*WaveK(i))+             & 
                  1d0/3d0*exp( 1.d0*cunit*WaveK(i))           )
    end do
    !
    open(18,file='6cc.dat')
    write(18,"(3(1x,A18))")'K','MKreal','MKimag'
    write(18,"(3(1x,E18.11))")(WaveK(i),real(WaveMK(i)),               &
                                                aimag(WaveMK(i)),i=0,im)
    close(18)
    print*,' << Out the 6cc.dat...done!'
    !
  end subroutine WA_6OCC
  !
  subroutine spectranalysi(nstart,nend,a,alfa,rek,imk,totd)
    !
    integer,intent(in) :: nstart,nend
    real(8),intent(in) :: a(nstart:nend)
    real(8),intent(in),optional :: alfa(1:2)
    !
    real(8),intent(out),optional :: rek(0:im),imk(0:im),totd
    complex(8) :: varc
    !
    integer :: i,n
    !
    WaveMK=0.d0
    !
    do i=0,im
      !
      do n=nstart,nend
        !
        WaveMK(i)=WaveMK(i)+a(n)*exp(dble(n)*cunit*WaveK(i))
        !
      enddo
      WaveMK(i)=WaveMK(i)*(0.d0,-1.d0)
      !
    end do
    !
    if(present(alfa)) then
      do i=0,im
        !
        varc=1.d0
        do n=1,2
          !
          varc=varc+alfa(n)*exp(dble(n*2-3)*cunit*WaveK(i))
          !
        enddo
        WaveMK(i)=WaveMK(i)/varc
        !
      end do
    endif
    !
    do i=0,im
      if(present(rek)) rek(i)=real(WaveMK(i))
      if(present(imk)) imk(i)=aimag(WaveMK(i))
    enddo
    !
    if(present(totd)) then
      totd=0.d0
      do i=1,im
        totd=totd+(WaveK(i)-WaveK(i-1))*aimag(WaveMK(i))
      end do
    endif
    !
    return
    !
  end subroutine spectranalysi
  !
  subroutine laptest
    !
    ! implicit none (type, external)
    ! external :: sgesv
    ! real     :: a(2, 2)  ! Matrix A.
    ! real     :: b(2)     ! Vector b/x.
    ! real     :: pivot(2) ! Pivot indices (list of swap operations).
    ! integer  :: rc       ! Return code.

    ! ! a = reshape([ 2., 3., 1., 1. ], [ 2, 2 ])
    ! a(1,1)=2.d0
    ! a(1,2)=1.d0
    ! a(2,1)=3.d0
    ! a(2,2)=1.d0
    ! b = [ 5., 6. ]
    ! print*, a
    ! print*, b

    ! call sgesv(2, 1, a, 2, pivot, b, 2, rc)


    ! if (rc /= 0) then
    !     print '(a, i0)', 'Error: ', rc
    !     stop
    ! end if

    ! print '("Solution (x1, x2): ", f0.4, ", ", f0.4)', b
  end subroutine laptest
  !
  function fracnum(r) result(f)
    !
    real(8),intent(in) :: r
    integer :: f(2)
    !
    integer :: i,ifst,ints
    real(8) :: resd,ri,resdmin
    !
    ifst=int(1.d0/r-1)
    ifst=max(ifst,1)
    i=ifst
    resd=10.d0
    resdmin=10.d0
    do while(.true.)
      !
      ri=r*dble(i)
      ints=int(ri)
      resd=abs(dble(ints)-ri)
      !
      if(resd<1.d-8 .or. abs(resd-1.d0)<1.d-8) exit
      !
      i=i+1
      !
      resdmin=min(resdmin,resd)
      !
    enddo
    !
    f(1)=nint(ri)
    f(2)=i
    !
    return
    !
  end  function fracnum
  !
  subroutine schemeoptimise_compact(norder,nstart,nend,weight,spect_r,spect_i)
    !
    integer,intent(in) ::  norder,nstart,nend
    real(8),intent(in) ::  weight(:)
    real(8),intent(out) :: spect_r(0:im),spect_i(0:im)
    !
    integer :: n,j,ns,ne,no
    real(8) :: d0
    integer,allocatable :: next(:)
    real(8),allocatable :: a(:,:),ainv(:,:),ao(:,:),b(:),c(:),d(:),hh(:)
    integer,allocatable :: rhs(:,:),lhs(:,:)
    !
    real(8) :: dissip,var1
    !
    ns=-(norder+1)/2+1
    ne=ns+norder-2
    !
    ! print*,' ** the normal stencile:',ns,ne
    ! !
    ! print*,' ** the extended stencile:',nstart,nend
    !
    no=ns-nstart+nend-ne
    ! print*,' ** number of extra node:',no
    !
    allocate(a(1:norder,1:norder),ainv(1:norder,1:norder),c(1:norder),d(1:norder),b(nstart:nend))
    allocate(ao(1:norder,1:no),next(1:no))
    allocate(rhs(nstart:nend,1:2),lhs(-1:1,1:2))
    !
    j=0
    do n=ns,ne
      if(n==0) cycle
      j=j+1
      a(:,j)=TaylorSeriesCoef(n,norder,'f')
    enddo
    !
    do n=-1,1,2
      j=j+1
      a(:,j)=TaylorSeriesCoef(n,norder,'df')
    enddo
    ainv=matinv(a)
    !
    c=0
    c(1)=1.d0
    !
    ! for extra nodes
    j=0
    do n=nstart,ns-1
      j=j+1
      ao(:,j)=TaylorSeriesCoef(n,norder,'f')
    enddo 
    do n=ne+1,nend
      j=j+1
      ao(:,j)=TaylorSeriesCoef(n,norder,'f')
    enddo
    !
    ! weight=0.5d0/36.d0
    !
    do n=1,no
      do j=1,norder
        c(j)=c(j)-weight(n)*ao(j,n)
      enddo
    enddo
    !
    d=matmul(ainv,c)
    d0=-sum(d(1:norder-2))
    !
    ! print*,' ** first derive'
    ! var1=0
    ! do n=1,norder
    !   var1=var1+a(1,n)*d(n)
    !   print*,n,'|',d(n),var1+weight*ao(1,1)
    ! enddo
    ! print*,' ** second derive'
    ! var1=0
    ! do n=1,norder
    !   var1=var1+a(2,n)*d(n)
    !   print*,n,'|',d(n),var1+weight*ao(2,1)
    ! enddo
    !
    j=0
    do n=ns,ne
      if(n==0) cycle
      j=j+1

      b(n)=d(j)
    enddo
    b(0)=d0
    !
    j=0
    !
    do n=nstart,ns-1
      j=j+1
      b(n)=weight(j)
    enddo 
    do n=ne+1,nend
      j=j+1
      b(n)=weight(j)
    enddo
    !
    ! do n=nstart,nend
    !   rhs(n,:)=fracnum(b(n))
    ! enddo
    ! !
    ! lhs(-1,:)=fracnum(-d(norder-1))
    ! lhs( 1,:)=fracnum(-d(norder))
    ! !
    ! write(*,'(A,I0,A)')'  ** The optimised ',norder,'th-order compact scheme:'
    ! write(*,'(A,4(I0,A))')'  ** ',lhs(-1,1),'/',lhs(-1,2),'*df(i-1) + df(i) + ', &
    !                               lhs(1,1),'/',lhs(1,2),'*df(i+1) = '
    ! do n= nstart,nend
    !   do j=0,n-nstart
    !     write(*,'(10X)',advance='no')
    !   enddo
    !   if(n==0)  write(*,'(I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i) + '
    !   if(n>0)   write(*,'(I0,A,I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i+',n,') + '
    !   if(n<0)   write(*,'(I0,A,I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i',n,') + '
    ! enddo
    ! write(*,*)
    !
    call spectranalysi(nstart=nstart,nend=nend,a=b(nstart:nend), &
                       alfa=-d(norder-1:norder), &
                       rek=spect_r,imk=spect_i,totd=dissip)
    !
    ! open(18,file='spectra.dat')
    ! write(18,"(3(1x,A18))")'K','MKreal','MKimag'
    ! write(18,"(3(1x,E18.11))")(WaveK(i),spect_r(i),spect_i(i),i=0,im)
    ! ! write(18,"(1x,A18)")'dissp'
    ! ! write(18,*)dissip
    ! close(18)
    ! print*,' << Out the spectra.dat...done!'
    !
    !
  end subroutine schemeoptimise_compact
  !
  subroutine lowdiss_compactscheme(norder,nstart,nend)
    !
    integer,intent(in) ::  norder,nstart,nend
    !
    integer :: n,j,ns,ne,no
    real(8) :: d0
    integer,allocatable :: next(:)
    real(8),allocatable :: a(:,:),ainv(:,:),ao(:,:),b(:),c(:),d(:),hh(:)
    integer,allocatable :: rhs(:,:),rhs2(:,:),lhs(:,:),lhs2(:,:)
    real(8),allocatable :: alfa(:),bb(:)
    real(8) :: spect_r(0:im),spect_i(0:im)
    real(8) ::  weight
    !
    real(8) :: dissip,var1
    character(len=2) :: oname
    !
    ns=-(norder+1)/2+1
    ne=ns+norder-2
    !
    print*,' ** the normal stencile:',ns,ne
    !
    print*,' ** the extended stencile:',nstart,nend
    !
    no=ns-nstart+nend-ne
    print*,' ** number of extra node:',no
    !
    ! conventional upwind scheme
    allocate(alfa(1:2),bb(nstart+1:nend-1))
    allocate(rhs2(nstart:nend,1:2),lhs2(-1:1,1:2))
    !
    call schemederive_compact(norder,alfa=alfa,bb=bb)
    !
    rhs2=0
    lhs2=0
    ! do n=nstart+1,nend-1
    !   rhs2(n,:)=fracnum(bb(n))
    !   print*,n,'|',rhs2(n,:)
    ! enddo
    !
    allocate(a(1:norder,1:norder),ainv(1:norder,1:norder),c(1:norder),d(1:norder),b(nstart:nend))
    allocate(ao(1:norder,1:no),next(1:no))
    allocate(rhs(nstart:nend,1:2),lhs(-1:1,1:2))
    !
    j=0
    do n=ns,ne
      if(n==0) cycle
      j=j+1
      a(:,j)=TaylorSeriesCoef(n,norder,'f')
    enddo
    !
    do n=-1,1,2
      j=j+1
      a(:,j)=TaylorSeriesCoef(n,norder,'df')
    enddo
    ainv=matinv(a)
    !
    c=0
    c(1)=1.d0
    !
    ! for extra nodes
    j=0
    do n=nstart,ns-1
      j=j+1
      ao(:,j)=TaylorSeriesCoef(n,norder,'f')
    enddo 
    do n=ne+1,nend
      j=j+1
      ao(:,j)=TaylorSeriesCoef(n,norder,'f')
    enddo
    !
    weight=-1.d0/480.d0
    !
    do j=1,norder
      c(j)=c(j)-weight*ao(j,1)
    enddo
    !
    d = matmul(ainv,c)
    d0=-sum(d(1:norder-2))
    !
    j=0
    do n=ns,ne
      if(n==0) cycle
      j=j+1

      b(n)=d(j)
    enddo
    b(0)=d0
    !
    j=0
    !
    do n=nstart,ns-1
      j=j+1
      b(n)=weight
    enddo 
    do n=ne+1,nend
      j=j+1
      b(n)=weight
    enddo
    !
    do n=nstart,nend
      rhs(n,:)=fracnum(b(n))
    enddo
    !
    ! lhs(-1,:)=fracnum((-d(norder-1)))
    ! lhs( 1,:)=fracnum((-d(norder)  ))
    !
    lhs(-1,:)=fracnum(alfa(1))
    lhs( 1,:)=fracnum(alfa(2))
    ! !
    lhs2(-1,:)=fracnum((-d(norder-1)-alfa(1))/weight)
    lhs2( 1,:)=fracnum((-d(norder)-alfa(2))/weight)
    !
    write(*,'(A,I0,A)')'  ** The optimised ',norder,'th-order low-dissipation compact scheme:'
    write(*,'(A,8(I0,A))')'  ** (',lhs(-1,1),'/',lhs(-1,2),'+',lhs2(-1,1),'/',lhs2(-1,2),'*b)df(i-1) + df(i) + (', &
                                  lhs(1,1),'/',lhs(1,2),'+',lhs2( 1,1),'/',lhs2( 1,2),'*b)df(i+1) = '
    do n= nstart,nend
      do j=0,n-nstart
        write(*,'(10X)',advance='no')
      enddo
      if(n==0)  write(*,'(I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i) + '
      if(n>0)   write(*,'(I0,A,I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i+',n,') + '
      if(n<0)   write(*,'(I0,A,I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i',n,') + '
    enddo
    ! write(*,*)
    !
    deallocate(rhs)
    !
    ! stop
    !
    allocate(hh(nstart+1:nend),rhs(nstart+1:nend,1:2))
    ! bb(nstart+1:nend-1)
    ! interface scheme
    hh(nstart+1)=-b(nstart)
    do n=nstart+2,nend-1
      hh(n)=hh(n-1)-b(n-1)
    enddo
    hh(nend)=b(nend)
    !
    do n=nstart+1,nend
      rhs(n,:)=fracnum(hh(n))
    enddo
    !
    ! do n=nstart+1,nend-1
    !   rhs(n,:)=fracnum(bb(n))
    ! enddo
    ! rhs(nend,:)=fracnum(0.d0)
    ! !
    ! do n=nstart+1,nend-1
    !   rhs2(n,:)=fracnum((hh(n)-bb(n))/weight)
    ! enddo
    ! rhs2(nend,:)=fracnum( hh(n)/weight)
    !
    write(*,'(A,I0,A)')'  ** The optimised ',norder,'th-order low-dissipation compact scheme:'
    write(*,*)'  ** The value of the weight ',weight
    write(*,'(A,8(I0,A))')'  ** (',lhs(-1,1),'/',lhs(-1,2),'+',lhs2(-1,1),'/',lhs2(-1,2),'*b)*f(i-1/2) + df(i+1/2) + (', &
                                  lhs(1,1),'/',lhs(1,2),'+',lhs2(1,1),'/',lhs2(1,2),'*b)*f(i+3/2) = '
    do n= nstart+1,nend
      do j=1,n-nstart
        write(*,'(10X)',advance='no')
      enddo
      if(n==0)  write(*,'(I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i) + '
      if(n>0)   write(*,'(I0,A,I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i+',n,') + '
      if(n<0)   write(*,'(I0,A,I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i',n,') + '
      ! if(n==0)  write(*,'(4(A,I0),A)')'(',rhs(n,1),'/',rhs(n,2),'+',rhs2(n,1),'/',rhs2(n,2),'*b)*f(i) + '
      ! if(n>0)   write(*,'(5(A,I0),A)')'(',rhs(n,1),'/',rhs(n,2),'+',rhs2(n,1),'/',rhs2(n,2),'*b)*f(i+',n,') + '
      ! if(n<0)   write(*,'(5(A,I0),A)')'(',rhs(n,1),'/',rhs(n,2),'+',rhs2(n,1),'/',rhs2(n,2),'*b)*f(i',n,') + '
    enddo
    write(*,*)
    !
    call spectranalysi(nstart=nstart,nend=nend,a=b(nstart:nend), &
                       alfa=-d(norder-1:norder), &
                       rek=spect_r,imk=spect_i,totd=dissip)
    !
    write(oname,'(i2.2)')norder
    open(18,file='spectra_'//oname//'order-LD_compact.dat')
    write(18,"(3(1x,A18))")'K','MKreal','MKimag'
    write(18,"(3(1x,E18.11))")(WaveK(i),spect_r(i),spect_i(i),i=0,im)
    close(18)
    print*,' << Out the spectra_',oname,'order-LD_compact.dat...done!'
    !
    !
  end subroutine lowdiss_compactscheme
  !
  !+--------------------------------------------------------------------+
  !| this subroutine is to derive scheme of compact finite-difference   |
  !| scheme by a given the order of truncation error.                   | 
  !+--------------------------------------------------------------------+ 
  !| CHANGE RECORD                                                      |
  !| -------------                                                      |
  !| 05-02-2022: Created by J. Fang @ Warrington                        |
  !+--------------------------------------------------------------------+ 
  subroutine schemederive_compact(norder,alfa,bb)
    !
    real(8),intent(out),allocatable,optional :: alfa(:),bb(:)
    !
    integer,intent(in) :: norder
    integer :: nstart,nend,n,j
    real(8) :: d0
    real(8),allocatable :: a(:,:),b(:),c(:),d(:),hh(:)
    integer,allocatable :: rhs(:,:),lhs(:,:)
    real(8) :: bwreso(0:im),bwdisp(0:im)
    character(len=2) :: oname
    !
    nstart=-(norder+1)/2+1
    nend=nstart+norder-2
    !
    write(*,'(2(A,I0),A)')'  ** stencile is from: [i',nstart,' ~ i+',nend,']'
    !
    allocate(a(1:norder,1:norder),c(1:norder),d(1:norder),b(nstart:nend))
    allocate(rhs(nstart:nend,1:2),lhs(-1:1,1:2))
    !
    j=0
    do n=nstart,nend
      if(n==0) cycle
      j=j+1
      a(:,j)=TaylorSeriesCoef(n,norder,'f')
    enddo
    !
    do n=-1,1,2
      j=j+1
      a(:,j)=TaylorSeriesCoef(n,norder,'df')
    enddo
    !
    a=matinv(a)
    !
    c=0
    c(1)=1.d0
    !
    d=matmul(a,c)
    d0=-sum(d(1:norder-2))
    !
    j=0
    do n=nstart,nend
      if(n==0) cycle
      j=j+1

      rhs(n,:)=fracnum(d(j))
      b(n)=d(j)
    enddo
    rhs(0,:)=fracnum(d0)
    b(0)=d0
    !
    lhs(-1,:)=fracnum(-d(norder-1))
    lhs( 1,:)=fracnum(-d(norder))
    !
    write(*,'(A,I0,A)')'  ** The formula of the ',norder,'th-order compact scheme:'
    write(*,'(A,4(I0,A))')'     ',lhs(-1,1),'/',lhs(-1,2),'*df(i-1) + df(i) + ', &
                                  lhs(1,1),'/',lhs(1,2),'*df(i+1) = '
    do n= nstart,nend
      do j=0,n-nstart
        write(*,'(10X)',advance='no')
      enddo
      if(n==0)  write(*,'(I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i) + '
      if(n>0)   write(*,'(I0,A,I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i+',n,') + '
      if(n<0)   write(*,'(I0,A,I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i',n,') + '
    enddo
    write(*,*)
    !
    deallocate(rhs)
    !
    allocate(hh(nstart+1:nend),rhs(nstart+1:nend,1:2))
    !
    ! interface scheme
    hh(nstart+1)=-b(nstart)
    do n=nstart+2,nend-1
      hh(n)=hh(n-1)-b(n-1)
    enddo
    hh(nend)=b(nend)
    !
    do n=nstart+1,nend
      rhs(n,:)=fracnum(hh(n))
    enddo
    !
    if(present(alfa) .and. present(bb)) then
      allocate(alfa(1:2),bb(nstart+1:nend))
      !
      alfa(1)=-d(norder-1)
      alfa(2)=-d(norder)
      !
      bb(nstart+1:nend)=hh(nstart+1:nend)
      !
    endif
    !
    write(*,'(A,4(I0,A))')'     ',lhs(-1,1),'/',lhs(-1,2),'*f(i-1/2) + f(i+1/2) + ', &
                                  lhs(1,1),'/',lhs(1,2),'*f(i+3/2) = '
    do n= nstart+1,nend
      do j=1,n-nstart
        write(*,'(10X)',advance='no')
      enddo
      if(n==0)  write(*,'(I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i) + '
      if(n>0)   write(*,'(I0,A,I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i+',n,') + '
      if(n<0)   write(*,'(I0,A,I0,A,I0,A)')rhs(n,1),'/',rhs(n,2),'*f(i',n,') + '
    enddo
    write(*,*)
    !
    call spectranalysi(nstart=nstart,nend=nend,a=b(nstart:nend), &
                       alfa=-d(norder-1:norder), &
                       rek=bwreso,imk=bwdisp)
    !
    write(oname,'(i2.2)')norder
    open(18,file='spectra_'//oname//'order_compact.dat')
    write(18,"(3(1x,A18))")'K','MKreal','MKimag'
    write(18,"(3(1x,E18.11))")(WaveK(i),bwreso(i),bwdisp(i),i=0,im)
    close(18)
    print*,' << Out the spectra_',oname,'order_compact.dat...done!'
    !
    ! write(*,'(A)')' ** f(i+1/2) = '
    ! do n= nstart+1,nstart+norder
    !   do j=1,n-nstart
    !     write(*,'(10X)',advance='no')
    !   enddo
    !   if(n==0)  write(*,'(I0,A,I0,A)')he(n,1),'/',he(n,2),'*f(i) + '
    !   if(n>0)   write(*,'(I0,A,I0,A,I0,A)')he(n,1),'/',he(n,2),'*f(i+',n,') + '
    !   if(n<0)   write(*,'(I0,A,I0,A,I0,A)')he(n,1),'/',he(n,2),'*f(i',n,') + '
    ! enddo
    ! write(*,*)
    ! !
    ! do n=nstart,nend
    !   write(*,'(3(I0,A))')n,':',rhs(n,1),'/',rhs(n,2),''
    ! enddo
    ! do n=-1,1,2
    !   write(*,'(3(I0,A))')n,':',lhs(n,1),'/',lhs(n,2),''
    ! enddo
    !
  end subroutine schemederive_compact
  !+---------------------------------------------------------------------+
  !| The end of the subroutine schemederive.                             |
  !+---------------------------------------------------------------------+
  !
  !+--------------------------------------------------------------------+
  !| this subroutine is to derive scheme of explicit finite-difference  |
  !| scheme by a given the order of truncation error.                   | 
  !+--------------------------------------------------------------------+ 
  !| CHANGE RECORD                                                      |
  !| -------------                                                      |
  !| 05-02-2022: Created by J. Fang @ Warrington                        |
  !+--------------------------------------------------------------------+                                    
  subroutine schemederive_explicit(norder)
    !
    external :: sgesv
    !
    integer,intent(in) :: norder
    integer :: nstart,nend,n,j
    integer,allocatable :: e(:,:),he(:,:)
    real(8),allocatable :: a(:,:),b(:,:),c(:),d(:),dd(:),hh(:)
    real(8) :: bwreso(0:im),bwdisp(0:im)
    character(len=2) :: oname
    !
    nstart=-(norder+1)/2
    nend=nstart+norder
    !
    write(*,'(2(A,I0),A)')'  ** stencile is from: [i',nstart,' ~ i+',nend,']'
    !
    allocate(a(1:norder,nstart:nend),b(norder,norder),c(norder),d(norder), &
             dd(nstart:nend),e(nstart:nend,1:2),                  &
             hh(nstart+1:nend),he(nstart+1:nend,1:2))
    !
    do n=nstart,nend
      a(:,n)=TaylorSeriesCoef(n,norder,'f')
    enddo
    !
    b(:,        1:-nstart)=a(:,nstart:-1)
    b(:,-nstart+1:norder) =a(:,1:nend)
    !
    b=matinv(b)
    !
    c=0
    c(1)=1.d0
    !
    d=MATMUL(b,c)
    !
    ! difference scheme
    dd(nstart:-1)=d(1:-nstart)
    dd(1:nend)=d(-nstart+1:norder)
    dd(0)=-sum(d)
    !
    ! interface scheme
    hh(nstart+1)=-dd(nstart)
    do n=nstart+2,nend-1
      hh(n)=hh(n-1)-dd(n-1)
    enddo
    hh(nend)=dd(nend)
    !
    do n=nstart,nend
      e(n,:)=fracnum(dd(n))
    enddo
    do n=nstart+1,nend
      he(n,:)=fracnum(hh(n))
    enddo
    !
    write(*,'(A,I0,A)')'  ** The formula of the ',norder,'th-order scheme:'
    write(*,'(A)')'     df/dx(i) = '
    do n= nstart,nend
      do j=0,n-nstart
        write(*,'(10X)',advance='no')
      enddo
      if(n==0)  write(*,'(I0,A,I0,A)')e(n,1),'/',e(n,2),'*f(i) + '
      if(n>0)   write(*,'(I0,A,I0,A,I0,A)')e(n,1),'/',e(n,2),'*f(i+',n,') + '
      if(n<0)   write(*,'(I0,A,I0,A,I0,A)')e(n,1),'/',e(n,2),'*f(i',n,') + '
    enddo
    write(*,*)
    write(*,'(A)')'     f(i+1/2) = '
    do n= nstart+1,nend
      do j=1,n-nstart
        write(*,'(10X)',advance='no')
      enddo
      if(n==0)  write(*,'(I0,A,I0,A)')he(n,1),'/',he(n,2),'*f(i) + '
      if(n>0)   write(*,'(I0,A,I0,A,I0,A)')he(n,1),'/',he(n,2),'*f(i+',n,') + '
      if(n<0)   write(*,'(I0,A,I0,A,I0,A)')he(n,1),'/',he(n,2),'*f(i',n,') + '
    enddo
    write(*,*)
    !
    call spectranalysi(nstart=nstart,nend=nend,a=dd(nstart:nend), &
                       rek=bwreso,imk=bwdisp)
    !
    write(oname,'(i2.2)')norder
    open(18,file='spectra_'//oname//'order_explicit.dat')
    write(18,"(3(1x,A18))")'K','MKreal','MKimag'
    write(18,"(3(1x,E18.11))")(WaveK(i),bwreso(i),bwdisp(i),i=0,im)
    close(18)
    print*,' << Out the spectra_',oname,'order_explicit.dat...done!'
    !
  end subroutine schemederive_explicit
  !+---------------------------------------------------------------------+
  !| The end of the subroutine schemederive.                             |
  !+---------------------------------------------------------------------+
  !
  subroutine schemederive5
    !
    external :: sgesv
    !
    integer :: norder,nstart,n
    real(8),allocatable :: a(:,:),b(:,:),c(:),d(:)
    !
    norder=5
    nstart=-3
    allocate(a(1:norder,nstart:nstart+norder),b(norder,norder),c(norder),d(norder))
    !
    do n=nstart,nstart+norder
      a(:,n)=TaylorSeriesCoef(n,norder,'f')
    enddo
    !
    b(:,1)=a(:,-3)
    b(:,2)=a(:,-2)
    b(:,3)=a(:,-1)
    b(:,4)=a(:,1)
    b(:,5)=a(:,2)
    !
    b=matinv(b)
    !
    c=0
    c(1)=1.d0
    !
    d=MATMUL(b,c)
    !
    print*,d(1)
    print*,d(2)
    print*,d(3)
    print*,d(4)
    print*,d(5)
    !
    !
  end subroutine schemederive5
  !
  function matinv(A) result(Ainv)
    real(8), dimension(:,:), intent(in) :: A
    real(8), dimension(size(A,1),size(A,2)) :: Ainv

    real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! ! External procedures defined in LAPACK
    ! external DGETRF
    ! external DGETRI

    ! ! Store A in Ainv to prevent it from being overwritten by LAPACK
    ! Ainv = A
    ! n = size(A,1)

    ! ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! ! using partial pivoting with row interchanges.
    ! call DGETRF(n, n, Ainv, n, ipiv, info)

    ! if (info /= 0) then
    !    stop 'Matrix is numerically singular!'
    ! end if

    ! ! DGETRI computes the inverse of a matrix using the LU factorization
    ! ! computed by DGETRF.
    ! ! call DGETRI(n, Ainv, n, ipiv, work, n, info)

    ! if (info /= 0) then
    !    stop 'Matrix inversion failed!'
    ! end if
  end function matinv
  !
  function TaylorSeriesCoef(delta,n,char) result(c)
    !
    integer,intent(in) :: delta,n
    character(len=*),intent(in) ::  char
    real(8) :: c(1:n)
    !
    integer :: i
    !
    if(char=='f') then
      do i=1,n
        c(i)=dble(delta**i)/dfac(i)
      enddo
    elseif(char=='df') then
      c(1)=1.d0
      do i=2,n
        c(i)=dble(delta**(i-1))/dfac(i-1)
      enddo
    else
      stop ' !! error of char @ TaylorSeriesCoef'
    endif
    !
  end function TaylorSeriesCoef
  !
  function dfac(n)
    !
    real(8) :: dfac
    integer :: n
    !
    dfac=1.d0
    !
    do i=1,n
      dfac=dfac*i
    end do
    !
    return
    !
  end function dfac
  !
  real(8) function exx(x)
    !
    real(8),intent(in) :: x
    !
    exx=(exp(x)-exp(-x))/(exp(x)+exp(-x))
    !
  end function exx
  !
end module numerics
!