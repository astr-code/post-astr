!+---------------------------------------------------------------------+
!| This module contains subroutines to process specific flow case.     |
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2020-05-05                                     |
!+---------------------------------------------------------------------+
module userdefine
  !
  implicit none
  !
  contains
  !
  subroutine supplymodify
    !
    use h5readwrite
    !
    real(8) :: forcex,forcey,forcez,q0
    integer :: nabnorm,nfile,nfile2d,ninlet,nins,nnorm,npabnorm,       &
               nsamples,nsave,nsave1,nsave2,nsgs,nstep,nstepsave1,     &
               nstepsave2
    real(8) :: randomv(0:14),timesave1,timesave2,twall
    !
    call H5ReadArray(forcex,'forcex','Outdat/supply.h5')
    call H5ReadArray(forcey,'forcey','Outdat/supply.h5')
    call H5ReadArray(forcez,'forcez','Outdat/supply.h5')
    call H5ReadArray(q0,'q0','Outdat/supply.h5')
    call H5ReadArray(nabnorm,'nabnorm','Outdat/supply.h5')
    call H5ReadArray(nfile,'nfile','Outdat/supply.h5')
    call H5ReadArray(nfile2d,'nfile2d','Outdat/supply.h5')
    call H5ReadArray(ninlet,'ninlet','Outdat/supply.h5')
    call H5ReadArray(nins,'nins','Outdat/supply.h5')
    call H5ReadArray(nnorm,'nnorm','Outdat/supply.h5')
    call H5ReadArray(npabnorm,'npabnorm','Outdat/supply.h5')
    call H5ReadArray(nsamples,'nsamples','Outdat/supply.h5')
    call H5ReadArray(nsave,'nsave','Outdat/supply.h5')
    call H5ReadArray(nsave1,'nsave1','Outdat/supply.h5')
    call H5ReadArray(nsave2,'nsave2','Outdat/supply.h5')
    call H5ReadArray(nsgs,'nsgs','Outdat/supply.h5')
    call H5ReadArray(nstep,'nstep','Outdat/supply.h5')
    call H5ReadArray(nstepsave1,'nstepsave1','Outdat/supply.h5')
    call H5ReadArray(nstepsave2,'nstepsave2','Outdat/supply.h5')
    call H5ReadArray(timesave1,'timesave1','Outdat/supply.h5')
    call H5ReadArray(timesave2,'timesave2','Outdat/supply.h5')
    call H5ReadArray(twall,'twall','Outdat/supply.h5')
    call H5ReadArray(randomv,14,'randomv','Outdat/supply.h5')
    !
    ! nfile=139
    ! nfile2d=731
    ! nstep=150000
    ! nins=13190
    ! nsamples=1000
    nstep=45000
    !
    call H5WriteArray(forcex,'forcex','Results/supply.h5')
    call H5WriteArray(forcey,'forcey','Results/supply.h5')
    call H5WriteArray(forcez,'forcez','Results/supply.h5')
    call H5WriteArray(q0,'q0','Results/supply.h5')
    call H5WriteArray(nabnorm,'nabnorm','Results/supply.h5')
    call H5WriteArray(nfile,'nfile','Results/supply.h5')
    call H5WriteArray(nfile2d,'nfile2d','Results/supply.h5')
    call H5WriteArray(ninlet,'ninlet','Results/supply.h5')
    call H5WriteArray(nins,'nins','Results/supply.h5')
    call H5WriteArray(nnorm,'nnorm','Results/supply.h5')
    call H5WriteArray(npabnorm,'npabnorm','Results/supply.h5')
    call H5WriteArray(nsamples,'nsamples','Results/supply.h5')
    call H5WriteArray(nsave,'nsave','Results/supply.h5')
    call H5WriteArray(nsave1,'nsave1','Results/supply.h5')
    call H5WriteArray(nsave2,'nsave2','Results/supply.h5')
    call H5WriteArray(nsgs,'nsgs','Results/supply.h5')
    call H5WriteArray(nstep,'nstep','Results/supply.h5')
    call H5WriteArray(nstepsave1,'nstepsave1','Results/supply.h5')
    call H5WriteArray(nstepsave2,'nstepsave2','Results/supply.h5')
    call H5WriteArray(timesave1,'timesave1','Results/supply.h5')
    call H5WriteArray(timesave2,'timesave2','Results/supply.h5')
    call H5WriteArray(twall,'twall','Results/supply.h5')
    call H5WriteArray(randomv,14,'randomv','Results/supply.h5')
    !
    !
  end subroutine supplymodify
  !
  subroutine tgvstate(infile)
    !
    character(len=*),intent(in) :: infile
    !
    integer :: nlimt,nmax,nvar
    parameter (nlimt=10000000,nvar=14)
    integer :: nstep
    real(8) :: time(0:nlimt),var(0:nlimt,1:nvar),dekdt(0:nlimt)
    real(8) :: deltat
    integer :: io,j,n
    !
    print*,infile
    !
    open(12,file=infile)
    read(12,*)
    io=0
    do while(io==0)
      read(12,*,iostat=io)j,time(j),(var(j,n),n=1,nvar)
    enddo
    close(12)
    print*, ' >> ',infile
    !
    nmax=j
    !
    deltat=time(1)
    dekdt(0)=-(var(1,1)-var(0,1))/deltat
    do j=1,nmax-1
      dekdt(j)=-0.5d0*(var(j+1,1)-var(j-1,1))/deltat
    enddo
    dekdt(nmax)=-(var(nmax,1)-var(nmax-1,1))/deltat
    !
    open(18,file=infile)
    read(18,*)
    do j=0,nmax,1
      write(18,"(I7,1X,E15.8E2,14(1X,E20.13E2))")j,time(j),var(j,1),   &
                                   var(j,2),dekdt(j),(var(j,n),n=4,nvar)
    enddo
    close(18)
    print*, ' << ',infile
    !
    stop
    !
  end subroutine tgvstate
  !
  subroutine uprofilesym
    !
    integer :: jm
    !
    parameter(jm=256)
    !
    real(8) :: y(0:jm),u(0:jm)
    !
    real(8) :: uint,udiff,var1,var2
    integer :: j
    !
    open(12,file='uprofile.dat')
    read(12,*)
    do j=0,jm
      read(12,*)y(j),u(j)
    enddo
    close(12)
    print*,' >> uprofile'
    !
    uint=0.d0
    do j=1,jm
      uint=uint+0.5d0*(y(j)-y(j-1))*(u(j)+u(j-1))
    enddo
    print*,'uint=',uint
    !
    udiff=0.d0
    do j=1,jm/2
      var1=abs(u(j-1)-u(jm-j+1))
      var2=abs(u(j)-u(jm-j))
      !
      udiff=udiff+0.5d0*(y(j)-y(j-1))*(var1+var2)
    enddo
    print*,'udiff=',udiff,2.d0*udiff/uint
    !
  end subroutine uprofilesym
  !
  subroutine dis_streamline_wall
    !
    integer :: im,i,j
    real(8),allocatable :: x(:),y(:),ywall(:),dis(:)
    real(8) :: xec,yec,xcc,ycc,var1,var2,a,b,c
    !
    im=16866
    !
    allocate(x(im),y(im),ywall(im),dis(im))
    !
    open(11,file='data/stream_line1.dat')
    do i=1,17
        read(11,*)
    enddo
    do i=1,im
        read(11,*)x(i),y(i)
    enddo
    close(11)
    print*,' >> data/stream_line1.dat'
    !
    xcc=0.d0
    ycc=0.d0
    xec=1.d0
    yec=1.d0
    !
    do i=1,im
        if(x(i)<xcc) then
            ywall(i)=ycc
        elseif(x(i)>xec) then
            ywall(i)=yec
        else
            ywall(i)=x(i)-xcc
        endif
    enddo
    !
    do i=1,im
      !
      dis(i)=1.d10
      do j=1,im
        var1=sqrt((x(i)-x(j))**2+(y(i)-ywall(j))**2)
        dis(i)=min(dis(i),var1)
      enddo
      !
    enddo
    !
    open(16,file='dis_streamline_wall.dat')
    do i=1,im
    write(16,*)x(i),y(i),dis(i)
    enddo
    close(16)
    print*,' << dis_streamline_wall.dat'
    !
  end subroutine dis_streamline_wall
  !
  !+-------------------------------------------------------------------+
  !| This function is to calculate the area of a triangle using        |
  !| Heron's formula.                                                  | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  pure real(8) function areatriangle(a,b,c)
    !
    ! arguments
    real(8),intent(in) :: a,b,c
    !
    ! local data
    real(8) :: var1
    !
    var1=0.5d0*(a+b+c)
    !
    areatriangle=sqrt(var1*(var1-a)*(var1-b)*(var1-c))
    !
    return
    !
  end function areatriangle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine areatriangle.                           |
  !+-------------------------------------------------------------------+
  !
  subroutine expdatint
    !
    use interpolation, only : regularlinearinterp
    !
    integer :: in1,in2,in3,in4,in
    integer :: i
    !
    real(8),allocatable,dimension(:) :: y1,y2,y3,y4,y,ro,uu,vv,ww,u11,u22,u33,roi
    real(8) :: twall,rowall,var1
    !
    in1=17
    in2=38
    in3=21
    in4=26
    !
    allocate(y1(0:in1),ro(0:in1))
    allocate(y2(0:in2),uu(0:in2))
    allocate(y3(0:in3),vv(0:in3))
    allocate(y4(0:in4),ww(0:in4))
    !
    open(13,file='PMdata.dat')
    do i=0,in1
      read(13,*)y1(i),ro(i)
    enddo
    read(13,*)
    do i=0,in2
      read(13,*)y2(i),uu(i)
    enddo
    read(13,*)
    do i=0,in3
      read(13,*)y3(i),vv(i)
    enddo
    read(13,*)
    do i=0,in4
      read(13,*)y4(i),ww(i)
    enddo
    close(13)
    print*,' >> PMdata.dat'
    !
    in=100
    allocate(y(0:in),u11(0:in),u22(0:in),u33(0:in),roi(0:in))
    !
    do i=0,in
      y(i)=1.5d0/in*i
    enddo
    !
    twall=10.254372d0*0.5d0
    rowall=1.d0/twall
    !
    call regularlinearinterp(y1,ro,in1, y,roi,in)
    call regularlinearinterp(y2,uu,in2, y,u11,in)
    call regularlinearinterp(y3,vv,in2, y,u22,in)
    call regularlinearinterp(y4,ww,in2, y,u33,in)
    !
    open(16,file='density_scale_PMdata.dat')
    do i=1,in
      var1=sqrt(roi(i)/rowall)
      ! var1=1.d0
      write(16,'(5(1X,E14.7E2))')y(i),roi(i),u11(i)*var1,u22(i)*var1,u33(i)*var1
    enddo
    close(16)
    print*,' << density_scale_dns.dat'
    !
  end subroutine expdatint
  !
  subroutine mijgen
    !
    integer :: i,j,k,l
    character(len=4) :: name1,name2,name3,name4
    !
    do i=1,3
    do j=1,3
    do k=1,3
        !

      write(name1,'(a,i1,i1)')'du',i,j
      write(name2,'(a,i1,i1)')'du',j,i
      write(name3,'(a,i1,i1)')'du',k,k
      ! write(name4,'(a,i1,i1)')'du',l,l
      write(*,'(A)',advance='no')name1//'*'//name2//'*'//name3//'+'
    enddo
    write(*,*)'&'
    enddo
    enddo
    ! enddo
    !
  end subroutine mijgen
  !
  subroutine momread
    !
    real(8) :: time
    integer :: record,i
    real(8) :: m11m11,m22m22,m33m33,m11m22,m11m33,m22m33, &
               m12m21,m13m31,m23m32,m12m12,m13m13,m23m23, &
               m21m21,m31m31,m32m32,miimjj,mijmij,mijmji
    !
    open(43,file='mom2nd.bin',form='unformatted',access='direct',recl=8*19)
    open(44,file='mom3rd.bin',form='unformatted',access='direct',recl=8*36)
    open(45,file='mom4th.bin',form='unformatted',access='direct',recl=8*9)
    !
    open(16,file='dudx2.dat')
    !
    record=0
    do i=1,10000
      record=record+1
      read(43,rec=record)time,m11m11,m22m22,m33m33,m11m22,m11m33,m22m33,m12m21,m13m31,m23m32,  &
                                 m12m12,m13m13,m23m23,m21m21,m31m31,m32m32,miimjj,mijmij,mijmji
      ! read(44,rec=record)time,m11m11m11,m22m22m22,m33m33m33,m11m11m22,   &
      !                         m11m11m33,m22m22m33,m22m22m11,m33m33m11,   &
      !                         m33m33m22,m11m12m12,m11m13m13,m22m23m23,   &
      !                         m22m21m21,m33m32m32,m33m31m31,m11m12m21,   &
      !                         m11m13m31,m22m23m32,m22m21m12,m33m32m23,   &
      !                         m33m31m13,m11m23m23,m11m32m32,m22m13m13,   &
      !                         m22m31m31,m33m12m12,m33m21m21,m11m23m32,   &
      !                         m22m13m31,m33m12m21
      ! read(45,rec=record)time,miimjjmkkmll,mijmjkmklmli,mijmijmklmkl,    &
      !                         mijmjimklmlk,mijmjimklmkl,mijmjimkkmll,    &
      !                         mijmijmkkmll,mijmjkmkimll
      write(16,*)i,time,m11m11+m22m22+m33m33
    enddo
    close(43)
    close(44)
    close(45)
    !
    close(16)
    print*,' << dudx2.dat'
    !
  end subroutine momread
  !
  subroutine canteratest
    !
    ! use cantera
    ! !
    ! type(phase_t) :: mixture
    ! logical :: lcomb
    ! character(len=255) :: chemfile
    ! character(len=3) :: odetype
    ! character(len=10) :: phase_id
    ! !
    ! phase_id='gas'
    ! chemfile='ES80_H2-8-16.xml'
    ! !
    ! !---CANTERA---
    ! mixture=importPhase(trim(chemfile),trim(phase_id))
    ! print*,mixture,phase_id,trim(chemfile)
    !
  end subroutine canteratest
  !
  subroutine writevelofied(fileinput,fileoutput)
    !
    use commvardefine,only: im,jm,km,u1,u2,u3
    use h5readwrite
    !
    character(len=*),intent(in) :: fileinput,fileoutput
    !
    real(8) :: time
    !
    allocate(u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),                    &
             u3(0:im,0:jm,0:km))
    !
    call H5ReadArray(time,'time',trim(fileinput))
    call H5ReadArray(u1,im,jm,km,'u1',trim(fileinput))
    call H5ReadArray(u2,im,jm,km,'u2',trim(fileinput))
    call H5ReadArray(u3,im,jm,km,'u3',trim(fileinput))
    !
    open(16,file=trim(fileoutput),form='unformatted')
    write(16)im,jm,km
    write(16)time
    write(16)u1(1:im,1:jm,1:km)
    write(16)u2(1:im,1:jm,1:km)
    write(16)u3(1:im,1:jm,1:km)
    close(16)
    print*,' <<',trim(fileoutput)
    !
  end subroutine writevelofied
  !
  subroutine readvelofied(fileinput)
    !
    use commvardefine,only: im,jm,km,u1,u2,u3
    use h5readwrite
    !
    character(len=*),intent(in) :: fileinput
    !
    real(8) :: time
    !
    open(18,file=trim(fileinput),form='unformatted')
    read(18)im,jm,km
    print*,' ** im,jm,km:',im,jm,km
    read(18)time
    print*,' ** time=',time
    allocate(u1(1:im,1:jm,1:km),u2(1:im,1:jm,1:km),                    &
             u3(1:im,1:jm,1:km))
    !
    read(18)u1
    read(18)u2
    read(18)u3
    close(18)
    print*,' >>',trim(fileinput)
    !
  end subroutine readvelofied
  !
  subroutine skewness_cal(fileinput)
    !
    use commvardefine,only: im,jm,km,u1,u2,u3,ro,p,t,mach,x,y,z
    use basicfunction,only: rns3d
    use gradsolver
    use h5readwrite
    use readwrite
    !
    character(len=*),intent(in) :: fileinput
    real(8),allocatable,dimension(:,:,:,:) :: du1,du2,du3
    !
    real(8) :: skew1,skew2,skew3
    !
    logical :: lfex
    integer :: i,j,k
    !
    print*,' ** to view 3d flow field structures.'
    !
    call H5_ReadGrid
    !
    allocate(u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),                    &
             u3(0:im,0:jm,0:km))  
    !
    call H5ReadArray(u1,im,jm,km,'u1',trim(fileinput))
    call H5ReadArray(u2,im,jm,km,'u2',trim(fileinput))
    call H5ReadArray(u3,im,jm,km,'u3',trim(fileinput))
    !
    allocate(du1(1:3,0:im,0:jm,0:km),du2(1:3,0:im,0:jm,0:km),          &
             du3(1:3,0:im,0:jm,0:km)  )
    !
    du1=grad_3d(u1,x,y,z)
    du2=grad_3d(u2)
    du3=grad_3d(u3)
    !
    skew1=0.d0
    skew2=0.d0
    skew3=0.d0
    do k=1,km
    do j=1,jm
    do i=1,im
      skew1=skew1+du1(1,i,j,k)**3
      skew2=skew2+du2(2,i,j,k)**3
      skew3=skew3+du3(3,i,j,k)**3
    enddo
    enddo
    enddo
    skew1=skew1/dble(im*jm*km)
    skew2=skew2/dble(im*jm*km)
    skew3=skew3/dble(im*jm*km)
    !
    print*,' skew1= ',skew1
    print*,' skew2= ',skew2
    print*,' skew3= ',skew3
    !
  end subroutine skewness_cal
  !
  subroutine filter_test
    !
    integer,parameter :: jm=128
    !
    real(8) :: y(0:jm),u(0:jm),v(0:jm),t(0:jm),tke(0:jm),omega(0:jm),miut(0:jm)
    real(8) :: miuf(0:jm)
    integer :: j
    !
    open(12,file='profile_k-omega-6cc.dat')
    read(12,*)
    do j=0,jm
      read(12,*)y(j),u(j),v(j),t(j),tke(j),omega(j),miut(j)
    enddo
    close(12)
    print*,' >> profile_k-omega-6cc.dat'
    !
    miuf=omega
    do j=1,jm-1
      miuf(j)=0.125d0*(omega(j-1)+omega(j+1))+0.75d0*omega(j)
    enddo
    !
    open(16,file='filter_results.dat')
    do j=0,jm
      write(16,*)y(j),omega(j),miuf(j)
    enddo
    close(16)
    print*,' << filter_results'
    !
  end subroutine filter_test
  !
  subroutine tgverror(fileinput)
    !
    use interpolation
    !
    character(len=*),intent(in) :: fileinput
    !
    integer :: nref,nmax
    real(8),allocatable,dimension(:) :: tim_ref,tke_ref,ens_ref,dis_ref
    real(8),allocatable,dimension(:) :: tim_res,tke_res,ens_res,dis_res
    real(8),allocatable,dimension(:) :: tim,tke,ens,dis
    !
    integer :: i,n
    real(8) :: tkemax,ensmax,dismax
    real(8) :: l2error_tke,l2error_ens,l2error_dis
    !
    nref=2000
    allocate(tim_ref(1:nref),tke_ref(1:nref),ens_ref(1:nref),dis_ref(1:nref))
    !
    tkemax=0.d0
    ensmax=0.d0
    dismax=0.d0
    open(12,file='refdat.dat')
    do i=1,nref
      read(12,*)tim_ref(i),tke_ref(i),dis_ref(i),ens_ref(i)
      tkemax=max(tkemax,tke_ref(i))
      ensmax=max(ensmax,ens_ref(i))
      dismax=max(dismax,dis_ref(i))
    enddo
    close(12)
    print*,' >> refdat.dat'
    print*,' ** max kinetic energy= ',tkemax
    print*,' ** max enstrophy     = ',ensmax
    print*,' ** max dissipation   = ',dismax
    !
    nmax=40000
    !
    allocate(tim(0:nmax),tke(0:nmax),ens(0:nmax),dis(0:nmax))
    allocate(tim_res(0:nmax),tke_res(0:nmax),ens_res(0:nmax),dis_res(0:nmax))
    !
    open(12,file=fileinput)
    read(12,*)
    do i=0,nmax
      read(12,*)n,tim(i),tke(i),ens(i),dis(i)
    enddo
    close(12)
    print*,' >> ',fileinput
    print*,' ** max time= ',tim(nmax)
    !
    tim_res=tim
    !
    call regularlinearinterp(tim_ref,tke_ref,ens_ref,dis_ref,nref-1, &
                             tim_res,tke_res,ens_res,dis_res,nmax)
    !
    ! calculate the error
    l2error_tke=0.d0
    l2error_ens=0.d0
    l2error_dis=0.d0
    do i=1,nmax
      l2error_tke=l2error_tke+(tke_res(i)-tke(i))**2
      l2error_ens=l2error_ens+(ens_res(i)-ens(i))**2
      l2error_dis=l2error_dis+(dis_res(i)-dis(i))**2
    enddo
    l2error_tke=sqrt(l2error_tke/dble(nmax))/tkemax
    l2error_ens=sqrt(l2error_ens/dble(nmax))/ensmax
    l2error_dis=sqrt(l2error_dis/dble(nmax))/dismax
    print*,' ** L2 error of kinetic energy= ',l2error_tke
    print*,' ** L2 error of enstrophy     = ',l2error_ens
    print*,' ** L2 error of dissipation   = ',l2error_dis
    !
  end subroutine tgverror
  !
  subroutine ijk2intone
    !
    integer(8) :: ijk
    integer :: i,j,k
    integer, parameter :: i16 = selected_int_kind(16)
    !
    i=9996
    j=1024
    k=5534
    !
    print*,i,j,k
    !
    ijk=i+j*1000000+k*1000000000000_i16
    !
    print*,ijk
    !
    k=ijk/1000000000000_i16
    j=ijk/1000000-k*1000000
    i=ijk-k*1000000000000_i16-j*1000000
    !
    print*,i,j,k
    !
  end subroutine ijk2intone
  !
  subroutine particle_gen
    !
    use commvardefine, only: pi
    use h5readwrite
    !
    integer :: numpart,im,jm,km,i,j,k,n,itime
    real(8) :: xmin,xmax,ymin,ymax,zmin,zmax,time
    real(8),allocatable,dimension(:) :: x,y,z
    !
    xmin=0.01d0
    xmax=4.d0*pi-0.01d0
    ymin=0.01d0
    ymax=2.d0-0.01d0
    zmin=0.01d0
    zmax=2.d0*pi-0.01d0
    !
    im=100
    jm=100
    km=100
    !
    numpart=im*jm*km
    !
    allocate(x(numpart),y(numpart),z(numpart))
    !
    n=0
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      n=n+1
      !
      x(n)=(xmax-xmin)/dble(im-1)*dble(i-1)+xmin
      ! x(n)=xmin
      y(n)=(ymax-ymin)/dble(jm-1)*dble(j-1)+ymin
      z(n)=(zmax-zmin)/dble(km-1)*dble(k-1)+zmin
      !
    enddo
    enddo
    enddo
    !
    itime=0
    time=0.d0
    !
    call H5WriteArray(varin=itime,vname='itime',fname='particle00000.h5')
    call H5WriteArray(varin=time, vname='time', fname='particle00000.h5')
    call H5WriteArray(varin=x,    vname='x',    fname='particle00000.h5', dim1=numpart-1)
    call H5WriteArray(varin=y,    vname='y',    fname='particle00000.h5', dim1=numpart-1)
    call H5WriteArray(varin=z,    vname='z',    fname='particle00000.h5', dim1=numpart-1)
    !
  end subroutine particle_gen
  !
  subroutine particle_stat(nsta,nend)
    !
    use h5readwrite
    !
    integer,intent(in) :: nsta,nend
    !
    integer :: numpart,im,jm,km,i,j,k,n,itime,jp
    real(8) :: time,var1
    character(len=5) :: fnum
    character(len=32) :: fname
    real(8),allocatable :: y(:),yc(:),uc(:),vc(:),wc(:)
    integer,allocatable :: couter(:)
    real(8),allocatable,dimension(:) :: xp,yp,zp,up,vp,wp,pdf
    !
    jm=200
    !
    allocate(y(0:jm),yc(1:jm),couter(1:jm))
    allocate(uc(1:jm),vc(1:jm),wc(1:jm))
    allocate(pdf(1:jm))
    !
    do j=0,jm
      y(j)=2.d0/dble(jm)*dble(j)
    enddo
    do j=1,jm
      yc(j)=0.5d0*(y(j)+y(j-1))
    enddo
    !
    couter=0
    uc=0.d0
    vc=0.d0
    wc=0.d0
    !
    do n=nsta,nend
      !
      write(fnum,'(I5.5)')n
      fname='./data/particle'//fnum//'.h5'
      !
      numpart=h5getdim1d(varname='x',filenma=trim(fname))
      !
      allocate(xp(numpart),yp(numpart),zp(numpart),up(numpart),vp(numpart),wp(numpart))
      !
      call H5ReadArray(varin=itime,vname='itime',fname=trim(fname))
      call H5ReadArray(varin=time, vname='time', fname=trim(fname))
      call H5ReadArray(varin=xp,   vname='x',    fname=trim(fname), dim1=numpart-1)
      call H5ReadArray(varin=yp,   vname='y',    fname=trim(fname), dim1=numpart-1)
      call H5ReadArray(varin=zp,   vname='z',    fname=trim(fname), dim1=numpart-1)
      call H5ReadArray(varin=up,   vname='u',    fname=trim(fname), dim1=numpart-1)
      call H5ReadArray(varin=vp,   vname='v',    fname=trim(fname), dim1=numpart-1)
      call H5ReadArray(varin=wp,   vname='w',    fname=trim(fname), dim1=numpart-1)
      !
      do jp=1,numpart
        !
        do j=1,jm
          !
          if(yp(jp)>=y(j-1) .and. yp(jp)<=y(j)) then
            !
            couter(j)=couter(j)+1
            !
            uc(j)=uc(j)+up(jp)
            vc(j)=vc(j)+vp(jp)
            wc(j)=wc(j)+wp(jp)
            !
            exit
            !
          endif
          !
        enddo
        !
      enddo
      !
      deallocate(xp,yp,zp,up,vp,wp)
      !
    enddo
    !
    do j=1,jm
      uc(j)=uc(j)/dble(couter(j))
      vc(j)=vc(j)/dble(couter(j))
      wc(j)=wc(j)/dble(couter(j))
    enddo
    !
    var1=0.d0
    do j=1,jm
      var1=var1+dble(couter(j)) !/(y(j)-y(j-1))
    enddo
    do j=1,jm
      pdf(j)=dble(couter(j))/var1
    enddo
    !
    open(18,file='profile_particles.dat')
    write(18,'(5(1X,A14))')'y','u','v','w','n'
    do j=1,jm
      write(18,'(5(1X,E14.7E2))')yc(j),uc(j),vc(j),wc(j),pdf(j)
    enddo
    close(18)
    print*,' << profile.dat'
    !
  end subroutine particle_stat
  !
  function int_to_str(i)
    integer, intent(in) :: i
    character(len=(1 + int(log10(real(i))))) :: int_to_str

    ! print*,1 + int(log10(real(i)))

    write(int_to_str, "(I0)") i
  end function int_to_str
  !
  subroutine incompact3d_visu(nsta,nend,input,path)
    !
    use commvardefine,only: im,jm,km
    use h5readwrite
    use WriteTec
    use basicfunction
    use LinearAlgegra
    use gradsolver
    !
    integer,intent(in) :: nsta,nend
    character(len=*),intent(in) :: input,path
    !
    integer :: numpart,nx,ny,nz,i,j,k,n,itime,jp,rc,fu,id,id2,nsamples
    real(8) :: time,var1,var2,var3,re,re_tau,Ha,Stuart,Rem,xlx,yly,zlz
    character(len=1) :: cc
    character(len=4) :: fnum
    character(len=7) :: stanum
    character(len=64) :: fname
    character(len=64) :: str

    character(len=16) :: flowtype

    logical :: lfilext
    !
    real(8),allocatable,dimension(:,:) :: je
    real(8),allocatable,dimension(:,:,:) :: x,y,z,u,v,w,q,bx,by,bz,vort,db1,db2,db3,je3d
    real(8),allocatable,dimension(:,:,:,:) :: bm,dbx,dby,dbz
    !
    namelist /BasicParam/ nx,ny,nz
    !
    ! Open and read Namelist file.
    open(action='read', file=input, iostat=rc, newunit=fu)
    !
    do while(rc==0)
      read(fu,'(a)',iostat=rc)str
      !
      id=index(str,'nx',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+5),*)nx
        endif
      endif
      !
      id=index(str,'ny',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+5),*)ny
        endif
      endif
      !
      id=index(str,'nz',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+5),*)nz
        endif
      endif
      !
      id=index(str,'xlx',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+5),*)xlx
        endif
      endif
      !
      id=index(str,'yly',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+5),*)yly
        endif
      endif
      !
      id=index(str,'zlz',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+5),*)zlz
        endif
      endif
      !
      id=index(str,'re',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+17),*)re
        endif
      endif
      !
      id=index(str,'hartmann',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+17),*)Ha
        endif
      endif
      !
      id=index(str,'Rem',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+17),*)Rem
        endif
      endif
      !
      id=index(str,'Stuart',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+17),*)Stuart
        endif
      endif
      !
    enddo
    !
    write(*,*)'    re: ',re
    write(*,*)'    Ha: ',Ha
    write(*,*)'Stuart: ',Stuart
    write(*,*)'   Rem: ',Rem
    !
    flowtype='tgv'
    !
    if(trim(flowtype)=='channel') then
      !
      ny=ny-1
      !
      write(*,'(3(A4))')'nx','ny','nz'
      write(*,'(3(I4))')nx,ny,nz
      !
      im=nx
      jm=ny
      km=nz
      !
      ! mesh generation
      allocate(x(0:nx,0:ny,0:nz),y(0:nx,0:ny,0:nz),z(0:nx,0:ny,0:nz))
      do k=0,nz
      do j=0,ny
      do i=0,nx
        x(i,j,k)=8.d0/dble(nx)*dble(i)
        z(i,j,k)=4.d0/dble(nz)*dble(k)
      enddo
      enddo
      enddo
      !
      open(18,file='yp.dat')
      do j=0,ny
        read(18,*)y(1,j,1)
      enddo
      close(18)
      print*,' >> yp.dat'
      !
      do k=0,nz
      do i=0,nx
      do j=0,ny
        y(i,j,k)=y(1,j,1)
      enddo
      enddo
      enddo
      !
      allocate(u(0:nx,0:ny,0:nz),q(0:nx,0:ny,0:nz))
      do n=nsta,nend
        !
        write(fnum,'(I4.4)')n
        !
        fname=path//'/ux-'//int_to_str(n)//'.bin'
        open(16,file=trim(fname),access="stream")
        read(16)u(1:nx,0:ny,1:nz)
        close(16)
        print*,' << ',trim(fname)
        !
        call bcchan(u)
        !
        fname=path//'/critq-'//int_to_str(n)//'.bin'
        open(16,file=trim(fname),access="stream")
        read(16)q(1:nx,0:ny,1:nz)
        close(16)
        print*,' << ',trim(fname)
        !
        call bcchan(q)
        !
        call writetecbin('teccs'//fnum//'.plt',x,'x',y,'y',z,'z',u,'u',q,'Q',nx,ny,nz)
   
      enddo

    elseif(trim(flowtype)=='tgv') then
      !
      write(*,'(3(A4))')'nx','ny','nz'
      write(*,'(3(I4))')nx,ny,nz
      !
      im=nx
      jm=ny
      km=nz
      !
      print*,xlx,yly,zlz
      !
      allocate(x(0:nx,0:ny,0:nz),y(0:nx,0:ny,0:nz),z(0:nx,0:ny,0:nz))
      do k=0,nz
      do j=0,ny
      do i=0,nx
        x(i,j,k)=xlx/dble(nx)*dble(i)
        y(i,j,k)=yly/dble(ny)*dble(j)
        z(i,j,k)=zlz/dble(nz)*dble(k)
      enddo
      enddo
      enddo
      !
      allocate(u(0:nx,0:ny,0:nz),v(0:nx,0:ny,0:nz),w(0:nx,0:ny,0:nz),vort(0:nx,0:ny,0:nz),q(0:nx,0:ny,0:nz))
      allocate(bm(0:nx,0:ny,0:nz,1:3),db1(1:2,0:nx,0:ny),db2(1:2,0:nx,0:ny),db3(1:2,0:nx,0:ny),je(0:nx,0:ny))
      allocate(je3d(0:nx,0:ny,0:nz),dbx(1:3,0:nx,0:ny,0:nz),dby(1:3,0:nx,0:ny,0:nz),dbz(1:3,0:nx,0:ny,0:nz))
      do n=nsta,nend
        !
        write(fnum,'(I4.4)')n
        !
        ! fname='data/ux-'//int_to_str(n)//'.bin'
        ! open(16,file=trim(fname),access="stream")
        ! read(16)u(1:nx,1:ny,1:nz)
        ! close(16)
        ! print*,' << ',trim(fname)
        ! !
        ! call bchomo(u)
        ! !
        ! fname='data/uy-'//int_to_str(n)//'.bin'
        ! open(16,file=trim(fname),access="stream")
        ! read(16)v(1:nx,1:ny,1:nz)
        ! close(16)
        ! print*,' << ',trim(fname)
        ! !
        ! call bchomo(v)
        ! !
        ! fname='data/uz-'//int_to_str(n)//'.bin'
        ! open(16,file=trim(fname),access="stream")
        ! read(16)w(1:nx,1:ny,1:nz)
        ! close(16)
        ! print*,' << ',trim(fname)
        ! !
        ! call bchomo(w)
        !
        fname='data/critq-'//int_to_str(n)//'.bin'
        open(16,file=trim(fname),access="stream")
        read(16)q(1:nx,1:ny,1:nz)
        close(16)
        call bchomo(q)
        print*,' << ',trim(fname)
        !
        fname='data/vort-'//int_to_str(n)//'.bin'
        open(16,file=trim(fname),access="stream")
        read(16)vort(1:nx,1:ny,1:nz)
        close(16)
        call bchomo(vort)
        print*,' << ',trim(fname)
        !
        call bchomo(w)
        !
        fname='data/B_x-'//int_to_str(n)//'.bin'
        !
        inquire(file=trim(fname),exist=lfilext)
        !
        if(lfilext) then
          !
          open(16,file=trim(fname),access="stream")
          read(16)bm(1:nx,1:ny,1:nz,1)
          close(16)
          print*,' << ',trim(fname)
          !
          call bchomo(bm(:,:,:,1))
          !
          fname='data/B_y-'//int_to_str(n)//'.bin'
          open(16,file=trim(fname),access="stream")
          read(16)bm(1:nx,1:ny,1:nz,2)
          close(16)
          print*,' << ',trim(fname)
          !
          call bchomo(bm(:,:,:,2))
          !
          fname='data/B_z-'//int_to_str(n)//'.bin'
          open(16,file=trim(fname),access="stream")
          read(16)bm(1:nx,1:ny,1:nz,3)
          close(16)
          print*,' << ',trim(fname)
          !
          call bchomo(bm(:,:,:,3))
          !
          db1=grad_xy(bm(:,:,0,1),x(:,:,0),y(:,:,0))
          db2=grad_xy(bm(:,:,0,2))
          db3=grad_xy(bm(:,:,0,3))
          !
          je=db2(1,:,:)-db1(2,:,:)
          !
          ! dbx=grad_3d(bm(:,:,:,1),x,y,z)
          ! dby=grad_3d(bm(:,:,:,2))
          ! dbz=grad_3d(bm(:,:,:,3))
          ! !
          ! je3d=omega(dbx,dby,dbz,'m')
          !
        endif
        !
        call writetecbin('tec2d'//fnum//'.plt',x(:,:,0),'x', &
                                               y(:,:,0),'y', &
                                               u(:,:,0),'u', &
                                               v(:,:,0),'v', &
                                           bm(:,:,0,1),'bx', &
                                           bm(:,:,0,2),'by', &
                                            vort(:,:,0),'vort', &
                                            je,'j',nx,ny)
        ! call writetecbin('tec3d'//fnum//'.plt',x,'x',y,'y',z,'z', &
        !                                    vort,'vort',q,'q',je3d,'j',nx,ny,nz)
        !
      enddo
      !
      stop
      !
    else
      print*,trim(flowtype)
      stop ' flowtype not defined '
    endif
    !
  end subroutine incompact3d_visu
  !
  subroutine incompact3d(nsta,nend,input,mode)
    !
    use commvardefine,only: im,jm,km
    use h5readwrite
    use WriteTec
    use basicfunction
    use LinearAlgegra
    use gradsolver
    !
    integer,intent(in) :: nsta,nend
    character(len=*),intent(in) :: input,mode
    !
    integer :: numpart,nx,ny,nz,i,j,k,n,itime,jp,rc,fu,id,id2,nsamples,istret
    real(8) :: time,var1,re,re_tau,Ha,Stuart,Rem,xlx,yly,zlz
    character(len=1) :: cc
    character(len=4) :: fnum
    character(len=7) :: stanum
    character(len=32) :: fname
    character(len=64) :: str
    real(8),allocatable,dimension(:,:,:) :: x,y,z,p,u,v,w,u11,u12,u13,u22,u23,u33
    real(8),allocatable,dimension(:,:) :: uyz,vyz,wyz
    real(8),allocatable,dimension(:) :: pm,um,vm,wm,uplus,yplus,uha,jxm,jym,jzm
    real(8),allocatable,dimension(:) :: uu,vv,ww,uv,pp
    real(8),allocatable,dimension(:,:) :: force,elec
    real(8),allocatable :: db2d(:,:,:),f2d(:,:,:),du(:,:,:),dv(:,:,:),omegaz(:,:),babs(:,:)
    real(8),allocatable,dimension(:,:,:,:) :: je,bm 
    real(8) :: magen(3),vec1(3),vec2(3)
    real(8) :: utaw,Re_bulk,u_bulk,ez_avg
    logical :: lfilext,lfilextj
    !
    ! Open and read Namelist file.
    open(action='read', file=input, iostat=rc, newunit=fu)
    !
    do while(rc==0)
      read(fu,'(a)',iostat=rc)str
      !
      id=index(str,'nx',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+5),*)nx
        endif
      endif
      !
      id=index(str,'ny',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+5),*)ny
        endif
      endif
      !
      id=index(str,'nz',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+5),*)nz
        endif
      endif
      !
      id=index(str,'xlx',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+5),*)xlx
        endif
      endif
      !
      id=index(str,'yly',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+5),*)yly
        endif
      endif
      !
      id=index(str,'zlz',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+5),*)zlz
        endif
      endif
      !
      id=index(str,'istret',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+5),*)istret
        endif
      endif
      !
      id=index(str,'re',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+17),*)re
        endif
      endif
      !
      id=index(str,'hartmann',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+17),*)Ha
        endif
      endif
      !
      id=index(str,'Rem',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+17),*)Rem
        endif
      endif
      !
      id=index(str,'Stuart',back=.false.)
      if(id==1) then
        id2=index(str,'=',back=.false.)
        if(id2>id) then
          read(str(id2+1:id2+17),*)Stuart
        endif
      endif
      !
    enddo
    !
    ny=ny-1
    !
    write(*,'(3(A4))')'nx','ny','nz'
    write(*,'(3(I4))')nx,ny,nz
    !
    im=nx
    jm=ny
    km=nz
    !
    write(*,*)'    re: ',re
    write(*,*)'    Ha: ',Ha
    write(*,*)'Stuart: ',Stuart
    write(*,*)'   Rem: ',Rem
    ! nx=256
    ! ny=128
    ! nz=256
    !
    ! mesh generation
    allocate(x(0:nx,0:ny,0:nz),y(0:nx,0:ny,0:nz),z(0:nx,0:ny,0:nz))
    do k=0,nz
    do j=0,ny
    do i=0,nx
      x(i,j,k)=xlx/dble(nx)*dble(i)
      z(i,j,k)=zlz/dble(nz)*dble(k)
    enddo
    enddo
    enddo
    !
    if(istret > 0) then
      open(18,file='yp.dat')
      do j=0,ny
        read(18,*)y(1,j,1)
      enddo
      close(18)
      print*,' >> yp.dat'
      !
      do k=0,nz
      do i=0,nx
      do j=0,ny
        y(i,j,k)=y(1,j,1)
      enddo
      enddo
      enddo
    else
      do k=0,nz
      do i=0,nx
      do j=0,ny
        y(i,j,k)=yly/dble(ny)*dble(j)
      enddo
      enddo
      enddo
    endif
    ! do k=0,nz
    ! do j=0,ny
    ! do i=0,nx
    !   x(i,j,k)=xlx/dble(nx)*dble(i)
    !   y(i,j,k)=yly/dble(ny)*dble(j)
    !   z(i,j,k)=zlz/dble(nz)*dble(k)
    ! enddo
    ! enddo
    ! enddo
    !
    allocate(p(0:nx,0:ny,0:nz),u(0:nx,0:ny,0:nz),v(0:nx,0:ny,0:nz),w(0:nx,0:ny,0:nz))
    allocate(u11(0:nx,0:ny,0:nz),u12(0:nx,0:ny,0:nz),u13(0:nx,0:ny,0:nz),u22(0:nx,0:ny,0:nz), &
             u23(0:nx,0:ny,0:nz),u33(0:nx,0:ny,0:nz))
    allocate(pm(0:ny),um(0:ny),vm(0:ny),wm(0:ny))
    allocate(uu(0:ny),vv(0:ny),ww(0:ny),uv(0:ny),pp(0:ny))
    !
    allocate(je(0:nx,0:ny,0:nz,1:3))
    allocate(bm(0:nx,0:ny,0:nz,1:3))
    allocate(jxm(0:ny),jym(0:ny),jzm(0:ny))
    !
    allocate(db2d(1:2,0:nx,0:ny),f2d(0:nx,0:ny,1:3),du(1:2,0:nx,0:ny),dv(1:2,0:nx,0:ny),omegaz(0:nx,0:ny),babs(0:nx,0:ny))
    !
    allocate(uyz(0:ny,0:nz),vyz(0:ny,0:nz),wyz(0:ny,0:nz))
    !
    pm=0.d0
    !
    um=0.d0
    vm=0.d0
    wm=0.d0
    !
    pp=0.d0
    uu=0.d0
    vv=0.d0
    ww=0.d0
    uv=0.d0
    !
    jxm=0.d0
    jym=0.d0
    jzm=0.d0
    !
    uyz=0.d0
    vyz=0.d0
    wyz=0.d0
    !
    if(trim(mode)=='readinst') then
      !
      do n=nsta,nend
        !
        write(fnum,'(I4.4)')n
        !
        fname='data/ux-'//int_to_str(n)//'.bin'
        open(16,file=trim(fname),access="stream")
        read(16)u(1:nx,:,1:nz)
        close(16)
        print*,' << ',trim(fname)
        !
        call bcchan(u)
        !
        fname='data/uy-'//int_to_str(n)//'.bin'
        open(16,file=trim(fname),access="stream")
        read(16)v(1:nx,:,1:nz)
        close(16)
        print*,' << ',trim(fname)
        !
        call bcchan(v)
        !
        fname='data/uz-'//int_to_str(n)//'.bin'
        open(16,file=trim(fname),access="stream")
        read(16)w(1:nx,:,1:nz)
        close(16)
        print*,' << ',trim(fname)
        !
        call bcchan(w)
        !
        fname='data/pp-'//int_to_str(n)//'.bin'
        open(16,file=trim(fname),access="stream")
        read(16)p(1:nx,:,1:nz)
        close(16)
        print*,' << ',trim(fname)
        !
        call bcchan(p)
        !
        fname='data/J_x-'//int_to_str(n)//'.bin'
        !
        inquire(file=trim(fname),exist=lfilextj)
        !
        if(lfilextj) then
          !
          open(16,file=trim(fname),access="stream")
          read(16)je(1:nx,:,1:nz,1)
          close(16)
          print*,' << ',trim(fname)
          !
          call bcchan(je(:,:,:,1))
          !
          fname='data/J_y-'//int_to_str(n)//'.bin'
          open(16,file=trim(fname),access="stream")
          read(16)je(1:nx,:,1:nz,2)
          close(16)
          print*,' << ',trim(fname)
          !
          call bcchan(je(:,:,:,2))
          !
          fname='data/J_z-'//int_to_str(n)//'.bin'
          open(16,file=trim(fname),access="stream")
          read(16)je(1:nx,:,1:nz,3)
          close(16)
          print*,' << ',trim(fname)
          !
          call bcchan(je(:,:,:,3))
          !
        endif
        !
        fname='data/B_x-'//int_to_str(n)//'.bin'
        !
        inquire(file=trim(fname),exist=lfilext)
        !
        if(lfilext) then
          !
          open(16,file=trim(fname),access="stream")
          read(16)bm(1:nx,:,1:nz,1)
          close(16)
          print*,' << ',trim(fname)
          !
          call bcchan(bm(:,:,:,1))
          !
          fname='data/B_y-'//int_to_str(n)//'.bin'
          open(16,file=trim(fname),access="stream")
          read(16)bm(1:nx,:,1:nz,2)
          close(16)
          print*,' << ',trim(fname)
          !
          call bcchan(bm(:,:,:,2))
          !
          fname='data/B_z-'//int_to_str(n)//'.bin'
          open(16,file=trim(fname),access="stream")
          read(16)bm(1:nx,:,1:nz,3)
          close(16)
          print*,' << ',trim(fname)
          !
          call bcchan(bm(:,:,:,3))
          !
          db2d=grad_xy(bm(:,:,0,1),x(:,:,0),y(:,:,0))
          !
          do j=0,ny
          do i=0,nx
            je(i,j,0,2)= 0.d0
            je(i,j,0,3)= - db2d(2,i,j)
          enddo
          enddo
          !
          db2d=grad_xy(bm(:,:,0,2))
          !
          do j=0,ny
          do i=0,nx
            je(i,j,0,1)= 0.d0
            je(i,j,0,3)= je(i,j,0,3) + db2d(1,i,j)
          enddo
          enddo
          !
          do j=0,ny
          do i=0,nx
            f2d(i,j,:)=cross_product(je(i,j,0,:),bm(i,j,0,:))*Stuart
          enddo
          enddo
          !
          du=grad_xy(u(:,:,0))
          dv=grad_xy(v(:,:,0))
          omegaz=dv(1,:,:)-du(2,:,:)
          !
          babs=sqrt(bm(:,:,0,1)**2+bm(:,:,0,2)**2)
          !
          call writetecbin('tecxy'//fnum//'.plt',x(:,:,0),'x',y(:,:,0),'y',u(:,:,0),'u',  &
                                                                           v(:,:,0),'v',  &
                                                                           p(:,:,0),'p',  &
                                                                           omegaz,'omegaz', &
                                                                           je(:,:,0,3),'Jz',&
                                                                           babs,'|B|',nx,ny)
          !
        endif
        !
        do k=1,nz
        do i=1,nx
        do j=0,ny
          !
          um(j)=um(j)+u(i,j,k)
          vm(j)=vm(j)+v(i,j,k)
          wm(j)=wm(j)+w(i,j,k)
          !
          pm(j)=pm(j)+p(i,j,k)
          !
          uu(j)=uu(j)+u(i,j,k)*u(i,j,k)
          vv(j)=vv(j)+v(i,j,k)*v(i,j,k)
          ww(j)=ww(j)+w(i,j,k)*w(i,j,k)
          uv(j)=uv(j)+u(i,j,k)*v(i,j,k)
          !
          pp(j)=pp(j)+p(i,j,k)*p(i,j,k)
        enddo
        enddo
        enddo
        !
        print*,' ** ww(0)',ww(0), w(nx/2,0,nz/2)
        !
        if(lfilextj) then
          !
          do k=1,nz
          do i=1,nx
          do j=0,ny
            !
            jxm(j)=jxm(j)+je(i,j,k,1)
            jym(j)=jym(j)+je(i,j,k,2)
            jzm(j)=jzm(j)+je(i,j,k,3)
            !
          enddo
          enddo
          enddo
          !
        endif
        !
      enddo
      !
      um=um/dble(nx*nz*(nend-nsta+1))
      vm=vm/dble(nx*nz*(nend-nsta+1))
      wm=wm/dble(nx*nz*(nend-nsta+1))
      pm=pm/dble(nx*nz*(nend-nsta+1))
      !
      uu=uu/dble(nx*nz*(nend-nsta+1))-um*um
      vv=vv/dble(nx*nz*(nend-nsta+1))-vm*vm
      ww=ww/dble(nx*nz*(nend-nsta+1))-wm*wm
      uv=uv/dble(nx*nz*(nend-nsta+1))-um*vm
      !
      pp=pp/dble(nx*nz*(nend-nsta+1))-pm*pm
      !
      jxm=jxm/dble(nx*nz*(nend-nsta+1))
      jym=jym/dble(nx*nz*(nend-nsta+1))
      jzm=jzm/dble(nx*nz*(nend-nsta+1))
      !
      ! call writetecbin('tecwall'//fnum//'.plt',x(:,0,:),'x',y(:,0,:),'y',z(:,0,:),'z',u(:,0,:),'u',v(:,0,:),'v',w(:,0,:),'w',nx,nz)
      !
    elseif(trim(mode)=='readmean') then
      !
      do n=nsta,nend,10000
        !
      write(stanum,'(I7.7)')n
      !
      fname='statistics/pmean.dat'//int_to_str(n)
      open(16,file=trim(fname),access="stream")
      read(16)p(1:nx,:,1:nz)
      close(16)
      print*,' << ',trim(fname)
      call bcchan(p)
      !
      fname='statistics/umean.dat'//int_to_str(n)
      open(16,file=trim(fname),access="stream")
      read(16)u(1:nx,:,1:nz)
      close(16)
      print*,' << ',trim(fname)
      call bcchan(u)
      !
      fname='statistics/vmean.dat'//int_to_str(n)
      open(16,file=trim(fname),access="stream")
      read(16)v(1:nx,:,1:nz)
      close(16)
      print*,' << ',trim(fname)
      call bcchan(v)
      !
      fname='statistics/wmean.dat'//int_to_str(n)
      open(16,file=trim(fname),access="stream")
      read(16)w(1:nx,:,1:nz)
      close(16)
      print*,' << ',trim(fname)
      call bcchan(w)
      !
      fname='statistics/uumean.dat'//int_to_str(n)
      open(16,file=trim(fname),access="stream")
      read(16)u11(1:nx,:,1:nz)
      close(16)
      print*,' << ',trim(fname)
      call bcchan(u11)
      !
      fname='statistics/uvmean.dat'//int_to_str(n)
      open(16,file=trim(fname),access="stream")
      read(16)u12(1:nx,:,1:nz)
      close(16)
      print*,' << ',trim(fname)
      call bcchan(u12)
      !
      fname='statistics/uwmean.dat'//int_to_str(n)
      open(16,file=trim(fname),access="stream")
      read(16)u13(1:nx,:,1:nz)
      close(16)
      print*,' << ',trim(fname)
      call bcchan(u13)
      !
      fname='statistics/vvmean.dat'//int_to_str(n)
      open(16,file=trim(fname),access="stream")
      read(16)u22(1:nx,:,1:nz)
      close(16)
      print*,' << ',trim(fname)
      call bcchan(u22)
      !
      fname='statistics/vwmean.dat'//int_to_str(n)
      open(16,file=trim(fname),access="stream")
      read(16)u23(1:nx,:,1:nz)
      close(16)
      print*,' << ',trim(fname)
      call bcchan(u23)
      !
      fname='statistics/wwmean.dat'//int_to_str(n)
      open(16,file=trim(fname),access="stream")
      read(16)u33(1:nx,:,1:nz)
      close(16)
      print*,' << ',trim(fname)
      call bcchan(u33)
      !
      nsamples=1
      ! open(12,file='nsamples.txt')
      ! read(12,*)nsamples
      ! close(12)
      ! print*,' ** nsamples=',nsamples
      !
      do k=1,nz
      do i=1,nx
      do j=0,ny
        !
        pm(j)=pm(j)+p(i,j,k)
        !
        um(j)=um(j)+u(i,j,k)
        vm(j)=vm(j)+v(i,j,k)
        wm(j)=wm(j)+w(i,j,k)
        !
        uu(j)=uu(j)+u11(i,j,k)
        vv(j)=vv(j)+u22(i,j,k)
        ww(j)=ww(j)+u33(i,j,k)
        uv(j)=uv(j)+u12(i,j,k)
        !
      enddo
      enddo
      enddo
      !
      pm=pm/dble(nsamples*nx*nz)
      !
      um=um/dble(nsamples*nx*nz)
      vm=vm/dble(nsamples*nx*nz)
      wm=wm/dble(nsamples*nx*nz)
      !
      uu=uu/dble(nsamples*nx*nz)-um*um
      vv=vv/dble(nsamples*nx*nz)-vm*vm
      ww=ww/dble(nsamples*nx*nz)-wm*wm
      uv=uv/dble(nsamples*nx*nz)-um*vm
      !
      do k=0,nz
      do j=0,ny
        do i=1,nx
          uyz(j,k)=uyz(j,k)+u(i,j,k)
          vyz(j,k)=vyz(j,k)+v(i,j,k)
          wyz(j,k)=wyz(j,k)+w(i,j,k)
        enddo
      enddo
      enddo
      !
      uyz=uyz/dble(nsamples*nx)
      vyz=vyz/dble(nsamples*nx)
      wyz=wyz/dble(nsamples*nx)
      !
      call writetecbin('tecyz'//stanum//'.plt',y(0,:,:),'y',z(0,:,:),'z',uyz,'u',vyz,'v',wyz,'w',ny,nz)
      !
      enddo
      !
    endif
    !
    open(18,file='uprofile.dat')
    write(18,"(5(1X,A15))")'y','u','v','w','p'
    do j=0,ny
      write(18,"(5(1X,E15.7E3))")y(0,j,0),um(j),vm(j),wm(j),pm(j)
    end do
    close(18)
    print*,' << uprofile.dat'
    !
    u_bulk=0.d0
    do j=1,ny
      u_bulk=u_bulk+0.5d0*(um(j)+um(j-1))*(y(1,j,1)-y(1,j-1,1))
    enddo
    u_bulk=u_bulk/2.d0
    print*,' ** u_bulk=',u_bulk
    !
    allocate(uplus(0:ny/2),yplus(0:ny/2))
    !
    call upluscal_imc(uplus,yplus,um(0:ny/2),y(0,0:ny/2,0),utaw=utaw,Re=re)
    !
    re_tau=utaw*re
    Re_bulk=u_bulk*re
    !
    print*,' ** Re      =',Re
    print*,' ** Re_bulk =',Re_bulk
    print*,' ** Re_tau =',re_tau
    !
    open(18,file='uplus.dat')
    write(18,"(2(1X,A15))")'yplus','uplus'
    write(18,"(2(1X,E15.7E3))")(yplus(j),uplus(j),j=0,ny/2)
    close(18)
    print*,' << uplus.dat ... done. '
    !
    open(18,file='resprofile.dat')
    write(18,"(6(1X,A15))")'y','uu','vv','ww','uv','pp'
    do j=0,ny/2
      write(18,"(6(1X,E15.7E3))")yplus(j),uu(j)/utaw/utaw,vv(j)/utaw/utaw,ww(j)/utaw/utaw,uv(j)/utaw/utaw,pp(j)
    end do
    close(18)
    print*,' << resprofile.dat'
    !
    ! stop
    !
    allocate(force(0:ny,1:3),elec(0:ny,1:3))
    !
    magen(1)=0.d0
    magen(2)=1.d0
    magen(3)=0.d0
    !
    do j=0,ny
      vec1(1)=jxm(j) !+ 0.d0
      vec1(2)=jym(j) !+ 0.d0
      vec1(3)=jzm(j) !- u_bulk*magen(2)
      !
      force(j,:)=cross_product(vec1,magen)*Ha**2/re
    enddo
    !
    do j=0,ny
      vec1(1)=jxm(j)
      vec1(2)=jym(j)
      vec1(3)=jzm(j)
      !
      vec2(1)=um(j)
      vec2(2)=vm(j)
      vec2(3)=wm(j)
      !
      elec(j,:)=vec1-cross_product(vec2,magen)
    enddo
    !
    open(18,file='eprofile.dat')
    write(18,"(4(1X,A15))")'y','Ex','Ey','Ez'
    do j=0,ny
      write(18,"(4(1X,E15.7E3))")y(0,j,0),elec(j,1),elec(j,2),elec(j,3)
    end do
    close(18)
    print*,' << eprofile.dat'
    !
    ez_avg=0.d0
    do j=1,ny
      ez_avg=ez_avg+0.5d0*(elec(j,3)+elec(j-1,3))*(y(1,j,1)-y(1,j-1,1))
    enddo
    print*,' ** Ez_avg =',ez_avg
    !
    open(18,file='jprofile.dat')
    write(18,"(7(1X,A15))")'y','jx','jy','jz','fx','fy','fz'
    do j=0,ny
      write(18,"(7(1X,E15.7E3))")y(0,j,0),jxm(j),jym(j),jzm(j),force(j,1),force(j,2),force(j,3)
    end do
    close(18)
    print*,' << jprofile.dat'
    !
    allocate(uha(0:ny))
    !
    print*,'                          Hatmann number:',Ha
    print*,' Coefficient for the analytical solution:',utaw*re_tau
    do j=0,ny
      uha(j)=utaw*re_tau/(Ha*tanh(Ha))*(1.d0-cosh(Ha*(y(0,j,0)-1.d0))/cosh(Ha))
    enddo
    !
    open(18,file='profile_Ha.dat')
    write(18,"(3(1X,A15),F12.6)")'y','u','ha=',Ha
    do j=0,ny
      write(18,"(2(1X,E15.7E3))")y(0,j,0),uha(j)
    end do
    close(18)
    print*,' << profile_Ha.dat'
    !
  end subroutine incompact3d
  !
  subroutine bcchan(var)
    !
    real(8),intent(inout) :: var(:,:,:)
    !
    integer :: imax,jmax,kmax
    !
    imax=size(var,1)
    jmax=size(var,2)
    kmax=size(var,3)
    !
    var(1,:,2:kmax)=var(imax,:,2:kmax)
    var(:,:,1)   =var(:,:,kmax)
    !
  end subroutine bcchan
  !
  subroutine bchomo(var)
    !
    real(8),intent(inout) :: var(:,:,:)
    !
    integer :: imax,jmax,kmax
    !
    imax=size(var,1)
    jmax=size(var,2)
    kmax=size(var,3)
    !
    var(1,:,:)=var(imax,:,:)
    var(:,1,:)=var(:,jmax,:)
    var(:,:,1)=var(:,:,kmax)

  end subroutine bchomo
  !
end module userdefine
!+---------------------------------------------------------------------+
!| The end of the module userdefine.                                   |
!+---------------------------------------------------------------------+
