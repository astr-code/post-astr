module WingDesign
  !
  implicit none
  !
  type :: triangle
    real(8) :: a(3),b(3),c(3),normdir(3),area,cen(3)
    !+-------------------+---------------------------------------------+
    !|            a,b,c  | coordinates of the 3 vertex of the triangle.|
    !|           normdir | normal direction of the face.               |
    !|              area | area of of the triangle.                    |
    !+-------------------+---------------------------------------------+
  end type triangle
  !
  type :: lsegment
    real(8) :: a(2),b(2),normdir(2),length,cen(2)
    !+-------------------+---------------------------------------------+
    !|            a,b,c  | coordinates of the 3 vertex of the triangle.|
    !|           normdir | normal direction of the face.               |
    !|              area | area of of the triangle.                    |
    !+-------------------+---------------------------------------------+
  end type lsegment
  !
  type :: solid
    !
    character(len=32) :: name
    real(8) :: xmin(3),xmax(3),xref(3),xcen(3)
    integer :: num_face,num_edge
    type(triangle),allocatable :: face(:)
    type(lsegment),allocatable :: edge(:)
    !
  end type solid
  !
  contains
  !
  subroutine wing2d_design
    !
    use WriteTec
    !
    integer,parameter :: im=120,km=32
    real(8) :: x(0:im,0:km),y(0:im,0:km),z(0:im,0:km)
    integer :: i,k
    open(12,file='geom.txt')
      read(12,'(//)')
      do i=0,im/2
        read(12,*)x(i,0),y(i,0)
      enddo
      do i=im,im/2,-1
        read(12,*)x(i,0),y(i,0)
      enddo
    close(12)
    print*,' >> geom.txt'
    do i=0,im
      print*,i,x(i,0),y(i,0)
    enddo
    !
    do k=0,km
    do i=0,im
      x(i,k)=x(i,0)
      y(i,k)=y(i,0)
      z(i,k)=0.2d0/dble(km)*dble(k)
    enddo
    enddo
    !
    call writetecbin('tecwing.plt',x,'x',y,'y',z,'z',im,km)
    !
  end subroutine wing2d_design

  !+-------------------------------------------------------------------+
  !| This function is to design a delta wing STL file.                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-Jun-2022: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine delta_wing_design
    !
    use commvardefine, only: pi
    use WriteTec
    use LinearAlgegra
    !
    integer,parameter :: m=31,sp=42,ncir=469,nspa=41
    ! real(8) :: x(m),xl(m),yl(m),xu(m),yu(m),yc(m)
    ! real(8) :: xfoil(1:2*m-1),yfoil(1:2*m-1)
    ! real(8) :: xdel(2*m-1,sp),ydel(2*m-1,sp),zdel(2*m-1,sp)
    real(8) :: x(m,sp),y(m,sp),z(m,sp)
    real(8) :: xdel(ncir,nspa),ydel(ncir,nspa),zdel(ncir,nspa)
    !
    integer :: i,k,nline
    real(8) :: scal,rad,norm1(3),var1,xc,yc,zc
    !
    integer :: nface,jface
    type(solid) :: stl_solid(1)
    !
    ! do i=1,m
    !   x(i)=1.d0/dble(m-1)*dble(i-1)
    ! enddo
    ! !
    ! call naca5digit(x,yc,xu,yu,xl,yl,'64004')
    ! !
    ! do i=1,m
    !   xfoil(i)=xu(i)
    !   yfoil(i)=yu(i)
    ! enddo
    ! !
    ! do i=m+1,2*m-1
    !   xfoil(i)=xl(2*m-i)
    !   yfoil(i)=yl(2*m-i)
    ! enddo
    ! !
    ! open(18,file='naca64004.dat')
    ! do i=1,2*m-1
    ! write(18,*)xfoil(i),yfoil(i)
    ! end do
    ! close(18)
    ! print*,' <<< naca64004.dat ... done.'
    ! !
    ! rad=65.d0/180.d0*pi
    ! !
    ! do k=1,sp
    !   z=2.d0*1.50115d0/dble(sp-1)*dble(k-1)-1.50115d0
    !   scal=abs(1.d0-abs(z)/tan(rad))
    !   xdel(:,k)=xfoil*scal+(1.d0-scal)
    !   ydel(:,k)=yfoil*scal
    !   zdel(:,k)=z
    !   !
    !   print*,k,z,scal
    !   !
    ! enddo
    ! open(12,file='grorig.dat')
    ! nline=0
    ! do k=1,sp
    !   do i=1,m
    !     read(12,*)x(i,k),y(i,k),z(i,k)
    !     nline=nline+1
    !   enddo
    !   read(12,*)
    !   nline=nline+1
    ! enddo
    ! close(12)
    ! print*,' >> grorig.dat'
    ! !
    ! do k=1,21
    !   !
    !   ! upper surface
    !   do i=1,m
    !     xdel(i,k+20)=x(i,k)
    !     ydel(i,k+20)=y(i,k)
    !     zdel(i,k+20)=z(i,k)
    !   enddo
    !   ! bottom surface
    !   do i=m+1,ncir
    !     xdel(i,k+20)=x(2*m-i,k+21)
    !     ydel(i,k+20)=y(2*m-i,k+21)
    !     zdel(i,k+20)=z(2*m-i,k+21)
    !   enddo
    !   !
    ! enddo
    ! !
    ! do k=1,20
    !   ! symmetrical half
    !   xdel(:,k)= xdel(:,nspa-(k-1))
    !   ydel(:,k)=-ydel(:,nspa-(k-1))
    !   zdel(:,k)= zdel(:,nspa-(k-1))
    !   !
    ! enddo
    !
    open(12,file='geom.dat')
      read(12,'(//)')
      do i=1,ncir
        read(12,*)xdel(i,1),ydel(i,1)
      enddo
      ! read(12,*)
      ! do i=ncir,ncir/2+1,-1
      !   read(12,*)xdel(i,1),ydel(i,1)
      !   print*,i,xdel(i,1),ydel(i,1)
      ! enddo
    close(12)
    print*,' >> geom.dat'
    !
    do k=1,nspa
    do i=1,ncir
      xdel(i,k)=xdel(i,1)
      ydel(i,k)=ydel(i,1)
      zdel(i,k)=0.2d0/dble(nspa)*dble(k)
    enddo
    enddo
    call writetecbin('tec_wing_surface.plt',xdel,'x',ydel,'y',zdel,'z',ncir-1,nspa-1)
    !
    ! convert to a STL file
    nface=(ncir-1)*(nspa-1)*2+(ncir-1)*2
    !
    print*,' ** number of faces:',nface
    !
    stl_solid(1)%name='double_delta_wing'
    allocate(stl_solid(1)%face(1:nface))
    !
    jface=0
    !
    k=1
    !
    xc=0.d0
    yc=0.d0
    zc=0.d0
    do i=1,ncir-1
      xc=xc+xdel(i,k)
      yc=yc+ydel(i,k)
      zc=zc+zdel(i,k)
    enddo
    xc=xc/dble(ncir-1)
    yc=yc/dble(ncir-1)
    zc=zc/dble(ncir-1)
    do i=1,ncir-1
      jface=jface+1
      stl_solid(1)%face(jface)%a(1)=xdel(i,k)
      stl_solid(1)%face(jface)%a(2)=ydel(i,k)
      stl_solid(1)%face(jface)%a(3)=zdel(i,k)
      !
      stl_solid(1)%face(jface)%b(1)=xdel(i+1,k)
      stl_solid(1)%face(jface)%b(2)=ydel(i+1,k)
      stl_solid(1)%face(jface)%b(3)=zdel(i+1,k)
      !
      stl_solid(1)%face(jface)%c(1)=xc
      stl_solid(1)%face(jface)%c(2)=yc
      stl_solid(1)%face(jface)%c(3)=zc
      !
      norm1=cross_product(stl_solid(1)%face(jface)%a-stl_solid(1)%face(jface)%b,         &
                          stl_solid(1)%face(jface)%c-stl_solid(1)%face(jface)%b )
      !
      var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
      !
      stl_solid(1)%face(jface)%normdir=norm1/var1
      !
    enddo
    !
    k=nspa
    !
    xc=0.d0
    yc=0.d0
    zc=0.d0
    do i=1,ncir-1
      xc=xc+xdel(i,k)
      yc=yc+ydel(i,k)
      zc=zc+zdel(i,k)
    enddo
    xc=xc/dble(ncir-1)
    yc=yc/dble(ncir-1)
    zc=zc/dble(ncir-1)
    do i=1,ncir-1
      jface=jface+1
      stl_solid(1)%face(jface)%b(1)=xdel(i,k)
      stl_solid(1)%face(jface)%b(2)=ydel(i,k)
      stl_solid(1)%face(jface)%b(3)=zdel(i,k)
      !
      stl_solid(1)%face(jface)%a(1)=xdel(i+1,k)
      stl_solid(1)%face(jface)%a(2)=ydel(i+1,k)
      stl_solid(1)%face(jface)%a(3)=zdel(i+1,k)
      !
      stl_solid(1)%face(jface)%c(1)=xc
      stl_solid(1)%face(jface)%c(2)=yc
      stl_solid(1)%face(jface)%c(3)=zc
      !
      norm1=cross_product(stl_solid(1)%face(jface)%a-stl_solid(1)%face(jface)%b,         &
                          stl_solid(1)%face(jface)%c-stl_solid(1)%face(jface)%b )
      !
      var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
      !
      stl_solid(1)%face(jface)%normdir=norm1/var1
      !
    enddo
    !
    do k=1,nspa-1
      do i=1,ncir-1
        !
        jface=jface+1
        !
        stl_solid(1)%face(jface)%a(1)=xdel(i,k)
        stl_solid(1)%face(jface)%a(2)=ydel(i,k)
        stl_solid(1)%face(jface)%a(3)=zdel(i,k)
        !
        stl_solid(1)%face(jface)%b(1)=xdel(i,k+1)
        stl_solid(1)%face(jface)%b(2)=ydel(i,k+1)
        stl_solid(1)%face(jface)%b(3)=zdel(i,k+1)
        !
        stl_solid(1)%face(jface)%c(1)=xdel(i+1,k)
        stl_solid(1)%face(jface)%c(2)=ydel(i+1,k)
        stl_solid(1)%face(jface)%c(3)=zdel(i+1,k)
        !
        norm1=cross_product(stl_solid(1)%face(jface)%a-stl_solid(1)%face(jface)%b,         &
                            stl_solid(1)%face(jface)%c-stl_solid(1)%face(jface)%b )
        !
        var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
        !
        stl_solid(1)%face(jface)%normdir=norm1/var1
        !
        jface=jface+1
        !
        stl_solid(1)%face(jface)%a(1)=xdel(i,k+1)
        stl_solid(1)%face(jface)%a(2)=ydel(i,k+1)
        stl_solid(1)%face(jface)%a(3)=zdel(i,k+1)
        !
        stl_solid(1)%face(jface)%b(1)=xdel(i+1,k+1)
        stl_solid(1)%face(jface)%b(2)=ydel(i+1,k+1)
        stl_solid(1)%face(jface)%b(3)=zdel(i+1,k+1)
        !
        stl_solid(1)%face(jface)%c(1)=xdel(i+1,k)
        stl_solid(1)%face(jface)%c(2)=ydel(i+1,k)
        stl_solid(1)%face(jface)%c(3)=zdel(i+1,k)
        !
        norm1=cross_product(stl_solid(1)%face(jface)%a-stl_solid(1)%face(jface)%b,         &
                            stl_solid(1)%face(jface)%c-stl_solid(1)%face(jface)%b )
        !
        var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
        !
        stl_solid(1)%face(jface)%normdir=norm1/var1
        !
      enddo
    enddo
    stl_solid(1)%num_face=jface
    !
    print*,' ** jface',jface
    !
    call stla_write('double_delta_wing.stl',stl_solid)
    !
    call tecsolid('tecsolid.plt',stl_solid)
    ! open(18,file='naca64004_chord.dat')
    ! do i=1,m
    ! write(18,*)x(i),yc(i)
    ! end do
    ! close(18)
    ! print*,' <<< naca64004_chord.dat ... done.'
    ! !
    ! open(18,file='naca64004_upper.dat')
    ! do i=1,m
    ! write(18,*)xu(i),yu(i)
    ! end do
    ! close(18)
    ! print*,' <<< naca64004_upper.dat ... done.'
    ! !
    ! open(18,file='naca64004_lower.dat')
    ! do i=1,m
    ! write(18,*)xl(i),yl(i)
    ! end do
    ! close(18)
    ! print*,' <<< naca64004_lower.dat ... done.'
    !
  end subroutine delta_wing_design
  !+-------------------------------------------------------------------+
  !| The end of the subroutine delta_wing_design.                      |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write a STL file.                      !
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 02-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine stla_write(out_file_name,asolid)
    !
    ! argument
    character(len=*),intent(in) :: out_file_name
    type(solid),intent(in) :: asolid(:)
    !
    ! local data
    integer :: iunit,ios,ierror,solidstate,lchar,i
    integer :: num_solid,nsolid,nface
    character(len=255) word1,word2,word3,text
    logical :: done
    real(8) :: dval
    !
    num_solid=size(asolid)
    !
    iunit=16
    !
    open(unit=iunit,file=out_file_name,iostat=ios)
    do nsolid=1,num_solid
      write(iunit,'(A,I0)')'solid ',nsolid
      do nface=1,asolid(nsolid)%num_face
        write(iunit,'(A,3(1X,E15.7E3))')'  facet normal ',             &
                                   asolid(nsolid)%face(nface)%normdir(:)
        write(iunit,'(A)')'    outer loop'
        write(iunit,'(A,3(1X,E15.7E3))')'      vertex ',              &
                                         asolid(nsolid)%face(nface)%a(:)
        write(iunit,'(A,3(1X,E15.7E3))')'      vertex ',              &
                                         asolid(nsolid)%face(nface)%b(:)
        write(iunit,'(A,3(1X,E15.7E3))')'      vertex ',              &
                                         asolid(nsolid)%face(nface)%c(:)
        write(iunit,'(A)')'    endloop'
        write(iunit,'(A)')'  endfacet'
      enddo
      write(iunit,'(A,I0)')'endsolid ',nsolid
    enddo
    close(iunit)
    print*,' << ',out_file_name
    !
    return
    !
  end subroutine stla_write
  !+-------------------------------------------------------------------+
  !| The end of the function stla_write.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is to generate a  4-digit NACA airfoil.             |
  !+-------------------------------------------------------------------+
  !| ref: https://en.wikipedia.org/wiki/NACA_airfoil                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine naca4digit(xin,x,y,name,surface)
    !
    ! arguments
    real(8),intent(in) :: xin
    real(8),intent(out) :: x,y
    character(len=4),intent(in) :: name
    character(len=5),intent(in) :: surface
    !
    ! local data 
    real(8) :: yt,t,yc,theter,mr,pr
    integer :: m,p
    !
    read(name(1:1),*)m
    read(name(2:2),*)p
    read(name(3:4),*)t
    !
    t=t/100.d0
    !
    yt=5.d0*t*(0.2969d0*sqrt(xin)-0.1260d0*xin-0.3516d0*xin*xin+       &
               0.2843d0*xin**3-0.1036d0*xin**4)
    !
    if(p==0 .and. m==0) then
      ! NACA-00XX, a symmetrical 4-digit NACA airfoil
      !
      if(surface=='upper') then
        x=xin
        y=yt
      elseif(surface=='lower') then
        x=xin
        y=-yt
      else
        print*,' surface=',surface
        stop '!! ERROR 1 @ naca4digit'
      endif
      !
    else
      !
      mr=0.01d0*dble(m)
      pr=0.1d0*dble(p)
      !
      if(xin>=0.d0 .and. xin<=pr) then
        yc=mr/pr/pr*(2.d0*pr*xin-xin*xin)
        theter=atan(2.d0*mr/pr/pr*(pr-xin))
      elseif(xin>=pr .and. xin<=1.d0) then
        yc=mr/(1-pr)**2*(1.d0-2.d0*pr+2.d0*pr*xin-xin**2)
        theter=atan(2.d0*mr/(1-pr)**2*(pr-xin))
      else
        stop '!! ERROR 2 @ naca4digit'
      endif
      !
      if(surface=='upper') then
        x=xin-yt*sin(theter)
        y= yc+yt*cos(theter)
      elseif(surface=='lower') then
        x=xin+yt*sin(theter)
        y= yc-yt*cos(theter)
      else
        print*,' surface=',surface
        stop '!! ERROR 3 @ naca4digit'
      endif
      !
    endif
    !
  end subroutine naca4digit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine naca4digit.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is to generate a  4-digit NACA airfoil.             |
  !+-------------------------------------------------------------------+
  !| ref: https://en.wikipedia.org/wiki/NACA_airfoil                   |
  !| http://airfoiltools.com/airfoil/naca5digit                        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine naca5digit(xin,yc,x_upper,y_upper,x_lower,y_lower,name)
    !
    ! arguments
    real(8),intent(in) :: xin(:)
    character(len=5),intent(in) :: name
    real(8),intent(out) :: yc(1:size(xin))
    real(8),intent(out) :: x_upper(1:size(xin)),y_upper(1:size(xin))
    real(8),intent(out) :: x_lower(1:size(xin)),y_lower(1:size(xin))
    !
    real(8) :: dyc(1:size(xin))
    !
    ! local data 
    real(8) :: xx,r,k1,k2,theter,yt,t
    integer :: l,p,q,im,i
    !
    read(name(1:1),*)l
    read(name(2:2),*)p
    read(name(3:3),*)q
    read(name(4:5),*)xx
    !
    !Digits Letter  Example Description
    !     1      L        2 This digit controls the camber. It indicates the designed coefficient of lift (Cl) multiplied by 3/20. In the examble L=2 so Cl=0.3
    !     2      P        3 The position of maximum camber divided by 20. In the examble P=3 so maximum camber is at 0.15 or 15% chord
    !     3      Q        0 0 = normal camber line, 1 = reflex camber line
    !   4-5     XX       12 The maximum thickness as percentage.In the examble XX=12 so the maximum thickness is 0.12 or 12% chord.
    write(*,'(A,A)')'NACA-',name
    write(*,'(A,I0)')'L: ',l
    write(*,'(A,I0)')'P: ',p
    write(*,'(A,I0)')'Q: ',Q
    write(*,'(A,F6.2)')'XX: ',xx
    !
    if(p==1 .and. q==0) then
      r=0.0580; k1=361.400 
    elseif(p==2 .and. q==0) then
      r=0.1260; k1=51.640  
    elseif(p==3 .and. q==0) then
      r=0.2025; k1=15.957  
    elseif(p==4 .and. q==0) then
      r=0.2900; k1=6.643   
    elseif(p==5 .and. q==0) then
      r=0.3910; k1=3.230   
    elseif(p==2 .and. q==1) then
      r=0.1300; k1=51.990;  k2=0.000764
    elseif(p==3 .and. q==1) then
      r=0.2170; k1=15.793;  k2=0.00677
    elseif(p==4 .and. q==1) then
      r=0.3180; k1=6.520 ;  k2=0.0303
    elseif(p==5 .and. q==1) then
      r=0.4410; k1=3.191 ;  k2=0.1355
    else
      print*,' digit out of table @ naca5digit'
      stop
    endif
    !
    t=xx/100.d0
    !
    im=size(xin)
    do i=1,im
      !
      ! print*,i,xin(i)
      if(xin(i)<r) then
         yc(i)=k1/6.d0*(xin(i)**3-3.d0*r*xin(i)**2+r**2*(3.d0-r)*xin(i))
        dyc(i)=k1/6.d0*(3.d0*xin(i)**2-6.d0*r*xin(i)+r**2*(3.d0-r))
      else
         yc(i)= k1*r**3/6.d0*(1.d0-xin(i))
        dyc(i)=-k1*r**3/6.d0
      endif
      !
      yt=5.d0*t*(0.2969d0*sqrt(xin(i))-0.1260d0*xin(i)-0.3516d0*xin(i)*xin(i)+       &
               0.2843d0*xin(i)**3-0.1036d0*xin(i)**4)
      !
      theter=atan(dyc(i))
      !
      x_upper(i)=xin(i)-yt*sin(theter)
      y_upper(i)=yc(i) +yt*cos(theter)
      !
      x_lower(i)=xin(i)+yt*sin(theter)
      y_lower(i)=yc(i) -yt*cos(theter)
      !
    enddo
    !
    ! tt=tt/100.d0
    ! !
    ! yt=5.d0*t*(0.2969d0*sqrt(xin)-0.1260d0*xin-0.3516d0*xin*xin+       &
    !            0.2843d0*xin**3-0.1036d0*xin**4)
    ! !
    ! if(p==0 .and. m==0) then
    !   ! The NACA five-digit series describes more complex airfoil shapes.
    !   !
    !   if(surface=='upper') then
    !     x=xin
    !     y=yt
    !   elseif(surface=='lower') then
    !     x=xin
    !     y=-yt
    !   else
    !     print*,' surface=',surface
    !     stop '!! ERROR 1 @ naca4digit'
    !   endif
    !   !
    ! else
    !   !
    !   if()
    !   mr=0.01d0*dble(m)
    !   pr=0.1d0*dble(p)
    !   !
    !   if(xin>=0.d0 .and. xin<=pr) then
    !     yc=mr/pr/pr*(2.d0*pr*xin-xin*xin)
    !     theter=atan(2.d0*mr/pr/pr*(pr-xin))
    !   elseif(xin>=pr .and. xin<=1.d0) then
    !     yc=mr/(1-pr)**2*(1.d0-2.d0*pr+2.d0*pr*xin-xin**2)
    !     theter=atan(2.d0*mr/(1-pr)**2*(pr-xin))
    !   else
    !     stop '!! ERROR 2 @ naca4digit'
    !   endif
    !   !
    !   if(surface=='upper') then
    !     x=xin-yt*sin(theter)
    !     y= yc+yt*cos(theter)
    !   elseif(surface=='lower') then
    !     x=xin+yt*sin(theter)
    !     y= yc-yt*cos(theter)
    !   else
    !     print*,' surface=',surface
    !     stop '!! ERROR 3 @ naca4digit'
    !   endif
    !   !
    ! endif
    !
  end subroutine naca5digit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine naca5digit.                             |
  !+-------------------------------------------------------------------+
  !!+-------------------------------------------------------------------+
  !| This subroutine is used to FE-Volume Tetrahedral Data.            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine tecsolid(filename,asolid,dim)
    !
    type(solid),intent(in) :: asolid(:)
    character(len=*),intent(in) :: filename
    integer,intent(in),optional :: dim
    !
    integer :: nsolid,jsolid,j,ndim,unitf
    !
    nsolid=size(asolid)
    !
    if(present(dim)) then
      ndim=dim
    else
      ndim=3
    endif
    !
    unitf=18
    !
    if(ndim==2) then
      !
      open(unitf,file=filename,form='formatted')
      write(unitf,'(A)')'TITLE = "FE-Volume TRIANGLE Data"'
      do jsolid=1,nsolid
        write(unitf,'(a)')'variables = "x", "y", "nx", "ny"'
        write(unitf,'(3(A,I0),A)')'ZONE T="P_',jsolid,                    &
                          '", DATAPACKING=BLOCK, NODES=',              &
                         asolid(jsolid)%num_edge*2,', ELEMENTS=',      &
                         asolid(jsolid)%num_edge,', ZONETYPE=FELINESEG'
        write(unitf,'(A)')'VarLocation=([3-4]=CellCentered)'
        write(unitf,'(2(1X,E15.7E3))')(asolid(jsolid)%edge(j)%a(1),       &
                                    asolid(jsolid)%edge(j)%b(1),       &
                                    j=1,asolid(jsolid)%num_edge)
        write(unitf,'(2(1X,E15.7E3))')(asolid(jsolid)%edge(j)%a(2),       &
                                    asolid(jsolid)%edge(j)%b(2),       &
                                    j=1,asolid(jsolid)%num_edge)
        write(unitf,'(8(1X,E15.7E3))')(asolid(jsolid)%edge(j)%normdir(1), &
                                            j=1,asolid(jsolid)%num_edge)
        write(unitf,'(8(1X,E15.7E3))')(asolid(jsolid)%edge(j)%normdir(2), &
                                            j=1,asolid(jsolid)%num_edge)
        write(unitf,'(2(1X,I8))')(j,j=1,asolid(jsolid)%num_edge*2)
        !
      enddo
      !
    elseif(ndim==3) then
      open(unitf,file=filename,form='formatted')
      write(unitf,'(A)')'TITLE = "FE-Volume TRIANGLE Data"'
      do jsolid=1,nsolid
        write(unitf,'(a)')'variables = "x","y","z","nx","ny","nz"'
        write(unitf,'(3(A,I0),A)')'ZONE T="P_',jsolid,                    &
                          '", DATAPACKING=BLOCK, NODES=',              &
                         asolid(jsolid)%num_face*3,', ELEMENTS=',      &
                         asolid(jsolid)%num_face,', ZONETYPE=FETRIANGLE'
        write(unitf,'(A)')'VarLocation=([4-6]=CellCentered)'
        write(unitf,'(3(1X,E15.7E3))')(asolid(jsolid)%face(j)%a(1),       &
                                    asolid(jsolid)%face(j)%b(1),       &
                                    asolid(jsolid)%face(j)%c(1),       &
                                    j=1,asolid(jsolid)%num_face)
        write(unitf,'(3(1X,E15.7E3))')(asolid(jsolid)%face(j)%a(2),       &
                                    asolid(jsolid)%face(j)%b(2),       &
                                    asolid(jsolid)%face(j)%c(2),       &
                                    j=1,asolid(jsolid)%num_face)
        write(unitf,'(3(1X,E15.7E3))')(asolid(jsolid)%face(j)%a(3),       &
                                    asolid(jsolid)%face(j)%b(3),       &
                                    asolid(jsolid)%face(j)%c(3),       &
                                    j=1,asolid(jsolid)%num_face)
        ! write(unitf,'(8(1X,E15.7E3))')(num1d3*(asolid(jsolid)%face(j)%a(1) + &
        !                                        asolid(jsolid)%face(j)%b(1) + &
        !                                        asolid(jsolid)%face(j)%c(1)),j=1,asolid(jsolid)%num_face)
        ! write(unitf,'(8(1X,E15.7E3))')(num1d3*(asolid(jsolid)%face(j)%a(2) + &
        !                                        asolid(jsolid)%face(j)%b(2) + &
        !                                        asolid(jsolid)%face(j)%c(2)),j=1,asolid(jsolid)%num_face)
        ! write(unitf,'(8(1X,E15.7E3))')(num1d3*(asolid(jsolid)%face(j)%a(3) + &
        !                                        asolid(jsolid)%face(j)%b(3) + &
        !                                        asolid(jsolid)%face(j)%c(3)),j=1,asolid(jsolid)%num_face)
        write(unitf,'(8(1X,E15.7E3))')(asolid(jsolid)%face(j)%normdir(1), &
                                            j=1,asolid(jsolid)%num_face)
        write(unitf,'(8(1X,E15.7E3))')(asolid(jsolid)%face(j)%normdir(2), &
                                            j=1,asolid(jsolid)%num_face)
        write(unitf,'(8(1X,E15.7E3))')(asolid(jsolid)%face(j)%normdir(3), &
                                            j=1,asolid(jsolid)%num_face)
        write(unitf,'(3(1X,I8))')(j,j=1,asolid(jsolid)%num_face*3)
        !
      enddo
      !
    endif
    !
    close(unitf)
    print*,' << ',filename
    !
    return
    !
  end subroutine tecsolid
  !+-------------------------------------------------------------------+
  !| The end of the subroutine tecsolid.                               |
  !+-------------------------------------------------------------------+
    !
end module WingDesign