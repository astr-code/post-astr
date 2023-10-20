!+---------------------------------------------------------------------+
!| This module contains subroutines used for interpolations.           |
!+---------------------------------------------------------------------+
module interpolation
  !
  implicit none
  !
  interface binterp
    module procedure binterp_2dmesh_1dmesh
    module procedure binterp_2dmesh_1dmesh_homz
  end interface binterp
  !
  interface interpolate
    module procedure qualinear_mesh2point_2d_3v_kspan
    module procedure qualinear_mesh2point_2d_4v_kspan
  end interface interpolate
  interface linear1d
    module procedure linear1d_s
    module procedure linear1d_a1
    module procedure linear1d_a2
  end interface linear1d
  !
  interface QuadraticSpline
    module procedure QuadraticSpline1d1v
  end interface QuadraticSpline
  !
  Interface regularlinearinterp
    !
    module procedure regularlinearinterp1d1v
    module procedure regularlinearinterp1d2v
    module procedure regularlinearinterp1d3v
    module procedure regularlinearinterp1d4v
    module procedure regularlinearinterp1d5v
    module procedure regularlinearinterp1d6v
    module procedure regularlinearinterp1d7v
    module procedure regularlinearinterp1d8v
    module procedure regularlinearinterp1d9v
    !
    module procedure regularlinearinterp2d1v
    module procedure regularlinearinterp2d2v
    module procedure regularlinearinterp2d3v
    module procedure regularlinearinterp2d4v
    module procedure regularlinearinterp2d5v
    module procedure regularlinearinterp2d5v_homz
    module procedure regularlinearinterp2d6v
    module procedure regularlinearinterp2d6v_homz
    module procedure regularlinearinterp2d7v
    !
    module procedure regularlinearinterp3d1v
    module procedure regularlinearinterp3d2v
    module procedure regularlinearinterp3d3v
    module procedure regularlinearinterp3d4v
    module procedure regularlinearinterp3d5v
    module procedure regularlinearinterp3d6v
    module procedure regularlinearinterp3d6v_straight
    module procedure regularlinearinterp3d6v_sp
    module procedure regularlinearinterp3d7v_sp
    !
  end Interface regularlinearinterp
  !
  interface InvdisInterp2D
    !
    module procedure InvdisInterp2D1
    module procedure InvdisInterp2D2
    module procedure InvdisInterp2D3
    module procedure InvdisInterp2D4
    module procedure InvdisInterp2D4_1_km
    module procedure InvdisInterp2D5
    module procedure InvdisInterp2D6
    !module procedure InvdisInterp2D7
    !module procedure InvdisInterp2D8
    !
  end interface InvdisInterp2D
  !
  interface InvdisInterp3D
    !
    module procedure InvdisInterp3D6
    !
  end interface InvdisInterp3D
  !
  contains
  !
  subroutine qualinear_mesh2point_2d_3v_kspan(x_src,y_src,v1_src,v2_src,v3_src, &
                                              x_dst,y_dst,v1_dst,v2_dst,v3_dst)
    !
    real(8),intent(in) :: x_src(:,:),y_src(:,:),x_dst,y_dst
    real(8),intent(in) :: v1_src(:,:,:),v2_src(:,:,:),v3_src(:,:,:)
    real(8),intent(out) :: v1_dst(:),v2_dst(:),v3_dst(:)
    !
    integer :: im,jm,km,i,j,k,is,js
    logical :: lin
    !
    real(8) :: var1,dismin
    integer :: perc
    real(8) :: wei(0:1,0:1),dis(0:1,0:1),temp(5)
    !
    im=size(v1_src,1)
    jm=size(v1_src,2)
    km=size(v1_src,3)
    print*,' ** dimension:',im,jm,km
    !
    lin=.false.
    looj: do j=1,jm-1
    looi: do i=1,im-1
      !
      lin=pointintriangle(x_src(i,j),x_src(i+1,j),x_src(i+1,j+1), &
                          y_src(i,j),y_src(i+1,j),y_src(i+1,j+1), &
                          x_dst,y_dst   )
      !
      if(lin) then
        is=i
        js=j
        exit looj
      endif
      !
      lin=pointintriangle(x_src(i+1,j+1),x_src(i,j+1),x_src(i,j), &
                          y_src(i+1,j+1),y_src(i,j+1),y_src(i,j), &
                          x_dst,y_dst   )
      !
      if(lin) then
        is=i
        js=j
        exit looj
      endif
      !
    enddo looi
    enddo looj
    !
    print*,' ** is,js:',is,js
    !
    var1=0.d0
    do j=0,1
    do i=0,1
      !
      dis(i,j)=(x_src(is+i,js+j)-x_dst)**2+                 &
               (y_src(is+i,js+j)-y_dst)**2
      !
      if(dis(i,j)<=1.d-16) then
        !
        v1_dst(1:km)=v1_src(is+i,js+j,1:km)
        v2_dst(1:km)=v2_src(is+i,js+j,1:km)
        v3_dst(1:km)=v3_src(is+i,js+j,1:km)
        !
        return
        !
      else
        wei(i,j)=1.d0/dis(i,j)
        var1=var1+wei(i,j)
      endif
      !
    enddo
    enddo
    !
    wei=wei/var1
    !
    v1_dst=0.d0
    v2_dst=0.d0
    v3_dst=0.d0
    do j=0,1
    do i=0,1
      !
      v1_dst(1:km)=v1_dst(1:km)+v1_src(is+i,js+j,1:km)*wei(i,j)
      v2_dst(1:km)=v2_dst(1:km)+v2_src(is+i,js+j,1:km)*wei(i,j)
      v3_dst(1:km)=v3_dst(1:km)+v3_src(is+i,js+j,1:km)*wei(i,j)
      !
    enddo
    enddo
    !
  end subroutine qualinear_mesh2point_2d_3v_kspan
  subroutine qualinear_mesh2point_2d_4v_kspan(x_src,y_src,v1_src,v2_src,v3_src, &
                                              v4_src, &
                                              x_dst,y_dst,v1_dst,v2_dst,v3_dst, &
                                              v4_dst)
    !
    real(8),intent(in) :: x_src(:,:),y_src(:,:),x_dst,y_dst
    real(8),intent(in) :: v1_src(:,:,:),v2_src(:,:,:),v3_src(:,:,:),v4_src(:,:,:)
    real(8),intent(out) :: v1_dst(:),v2_dst(:),v3_dst(:),v4_dst(:)
    !
    integer :: im,jm,km,i,j,k,is,js
    logical :: lin
    !
    real(8) :: var1,dismin
    integer :: perc
    real(8) :: wei(0:1,0:1),dis(0:1,0:1),temp(5)
    !
    im=size(v1_src,1)
    jm=size(v1_src,2)
    km=size(v1_src,3)
    print*,' ** dimension:',im,jm,km
    !
    lin=.false.
    is=-1
    js=-1
    looj: do j=1,jm-1
    looi: do i=1,im-1
      !
      lin=pointintriangle(x_src(i,j),x_src(i+1,j),x_src(i+1,j+1), &
                          y_src(i,j),y_src(i+1,j),y_src(i+1,j+1), &
                          x_dst,y_dst   )
      !
      if(lin) then
        is=i
        js=j
        exit looj
      endif
      !
      lin=pointintriangle(x_src(i+1,j+1),x_src(i,j+1),x_src(i,j), &
                          y_src(i+1,j+1),y_src(i,j+1),y_src(i,j), &
                          x_dst,y_dst   )
      !
      if(lin) then
        is=i
        js=j
        exit looj
      endif
      !
    enddo looi
    enddo looj
    !
    print*,' ** is,js:',is,js
    print*,' ** x_dst,y_dst:',x_dst,y_dst
    print*,' ** x_src,y_src:',x_src(is,js),y_src(is,js)
    !
    var1=0.d0
    do j=0,1
    do i=0,1
      !
      dis(i,j)=(x_src(is+i,js+j)-x_dst)**2+                 &
               (y_src(is+i,js+j)-y_dst)**2
      !
      if(dis(i,j)<=1.d-16) then
        !
        v1_dst(1:km)=v1_src(is+i,js+j,1:km)
        v2_dst(1:km)=v2_src(is+i,js+j,1:km)
        v3_dst(1:km)=v3_src(is+i,js+j,1:km)
        v4_dst(1:km)=v4_src(is+i,js+j,1:km)
        !
        return
        !
      else
        wei(i,j)=1.d0/dis(i,j)
        var1=var1+wei(i,j)
      endif
      !
    enddo
    enddo
    !
    wei=wei/var1
    !
    v1_dst=0.d0
    v2_dst=0.d0
    v3_dst=0.d0
    do j=0,1
    do i=0,1
      !
      v1_dst(1:km)=v1_dst(1:km)+v1_src(is+i,js+j,1:km)*wei(i,j)
      v2_dst(1:km)=v2_dst(1:km)+v2_src(is+i,js+j,1:km)*wei(i,j)
      v3_dst(1:km)=v3_dst(1:km)+v3_src(is+i,js+j,1:km)*wei(i,j)
      v4_dst(1:km)=v4_dst(1:km)+v4_src(is+i,js+j,1:km)*wei(i,j)
      !
    enddo
    enddo
    !
  end subroutine qualinear_mesh2point_2d_4v_kspan
  !
  function pointintriangle(xa,xb,xc,ya,yb,yc,xp,yp) result(lin)
    !
    ! arguments
    real(8),intent(in) :: xa,xb,xc,ya,yb,yc,xp,yp
    logical :: lin
    !
    ! local data
    real(8) :: za,zb,zc,zp
    real(8) :: pa(3),pb(3),pc(3),t1(3),t2(3),t3(3),v0(3),v1(3),v2(3)
    real(8) :: a1,a2,a3,error,var1,var2,dot00,dot01,dot02,dot11,dot12, &
               inverDeno,u,v
    real(8) :: epsilon
    !
    v0(1) = xc - xa; v1(1) = xb - xa; v2(1) = xp - xa
    v0(2) = yc - ya; v1(2) = yb - ya; v2(2) = yp - ya
    v0(3) = 0.d0;    v1(3) = 0.d0;    v2(3) = 0.d0;
    ! v0(3) = zc - za; v1(3) = zb - za; v2(3) = zp - za

    dot00 = dot_product(v0,v0)
    dot01 = dot_product(v0,v1)
    dot02 = dot_product(v0,v2)
    dot11 = dot_product(v1,v1)
    dot12 = dot_product(v1,v2)

    inverDeno = 1.d0 / (dot00 * dot11 - dot01 * dot01) 
    
    u = (dot11 * dot02 - dot01 * dot12) * inverDeno

    if(isnan(u)) then
      print*,xa 
      print*,xb
      print*,xc
      stop ' !! ERROR @ pointintriangle_nodes'
    endif
    ! print*,' ** u=',u
    !
    if (u < 0.d0 .or. u > 1.d0) then
      ! if u out of range, return directly
      lin=.false.
      return
    endif
    !
    v = (dot00 * dot12 - dot01 * dot02) * inverDeno
    ! print*,' ** v=',v
    if (v < 0.d0 .or. v > 1.d0) then
      ! if v out of range, return directly
      lin=.false.
      return 
    endif
    !
    if(u + v <= 1.d0) then
      lin=.true.
    else
      lin=.false.
    endif
    !
    return
    !
  end function pointintriangle
  !
  subroutine QuadraticSpline1d1v(x1,v11,dim1,                      &
                                 x2,v21,dim2)
    !
    integer,intent(in) :: dim1,dim2
    !
    real(8),intent(in) :: x1(0:dim1),v11(0:dim1)
    real(8),intent(in) :: x2(0:dim2)
    real(8),intent(out) :: v21(0:dim2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    real(8) :: zcurv(0:dim1)
    !
    zcurv(0)=0.d0
    do i=1,dim1
      zcurv(i)=-zcurv(i-1)+2.d0*(v11(i)-v11(i-1))/(x1(i)-x1(i-1))
    enddo
    !
    do i=0,dim2
      do i1=1,dim1
        if( x2(i)>=x1(i1-1) .and. x2(i)<=x1(i1) ) then
          v21(i)=v11(i1-1)+zcurv(i1-1)*(x2(i)-x1(i1-1))+               &
                 0.5d0*(zcurv(i1)-zcurv(i1-1))/(x1(i1)-x1(i1-1))*(x2(i)-x1(i1-1))**2
        endif
      enddo
    enddo
    !
  end subroutine QuadraticSpline1d1v
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to interpolate 3D flowfield between 2     |
  !| regular meshes with linear approach.                              |
  !+-------------------------------------------------------------------+
  subroutine regularlinearinterp1d1v(x1,v11,dim1,                      &
                                     x2,v21,dim2,mode)
    !
    integer,intent(in) :: dim1,dim2
    !
    real(8),intent(in) :: x1(0:dim1),v11(0:dim1)
    real(8),intent(in) :: x2(0:dim2)
    character(len=3),optional :: mode
    real(8),intent(out) :: v21(0:dim2)
    real(8) :: grad0,grad1,segd0,grad,gradm,gradm1,segdm,damp
    !
    integer :: i,j,k,i1,j1,k1,nper,i0m,irem
    !
    i0m=-100
    irem=-100
    !
    ! if(present(mode)) print*, '** intepolation mode:' mode
    !
    grad0=(v11(1)-v11(0))/(x1(1)-x1(0))
    grad1=(v11(2)-v11(0))/(x1(2)-x1(0))
    segd0=(grad1-grad0)/(x1(1)-x1(0))
    !
    gradm=(v11(dim1)-v11(dim1-1))/(x1(dim1)-x1(dim1-1))
    gradm1=(v11(dim1)-v11(dim1-2))/(x1(dim1)-x1(dim1-2))
    segdm=(gradm-gradm1)/(x1(dim1)-x1(dim1-1))
    !
    nper=0
    do i=0,dim2
      !
    	do i1=1,dim1
        !
    		if( x2(i)>=x1(i1-1) .and. x2(i)<=x1(i1) ) then
          !
          v21(i)=linear1d(x1(i1-1),x1(i1),v11(i1-1),v11(i1),x2(i))
          !
          exit
          !
        elseif(x2(i)<x1(0)) then
          !
          if( .not. present(mode) .or. mode(1:1)=='|') exit
          !
          ! if(i0m<0) i0m=i ! first node
          !
          if(mode(1:1)=='-') then
            v21(i)=v11(0)
          elseif(mode(1:1)=='/') then
            v21(i)=v11(0)+(x2(i)-x1(0))*grad0
          elseif(mode(1:1)=='~') then
            grad=grad0+(x2(i)-x1(0))*segd0
            v21(i)=v11(0)+(x2(i)-x1(0))*grad
          elseif(mode(1:1)=='<') then
            damp=exp(-1000000.d0*(x2(i)-x1(0))**2)
            ! print*,x2(i),'<',damp
            v21(i)=v11(dim1)*damp
          else
            stop ' mode not defined @ regularlinearinterp1d1v'
          endif
          !
          exit
          !
        elseif(x2(i)>x1(dim1)) then
          !
          if( .not. present(mode) .or. mode(3:3)=='|') exit
          !
          if(irem<0) irem=i ! first node
          !
          if(mode(3:3)=='-') then
            v21(i)=v11(dim1)
          elseif(mode(3:3)=='/') then
            v21(i)=v11(dim1)+(x2(i)-x1(dim1))*gradm
          elseif(mode(3:3)=='~') then
            grad=gradm+(x2(i)-x1(dim1))*segdm
            v21(i)=v11(dim1)+(x2(i)-x1(dim1))*grad
          elseif(mode(3:3)=='>') then
            ! damp=exp(-0.1d0*dble(i-irem))
            damp=exp(-1000000.d0*(x2(i)-x1(dim1))**2)
            ! print*,x2(i),'>',damp
            v21(i)=v11(dim1)*damp
            ! grad=gradm*exp(-0.1d0*dble(i-irem))
            ! v21(i)=v11(dim1)+(x2(i)-x1(dim1))*grad
          else
            stop ' mode not defined @ regularlinearinterp1d1v'
          endif
          !
          exit
          !
          !
    		end if
        !
    	end do
      !
      nper=nper+1
      !
      ! write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
      !                     '  ** Interpolating ... ',100*nper/(dim2+1),' %'
    end do
    !
    ! print*
    !
    ! stop
    !
  end subroutine regularlinearinterp1d1v
  !
  subroutine regularlinearinterp1d2v(x1,v11,v12,dim1,                  &
                                     x2,v21,v22,dim2)
    !
    integer,intent(in) :: dim1,dim2
    !
    real(8),intent(in) :: x1(0:dim1),v11(0:dim1),v12(0:dim1)
    real(8),intent(in) :: x2(0:dim2)
    real(8),intent(out) :: v21(0:dim2),v22(0:dim2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    !
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i)>=x1(i1-1) .and. x2(i)<=x1(i1) ) then
    		  !
    		  v21(i)=linear1d(x1(i1-1),x1(i1),v11(i1-1),v11(i1),x2(i))
    		  v22(i)=linear1d(x1(i1-1),x1(i1),v12(i1-1),v12(i1),x2(i))
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	                    '  ** Interpolating ... ',100*nper/(dim2+1),' %'
    end do
    !
    print*
    !
  end subroutine regularlinearinterp1d2v
  !
  subroutine regularlinearinterp1d3v(x1,v11,v12,v13,dim1,              &
                                     x2,v21,v22,v23,dim2)
    !
    integer,intent(in) :: dim1,dim2
    !
    real(8),intent(in) :: x1(0:dim1),v11(0:dim1),v12(0:dim1),v13(0:dim1)
    real(8),intent(in) :: x2(0:dim2)
    real(8),intent(out) :: v21(0:dim2),v22(0:dim2),v23(0:dim2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    !
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i)>=x1(i1-1) .and. x2(i)<=x1(i1) ) then
    		  !
    		  v21(i)=linear1d(x1(i1-1),x1(i1),v11(i1-1),v11(i1),x2(i))
    		  v22(i)=linear1d(x1(i1-1),x1(i1),v12(i1-1),v12(i1),x2(i))
    		  v23(i)=linear1d(x1(i1-1),x1(i1),v13(i1-1),v13(i1),x2(i))
    		  !
    		  exit
    		  !
        elseif(x2(i)>=x1(dim1)) then
          !
          v21(i)=v11(dim1)
          v22(i)=v12(dim1)
          v23(i)=v13(dim1)
          !
          exit
          !
        elseif(x2(i)<=x1(0)) then
          !
          v21(i)=v11(0)
          v22(i)=v12(0)
          v23(i)=v13(0)
          !
          exit
          !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	! write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	!                     '  ** Interpolating ... ',100*nper/(dim2+1),' %'
    end do
    !
    print*,' '
    !
  end subroutine regularlinearinterp1d3v
  !
  subroutine regularlinearinterp1d4v(x1,v11,v12,v13,v14,dim1,          &
                                     x2,v21,v22,v23,v24,dim2)
    !
    integer,intent(in) :: dim1,dim2
    !
    real(8),intent(in) :: x1(0:dim1),v11(0:dim1),v12(0:dim1),          &
                          v13(0:dim1),v14(0:dim1)
    real(8),intent(in) :: x2(0:dim2)
    real(8),intent(out) :: v21(0:dim2),v22(0:dim2),v23(0:dim2),        &
                           v24(0:dim2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    !
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i)>=x1(i1-1) .and. x2(i)<=x1(i1) ) then
    		  !
    		  v21(i)=linear1d(x1(i1-1),x1(i1),v11(i1-1),v11(i1),x2(i))
    		  v22(i)=linear1d(x1(i1-1),x1(i1),v12(i1-1),v12(i1),x2(i))
    		  v23(i)=linear1d(x1(i1-1),x1(i1),v13(i1-1),v13(i1),x2(i))
    		  v24(i)=linear1d(x1(i1-1),x1(i1),v14(i1-1),v14(i1),x2(i))
    		  !
    		  exit
    		  !
        elseif(x2(i)>=x1(dim1)) then
          !
          v21(i)=v11(dim1)
          v22(i)=v12(dim1)
          v23(i)=v13(dim1)
          v24(i)=v14(dim1)
          !
          exit
          !
        elseif(x2(i)<=x1(0)) then
          !
          v21(i)=v11(0)
          v22(i)=v12(0)
          v23(i)=v13(0)
          v24(i)=v14(0)
          !
          exit
          !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	                    '  ** Interpolating ... ',100*nper/(dim2+1),' %'
    end do
    !
    print*
    !
  end subroutine regularlinearinterp1d4v
  !
  subroutine regularlinearinterp1d5v(x1,v11,v12,v13,v14,v15,dim1,      &
                                     x2,v21,v22,v23,v24,v25,dim2)
    !
    integer,intent(in) :: dim1,dim2
    !
    real(8),intent(in) :: x1(0:dim1),v11(0:dim1),v12(0:dim1),          &
                          v13(0:dim1),v14(0:dim1),v15(0:dim1)
    real(8),intent(in) :: x2(0:dim2)
    real(8),intent(out) :: v21(0:dim2),v22(0:dim2),v23(0:dim2),        &
                           v24(0:dim2),v25(0:dim2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    real(8) :: dist,damp
    !
    !
    nper=0
    do i=0,dim2
      !
    	do i1=1,dim1
        !
    		if( x2(i)>=x1(i1-1) .and. x2(i)<=x1(i1) ) then
          !
          v21(i)=linear1d(x1(i1-1),x1(i1),v11(i1-1),v11(i1),x2(i))
          v22(i)=linear1d(x1(i1-1),x1(i1),v12(i1-1),v12(i1),x2(i))
          v23(i)=linear1d(x1(i1-1),x1(i1),v13(i1-1),v13(i1),x2(i))
          v24(i)=linear1d(x1(i1-1),x1(i1),v14(i1-1),v14(i1),x2(i))
          v25(i)=linear1d(x1(i1-1),x1(i1),v15(i1-1),v15(i1),x2(i))
          !
          exit
          !
        elseif(x2(i)>x1(dim1)) then
          !
          ! dist=x2(i)-x1(dim1)
          ! damp=dexp(-dist**2)
          ! !
          ! v21(i)=v11(dim1)*damp
          ! v22(i)=v12(dim1)*damp
          ! v23(i)=v13(dim1)*damp
          ! v24(i)=v14(dim1)*damp
          ! v25(i)=v15(dim1)*damp
          v21(i)=v11(dim1)
          v22(i)=v12(dim1)
          v23(i)=v13(dim1)
          v24(i)=v14(dim1)
          v25(i)=v15(dim1)
          !
        end if
        !
    	end do
      !
      nper=nper+1
      !
      ! write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
      !                     '  ** Interpolating ... ',100*nper/(dim2+1),' %'
    end do
    !
    ! print*
    !
  end subroutine regularlinearinterp1d5v
  !
  subroutine regularlinearinterp1d6v(x1,v11,v12,v13,v14,v15,v16,dim1,  &
                                     x2,v21,v22,v23,v24,v25,v26,dim2)
    !
    integer,intent(in) :: dim1,dim2
    !
    real(8),intent(in) :: x1(0:dim1),v11(0:dim1),v12(0:dim1),          &
                          v13(0:dim1),v14(0:dim1),v15(0:dim1),         &
                          v16(0:dim1)
    real(8),intent(in) :: x2(0:dim2)
    real(8),intent(out) :: v21(0:dim2),v22(0:dim2),v23(0:dim2),        &
                           v24(0:dim2),v25(0:dim2),v26(0:dim2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    !
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i)>=x1(i1-1) .and. x2(i)<=x1(i1) ) then
    		  !
    		  v21(i)=linear1d(x1(i1-1),x1(i1),v11(i1-1),v11(i1),x2(i))
    		  v22(i)=linear1d(x1(i1-1),x1(i1),v12(i1-1),v12(i1),x2(i))
    		  v23(i)=linear1d(x1(i1-1),x1(i1),v13(i1-1),v13(i1),x2(i))
    		  v24(i)=linear1d(x1(i1-1),x1(i1),v14(i1-1),v14(i1),x2(i))
    		  v25(i)=linear1d(x1(i1-1),x1(i1),v15(i1-1),v15(i1),x2(i))
    		  v26(i)=linear1d(x1(i1-1),x1(i1),v16(i1-1),v16(i1),x2(i))
    		  !
    		  exit
    		  !
            elseif(x2(i)>=x1(dim1) ) then
              !
              v21(i)=linear1d(x1(dim1-1),x1(dim1),v11(dim1-1),v11(dim1),x2(i))
              v22(i)=linear1d(x1(dim1-1),x1(dim1),v12(dim1-1),v12(dim1),x2(i))
              v23(i)=linear1d(x1(dim1-1),x1(dim1),v13(dim1-1),v13(dim1),x2(i))
              v24(i)=linear1d(x1(dim1-1),x1(dim1),v14(dim1-1),v14(dim1),x2(i))
              v25(i)=linear1d(x1(dim1-1),x1(dim1),v15(dim1-1),v15(dim1),x2(i))
              v26(i)=linear1d(x1(dim1-1),x1(dim1),v16(dim1-1),v16(dim1),x2(i))
              !
              exit
            elseif(x2(i)<=x1(0) ) then
              !
              v21(i)=linear1d(x1(0),x1(1),v11(0),v11(1),x2(i))
              v22(i)=linear1d(x1(0),x1(1),v12(0),v12(1),x2(i))
              v23(i)=linear1d(x1(0),x1(1),v13(0),v13(1),x2(i))
              v24(i)=linear1d(x1(0),x1(1),v14(0),v14(1),x2(i))
              v25(i)=linear1d(x1(0),x1(1),v15(0),v15(1),x2(i))
              v26(i)=linear1d(x1(0),x1(1),v16(0),v16(1),x2(i))
              !
              exit
            end if
    		!
    	end do
    	!
    	nper=nper+1
    	! !
    	! write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	!                     '  ** Interpolating ... ',100*nper/(dim2+1),' %'
    end do
    !
    ! print*
    !
  end subroutine regularlinearinterp1d6v
  !
  subroutine regularlinearinterp1d7v(x1,v11,v12,v13,v14,v15,v16,v17,   &
                                                                dim1,  &
                                     x2,v21,v22,v23,v24,v25,v26,v27,   &
                                                                dim2)
    !
    integer,intent(in) :: dim1,dim2
    !
    real(8),intent(in) :: x1(0:dim1),v11(0:dim1),v12(0:dim1),          &
                          v13(0:dim1),v14(0:dim1),v15(0:dim1),         &
                          v16(0:dim1),v17(0:dim1)
    real(8),intent(in) :: x2(0:dim2)
    real(8),intent(out) :: v21(0:dim2),v22(0:dim2),v23(0:dim2),        &
                           v24(0:dim2),v25(0:dim2),v26(0:dim2),        &
                           v27(0:dim2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    !
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i)>=x1(i1-1) .and. x2(i)<=x1(i1) ) then
    		  !
    		  v21(i)=linear1d(x1(i1-1),x1(i1),v11(i1-1),v11(i1),x2(i))
    		  v22(i)=linear1d(x1(i1-1),x1(i1),v12(i1-1),v12(i1),x2(i))
    		  v23(i)=linear1d(x1(i1-1),x1(i1),v13(i1-1),v13(i1),x2(i))
    		  v24(i)=linear1d(x1(i1-1),x1(i1),v14(i1-1),v14(i1),x2(i))
    		  v25(i)=linear1d(x1(i1-1),x1(i1),v15(i1-1),v15(i1),x2(i))
    		  v26(i)=linear1d(x1(i1-1),x1(i1),v16(i1-1),v16(i1),x2(i))
    		  v27(i)=linear1d(x1(i1-1),x1(i1),v17(i1-1),v17(i1),x2(i))
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	                    '  ** Interpolating ... ',100*nper/(dim2+1),' %'
    end do
    !
    print*
    !
  end subroutine regularlinearinterp1d7v
  !
  subroutine regularlinearinterp1d8v(x1,v11,v12,v13,v14,v15,v16,v17,   &
                                                            v18,dim1,  &
                                     x2,v21,v22,v23,v24,v25,v26,v27,   &
                                                            v28,dim2)
    !
    integer,intent(in) :: dim1,dim2
    !
    real(8),intent(in) :: x1(0:dim1),v11(0:dim1),v12(0:dim1),          &
                          v13(0:dim1),v14(0:dim1),v15(0:dim1),         &
                          v16(0:dim1),v17(0:dim1),v18(0:dim1)
    real(8),intent(in) :: x2(0:dim2)
    real(8),intent(out) :: v21(0:dim2),v22(0:dim2),v23(0:dim2),        &
                           v24(0:dim2),v25(0:dim2),v26(0:dim2),        &
                           v27(0:dim2),v28(0:dim2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    !
    !
    nper=0
    do i=0,dim2
      !
      do i1=1,dim1
        !
        if( x2(i)>=x1(i1-1) .and. x2(i)<=x1(i1) ) then
          !
          v21(i)=linear1d(x1(i1-1),x1(i1),v11(i1-1),v11(i1),x2(i))
          v22(i)=linear1d(x1(i1-1),x1(i1),v12(i1-1),v12(i1),x2(i))
          v23(i)=linear1d(x1(i1-1),x1(i1),v13(i1-1),v13(i1),x2(i))
          v24(i)=linear1d(x1(i1-1),x1(i1),v14(i1-1),v14(i1),x2(i))
          v25(i)=linear1d(x1(i1-1),x1(i1),v15(i1-1),v15(i1),x2(i))
          v26(i)=linear1d(x1(i1-1),x1(i1),v16(i1-1),v16(i1),x2(i))
          v27(i)=linear1d(x1(i1-1),x1(i1),v17(i1-1),v17(i1),x2(i))
          v28(i)=linear1d(x1(i1-1),x1(i1),v18(i1-1),v18(i1),x2(i))
          !
          exit
          !
        end if
        !
      end do
      !
      nper=nper+1
      !
      write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
                          '  ** Interpolating ... ',100*nper/(dim2+1),' %'
    end do
    !
    print*
    !
  end subroutine regularlinearinterp1d8v
  !
  subroutine regularlinearinterp1d9v(x1,v11,v12,v13,v14,v15,v16,v17,   &
                                                        v18,v19,dim1,  &
                                     x2,v21,v22,v23,v24,v25,v26,v27,   &
                                                        v28,v29,dim2)
    !
    integer,intent(in) :: dim1,dim2
    !
    real(8),intent(in) :: x1(0:dim1),v11(0:dim1),v12(0:dim1),          &
                          v13(0:dim1),v14(0:dim1),v15(0:dim1),         &
                          v16(0:dim1),v17(0:dim1),v18(0:dim1),         &
                          v19(0:dim1)
    real(8),intent(in) :: x2(0:dim2)
    real(8),intent(out) :: v21(0:dim2),v22(0:dim2),v23(0:dim2),        &
                           v24(0:dim2),v25(0:dim2),v26(0:dim2),        &
                           v27(0:dim2),v28(0:dim2),v29(0:dim2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    !
    !
    nper=0
    do i=0,dim2
      !
      do i1=1,dim1
        !
        if( x2(i)>=x1(i1-1) .and. x2(i)<=x1(i1) ) then
          !
          v21(i)=linear1d(x1(i1-1),x1(i1),v11(i1-1),v11(i1),x2(i))
          v22(i)=linear1d(x1(i1-1),x1(i1),v12(i1-1),v12(i1),x2(i))
          v23(i)=linear1d(x1(i1-1),x1(i1),v13(i1-1),v13(i1),x2(i))
          v24(i)=linear1d(x1(i1-1),x1(i1),v14(i1-1),v14(i1),x2(i))
          v25(i)=linear1d(x1(i1-1),x1(i1),v15(i1-1),v15(i1),x2(i))
          v26(i)=linear1d(x1(i1-1),x1(i1),v16(i1-1),v16(i1),x2(i))
          v27(i)=linear1d(x1(i1-1),x1(i1),v17(i1-1),v17(i1),x2(i))
          v28(i)=linear1d(x1(i1-1),x1(i1),v18(i1-1),v18(i1),x2(i))
          v29(i)=linear1d(x1(i1-1),x1(i1),v19(i1-1),v19(i1),x2(i))
          !
          exit
          !
        end if
        !
      end do
      !
      nper=nper+1
      !
      write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
                          '  ** Interpolating ... ',100*nper/(dim2+1),' %'
    end do
    !
    print*
    !
  end subroutine regularlinearinterp1d9v
  !
  !+-------------------------------------------------------------------+
  !| The end of the subroutine regularlinearinterp1d                   |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to interpolate 3D flowfield between 2     |
  !| regular meshes with linear approach.                              |
  !+-------------------------------------------------------------------+
  subroutine regularlinearinterp3d1v(x1,y1,z1,v11,dim1,djm1,dkm1,      &
                                     x2,y2,z2,v21,dim2,djm2,dkm2)
    !
    integer,intent(in) :: dim1,djm1,dkm1,dim2,djm2,dkm2
    !
    real(8),intent(in) :: x1(0:dim1,0:djm1,0:dkm1),                    &
                          y1(0:dim1,0:djm1,0:dkm1),                    &
                          z1(0:dim1,0:djm1,0:dkm1),                    &
                          v11(0:dim1,0:djm1,0:dkm1)
    real(8),intent(in) :: x2(0:dim2,0:djm2,0:dkm2),                    &
                          y2(0:dim2,0:djm2,0:dkm2),                    &
                          z2(0:dim2,0:djm2,0:dkm2)
    real(8),intent(out) :: v21(0:dim2,0:djm2,0:dkm2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    !
    real(8) :: vti(0:dim2,0:djm1,0:dkm1),vtj(0:dim2,0:djm2,0:dkm1)
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i,0,0)>=x1(i1-1,0,0) .and. x2(i,0,0)<=x1(i1,0,0) ) then
    		  !
    		  do k1=0,dkm1
    		  do j1=0,djm1
    		    vti(i,j1,k1)=linear1d(x1(i1-1,0,0), x1(i1,0,0),            &
    		                          v11(i1-1,j1,k1),v11(i1,j1,k1),       &
    		                          x2(i,0,0)                            )
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    !
    do j=0,djm2
    	!
    	do j1=1,djm1
    		!
    		if( y2(0,j,0)>=y1(0,j1-1,0) .and. y2(0,j,0)<=y1(0,j1,0) ) then
    		  do k1=0,dkm1
    		  do i=0,dim2
    		    vtj(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),            &
    		                         vti(i,j1-1,k1),vti(i,j1,k1),          &
    		                         y2(0,j,0))
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    !
    do k=0,dkm2
    	!
    	do k1=1,dkm1
    		!
    		if( z2(0,0,k)>=z1(0,0,k1-1) .and. z2(0,0,k)<=z1(0,0,k1) ) then
    		  do j=0,djm2
    		  do i=0,dim2
    		    v21(i,j,k)=linear1d(z1(0,0,k1-1),  z1(0,0,k1),            &
    		                         vtj(i,j,k1-1),vtj(i,j,k1),          &
    		                         z2(0,0,k))
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    print*
    !
  end subroutine regularlinearinterp3d1v
  !!
  subroutine regularlinearinterp3d2v(x1,y1,z1,                         &
                                     v11,v12,                          &
                                     dim1,djm1,dkm1,                   &
                                     x2,y2,z2,                         &
                                     v21,v22,                          &
                                     dim2,djm2,dkm2                    )
    !
    integer,intent(in) :: dim1,djm1,dkm1,dim2,djm2,dkm2
    !
    real(8),intent(in) :: x1(0:dim1,0:djm1,0:dkm1),                    &
                          y1(0:dim1,0:djm1,0:dkm1),                    &
                          z1(0:dim1,0:djm1,0:dkm1),                    &
                          v11(0:dim1,0:djm1,0:dkm1),                    &
                          v12(0:dim1,0:djm1,0:dkm1)
    real(8),intent(in) :: x2(0:dim2,0:djm2,0:dkm2),                    &
                          y2(0:dim2,0:djm2,0:dkm2),                    &
                          z2(0:dim2,0:djm2,0:dkm2)
    real(8),intent(out) :: v21(0:dim2,0:djm2,0:dkm2),                   &
                           v22(0:dim2,0:djm2,0:dkm2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    !
    real(8) :: vti1(0:dim2,0:djm1,0:dkm1),vtj1(0:dim2,0:djm2,0:dkm1),  &
               vti2(0:dim2,0:djm1,0:dkm1),vtj2(0:dim2,0:djm2,0:dkm1)
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i,0,0)>=x1(i1-1,0,0) .and. x2(i,0,0)<=x1(i1,0,0) ) then
    		  !
    		  do k1=0,dkm1
    		  do j1=0,djm1
    		    vti1(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v11(i1-1,j1,k1),v11(i1,j1,k1),       &
    		                          x2(i,0,0)                            )
    		    vti2(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v12(i1-1,j1,k1),v12(i1,j1,k1),       &
    		                          x2(i,0,0)                            )
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    !
    do j=0,djm2
    	!
    	do j1=1,djm1
    		!
    		if( y2(0,j,0)>=y1(0,j1-1,0) .and. y2(0,j,0)<=y1(0,j1,0) ) then
    		  do k1=0,dkm1
    		  do i=0,dim2
    		    vtj1(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti1(i,j1-1,k1),vti1(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj2(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti2(i,j1-1,k1),vti2(i,j1,k1),        &
    		                         y2(0,j,0))
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    !
    do k=0,dkm2
    	!
    	do k1=1,dkm1
    		!
    		if( z2(0,0,k)>=z1(0,0,k1-1) .and. z2(0,0,k)<=z1(0,0,k1) ) then
    		  do j=0,djm2
    		  do i=0,dim2
    		    v21(i,j,k)=linear1d(z1(0,0,k1-1),     z1(0,0,k1),          &
    		                         vtj1(i,j,k1-1),vtj1(i,j,k1),          &
    		                         z2(0,0,k))
    		    v22(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj2(i,j,k1-1),vtj2(i,j,k1),          &
    		                         z2(0,0,k))
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    print*
    !
  end subroutine regularlinearinterp3d2v
  !!
  subroutine regularlinearinterp3d3v(x1,y1,z1,                         &
                                     v11,v12,v13,                      &
                                     dim1,djm1,dkm1,                   &
                                     x2,y2,z2,                         &
                                     v21,v22,v23,                      &
                                     dim2,djm2,dkm2                    )
    !
    integer,intent(in) :: dim1,djm1,dkm1,dim2,djm2,dkm2
    !
    real(8),intent(in) :: x1(0:dim1,0:djm1,0:dkm1),                    &
                          y1(0:dim1,0:djm1,0:dkm1),                    &
                          z1(0:dim1,0:djm1,0:dkm1),                    &
                          v11(0:dim1,0:djm1,0:dkm1),                   &
                          v12(0:dim1,0:djm1,0:dkm1),                   &
                          v13(0:dim1,0:djm1,0:dkm1)
    real(8),intent(in) :: x2(0:dim2,0:djm2,0:dkm2),                    &
                          y2(0:dim2,0:djm2,0:dkm2),                    &
                          z2(0:dim2,0:djm2,0:dkm2)
    real(8),intent(out) :: v21(0:dim2,0:djm2,0:dkm2),                  &
                           v22(0:dim2,0:djm2,0:dkm2),                  &
                           v23(0:dim2,0:djm2,0:dkm2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    !
    real(8) :: vti1(0:dim2,0:djm1,0:dkm1),vtj1(0:dim2,0:djm2,0:dkm1),  &
               vti2(0:dim2,0:djm1,0:dkm1),vtj2(0:dim2,0:djm2,0:dkm1),  &
               vti3(0:dim2,0:djm1,0:dkm1),vtj3(0:dim2,0:djm2,0:dkm1)
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i,0,0)>=x1(i1-1,0,0) .and. x2(i,0,0)<=x1(i1,0,0) ) then
    		  !
    		  do k1=0,dkm1
    		  do j1=0,djm1
    		    vti1(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v11(i1-1,j1,k1),v11(i1,j1,k1),       &
    		                          x2(i,0,0)                            )
    		    vti2(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v12(i1-1,j1,k1),v12(i1,j1,k1),       &
    		                          x2(i,0,0)                            )
    		    vti3(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v13(i1-1,j1,k1),v13(i1,j1,k1),       &
    		                          x2(i,0,0)                           )
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    !
    do j=0,djm2
    	!
    	do j1=1,djm1
    		!
    		if( y2(0,j,0)>=y1(0,j1-1,0) .and. y2(0,j,0)<=y1(0,j1,0) ) then
    		  do k1=0,dkm1
    		  do i=0,dim2
    		    vtj1(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti1(i,j1-1,k1),vti1(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj2(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti2(i,j1-1,k1),vti2(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj3(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti3(i,j1-1,k1),vti3(i,j1,k1),        &
    		                         y2(0,j,0))
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    !
    do k=0,dkm2
    	!
    	do k1=1,dkm1
    		!
    		if( z2(0,0,k)>=z1(0,0,k1-1) .and. z2(0,0,k)<=z1(0,0,k1) ) then
    		  do j=0,djm2
    		  do i=0,dim2
    		    v21(i,j,k)=linear1d(z1(0,0,k1-1),     z1(0,0,k1),          &
    		                         vtj1(i,j,k1-1),vtj1(i,j,k1),          &
    		                         z2(0,0,k))
    		    v22(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj2(i,j,k1-1),vtj2(i,j,k1),          &
    		                         z2(0,0,k))
    		    v23(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj3(i,j,k1-1),vtj3(i,j,k1),          &
    		                         z2(0,0,k))
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    print*
    !
  end subroutine regularlinearinterp3d3v
  !!
  subroutine regularlinearinterp3d4v(x1,y1,z1,                         &
                                     v11,v12,v13,v14,                  &
                                     dim1,djm1,dkm1,                   &
                                     x2,y2,z2,                         &
                                     v21,v22,v23,v24,                  &
                                     dim2,djm2,dkm2                    )
    !
    integer,intent(in) :: dim1,djm1,dkm1,dim2,djm2,dkm2
    !
    real(8),intent(in) :: x1(0:dim1,0:djm1,0:dkm1),                    &
                          y1(0:dim1,0:djm1,0:dkm1),                    &
                          z1(0:dim1,0:djm1,0:dkm1),                    &
                          v11(0:dim1,0:djm1,0:dkm1),                   &
                          v12(0:dim1,0:djm1,0:dkm1),                   &
                          v13(0:dim1,0:djm1,0:dkm1),                   &
                          v14(0:dim1,0:djm1,0:dkm1)
    real(8),intent(in) :: x2(0:dim2,0:djm2,0:dkm2),                    &
                          y2(0:dim2,0:djm2,0:dkm2),                    &
                          z2(0:dim2,0:djm2,0:dkm2)
    real(8),intent(out) :: v21(0:dim2,0:djm2,0:dkm2),                  &
                           v22(0:dim2,0:djm2,0:dkm2),                  &
                           v23(0:dim2,0:djm2,0:dkm2),                  &
                           v24(0:dim2,0:djm2,0:dkm2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    !
    real(8) :: vti1(0:dim2,0:djm1,0:dkm1),vtj1(0:dim2,0:djm2,0:dkm1),  &
               vti2(0:dim2,0:djm1,0:dkm1),vtj2(0:dim2,0:djm2,0:dkm1),  &
               vti3(0:dim2,0:djm1,0:dkm1),vtj3(0:dim2,0:djm2,0:dkm1),  &
               vti4(0:dim2,0:djm1,0:dkm1),vtj4(0:dim2,0:djm2,0:dkm1)
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i,0,0)>=x1(i1-1,0,0) .and. x2(i,0,0)<=x1(i1,0,0) ) then
    		  !
    		  do k1=0,dkm1
    		  do j1=0,djm1
    		    vti1(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v11(i1-1,j1,k1),v11(i1,j1,k1),       &
    		                          x2(i,0,0)                            )
    		    vti2(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v12(i1-1,j1,k1),v12(i1,j1,k1),       &
    		                          x2(i,0,0)                            )
    		    vti3(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v13(i1-1,j1,k1),v13(i1,j1,k1),       &
    		                          x2(i,0,0)                            )
    		    vti4(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v14(i1-1,j1,k1),v14(i1,j1,k1),       &
    		                          x2(i,0,0)                            )
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    !
    do j=0,djm2
    	!
    	do j1=1,djm1
    		!
    		if( y2(0,j,0)>=y1(0,j1-1,0) .and. y2(0,j,0)<=y1(0,j1,0) ) then
    		  do k1=0,dkm1
    		  do i=0,dim2
    		    vtj1(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti1(i,j1-1,k1),vti1(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj2(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti2(i,j1-1,k1),vti2(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj3(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti3(i,j1-1,k1),vti3(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj4(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti4(i,j1-1,k1),vti4(i,j1,k1),        &
    		                         y2(0,j,0))
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    !
    do k=0,dkm2
    	!
    	do k1=1,dkm1
    		!
    		if( z2(0,0,k)>=z1(0,0,k1-1) .and. z2(0,0,k)<=z1(0,0,k1) ) then
    		  do j=0,djm2
    		  do i=0,dim2
    		    v21(i,j,k)=linear1d(z1(0,0,k1-1),     z1(0,0,k1),          &
    		                         vtj1(i,j,k1-1),vtj1(i,j,k1),          &
    		                         z2(0,0,k))
    		    v22(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj2(i,j,k1-1),vtj2(i,j,k1),          &
    		                         z2(0,0,k))
    		    v23(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj3(i,j,k1-1),vtj3(i,j,k1),          &
    		                         z2(0,0,k))
    		    v24(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj4(i,j,k1-1),vtj4(i,j,k1),          &
    		                         z2(0,0,k))
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    print*
    !
  end subroutine regularlinearinterp3d4v
  !!
  subroutine regularlinearinterp3d5v(x1,y1,z1,                         &
                                     v11,v12,v13,v14,v15,              &
                                     dim1,djm1,dkm1,                   &
                                     x2,y2,z2,                         &
                                     v21,v22,v23,v24,v25,              &
                                     dim2,djm2,dkm2                    )
    !
    integer,intent(in) :: dim1,djm1,dkm1,dim2,djm2,dkm2
    !
    real(8),intent(in) :: x1(0:dim1,0:djm1,0:dkm1),                    &
                          y1(0:dim1,0:djm1,0:dkm1),                    &
                          z1(0:dim1,0:djm1,0:dkm1),                    &
                          v11(0:dim1,0:djm1,0:dkm1),                   &
                          v12(0:dim1,0:djm1,0:dkm1),                   &
                          v13(0:dim1,0:djm1,0:dkm1),                   &
                          v14(0:dim1,0:djm1,0:dkm1),                   &
                          v15(0:dim1,0:djm1,0:dkm1)
    real(8),intent(in) :: x2(0:dim2,0:djm2,0:dkm2),                    &
                          y2(0:dim2,0:djm2,0:dkm2),                    &
                          z2(0:dim2,0:djm2,0:dkm2)
    real(8),intent(out) :: v21(0:dim2,0:djm2,0:dkm2),                  &
                           v22(0:dim2,0:djm2,0:dkm2),                  &
                           v23(0:dim2,0:djm2,0:dkm2),                  &
                           v24(0:dim2,0:djm2,0:dkm2),                  &
                           v25(0:dim2,0:djm2,0:dkm2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    !
    real(8) :: vti1(0:dim2,0:djm1,0:dkm1),vtj1(0:dim2,0:djm2,0:dkm1),  &
               vti2(0:dim2,0:djm1,0:dkm1),vtj2(0:dim2,0:djm2,0:dkm1),  &
               vti3(0:dim2,0:djm1,0:dkm1),vtj3(0:dim2,0:djm2,0:dkm1),  &
               vti4(0:dim2,0:djm1,0:dkm1),vtj4(0:dim2,0:djm2,0:dkm1),  &
               vti5(0:dim2,0:djm1,0:dkm1),vtj5(0:dim2,0:djm2,0:dkm1)
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i,0,0)>=x1(i1-1,0,0) .and. x2(i,0,0)<=x1(i1,0,0) ) then
    		  !
    		  do k1=0,dkm1
    		  do j1=0,djm1
    		    vti1(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v11(i1-1,j1,k1),v11(i1,j1,k1),       &
    		                          x2(i,0,0)                            )
    		    vti2(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v12(i1-1,j1,k1),v12(i1,j1,k1),       &
    		                          x2(i,0,0)                            )
    		    vti3(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v13(i1-1,j1,k1),v13(i1,j1,k1),       &
    		                          x2(i,0,0)                            )
    		    vti4(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v14(i1-1,j1,k1),v14(i1,j1,k1),       &
    		                          x2(i,0,0)                            )
    		    vti5(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v15(i1-1,j1,k1),v15(i1,j1,k1),       &
    		                          x2(i,0,0)                            )
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    !
    do j=0,djm2
    	!
    	do j1=1,djm1
    		!
    		if( y2(0,j,0)>=y1(0,j1-1,0) .and. y2(0,j,0)<=y1(0,j1,0) ) then
    		  do k1=0,dkm1
    		  do i=0,dim2
    		    vtj1(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti1(i,j1-1,k1),vti1(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj2(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti2(i,j1-1,k1),vti2(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj3(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti3(i,j1-1,k1),vti3(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj4(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti4(i,j1-1,k1),vti4(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj5(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti5(i,j1-1,k1),vti5(i,j1,k1),        &
    		                         y2(0,j,0))
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    !
    do k=0,dkm2
    	!
    	do k1=1,dkm1
    		!
    		if( z2(0,0,k)>=z1(0,0,k1-1) .and. z2(0,0,k)<=z1(0,0,k1) ) then
    		  do j=0,djm2
    		  do i=0,dim2
    		    v21(i,j,k)=linear1d(z1(0,0,k1-1),     z1(0,0,k1),          &
    		                         vtj1(i,j,k1-1),vtj1(i,j,k1),          &
    		                         z2(0,0,k))
    		    v22(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj2(i,j,k1-1),vtj2(i,j,k1),          &
    		                         z2(0,0,k))
    		    v23(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj3(i,j,k1-1),vtj3(i,j,k1),          &
    		                         z2(0,0,k))
    		    v24(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj4(i,j,k1-1),vtj4(i,j,k1),          &
    		                         z2(0,0,k))
    		    v25(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj5(i,j,k1-1),vtj5(i,j,k1),          &
    		                         z2(0,0,k))
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    print*
    !
  end subroutine regularlinearinterp3d5v
  !
  subroutine regularlinearinterp3d6v_cublinear(x1,y1,z1,                         &
                                     v11,v12,v13,v14,v15,v16,          &
                                     x2,y2,z2,                         &
                                     v21,v22,v23,v24,v25,v26 )
    !
    ! arguments
    real(8),intent(in), dimension(0:,0:,0:) :: x1,y1,z1,v11,v12,v13,v14,v15,v16
    real(8),intent(in), dimension(0:,0:,0:) :: x2,y2,z2
    real(8),intent(out),dimension(0:,0:,0:) :: v21,v22,v23,v24,v25,v26
    !
    ! local data
    integer :: dim1,djm1,dkm1,dim2,djm2,dkm2,i1,j1,k1,i2,j2,k2,nper,ic,jc,kc
    real(8) :: xd,yd,zd,c00,c01,c10,c11,c0,c1,rper
    real(8) :: p(3)
    logical :: lfound
    !
    !$ SAVE p
    !
    !$OMP THREADPRIVATE(p)
    !
    dim1=size(x1,1)-1
    djm1=size(x1,2)-1
    dkm1=size(x1,3)-1
    !
    dim2=size(x2,1)-1
    djm2=size(x2,2)-1
    dkm2=size(x2,3)-1
    !
    p(1)=x2(1,1,1)
    p(2)=y2(1,1,1)
    p(3)=z2(1,1,1)
    !
    call ijksearch(p,ic,jc,kc,lfound,initial=.true.,xin=x1,yin=y1,zin=z1)
    !
    print*,' ** interpolating'
    !
    nper=0
    !
    !$OMP parallel default(shared) private(i1,j1,k1,i2,j2,k2,xd,yd,zd,c00,c01,c10,c11,c0,c1,nper,rper)
    !$OMP do
    !
    do k2=0,dkm2
    do j2=0,djm2
    do i2=0,dim2
      !
      p(1)=x2(i2,j2,k2)
      p(2)=y2(i2,j2,k2)
      p(3)=z2(i2,j2,k2)
      !
      call ijksearch(p,i1,j1,k1,lfound,initial=.false.)
      !
      if(.not. lfound) then
        print*,' cant find a cell containg p at: ',i2,j2,k2
        stop
      endif
      !
      xd=(x2(i2,j2,k2)-x1(i1-1,j1,k1))/(x1(i1,j1,k1)-x1(i1-1,j1,k1))
      yd=(y2(i2,j2,k2)-y1(i1,j1-1,k1))/(y1(i1,j1,k1)-y1(i1,j1-1,k1))
      zd=(z2(i2,j2,k2)-z1(i1,j1,k1-1))/(z1(i1,j1,k1)-z1(i1,j1,k1-1))
      !
      c00=v11(i1-1,j1-1,k1-1)*(1.d0-xd)+v11(i1,j1-1,k1-1)*xd
      c01=v11(i1-1,j1,  k1-1)*(1.d0-xd)+v11(i1,j1,  k1-1)*xd
      c10=v11(i1-1,j1,  k1)  *(1.d0-xd)+v11(i1,j1,  k1)  *xd
      c11=v11(i1-1,j1-1,k1)  *(1.d0-xd)+v11(i1,j1-1,k1)  *xd
      !
      c0=c00*(1.d0-yd)+c10*yd
      c1=c01*(1.d0-yd)+c11*yd
      !
      v21(i2,j2,k2)=c0*(1.d0-zd)+c1*zd
      !
      c00=v12(i1-1,j1-1,k1-1)*(1.d0-xd)+v12(i1,j1-1,k1-1)*xd
      c01=v12(i1-1,j1,  k1-1)*(1.d0-xd)+v12(i1,j1,  k1-1)*xd
      c10=v12(i1-1,j1,  k1)  *(1.d0-xd)+v12(i1,j1,  k1)  *xd
      c11=v12(i1-1,j1-1,k1)  *(1.d0-xd)+v12(i1,j1-1,k1)  *xd
      !
      c0=c00*(1.d0-yd)+c10*yd
      c1=c01*(1.d0-yd)+c11*yd
      !
      v22(i2,j2,k2)=c0*(1.d0-zd)+c1*zd
      !
      c00=v13(i1-1,j1-1,k1-1)*(1.d0-xd)+v13(i1,j1-1,k1-1)*xd
      c01=v13(i1-1,j1,  k1-1)*(1.d0-xd)+v13(i1,j1,  k1-1)*xd
      c10=v13(i1-1,j1,  k1)  *(1.d0-xd)+v13(i1,j1,  k1)  *xd
      c11=v13(i1-1,j1-1,k1)  *(1.d0-xd)+v13(i1,j1-1,k1)  *xd
      !
      c0=c00*(1.d0-yd)+c10*yd
      c1=c01*(1.d0-yd)+c11*yd
      !
      v23(i2,j2,k2)=c0*(1.d0-zd)+c1*zd
      !
      c00=v14(i1-1,j1-1,k1-1)*(1.d0-xd)+v14(i1,j1-1,k1-1)*xd
      c01=v14(i1-1,j1,  k1-1)*(1.d0-xd)+v14(i1,j1,  k1-1)*xd
      c10=v14(i1-1,j1,  k1)  *(1.d0-xd)+v14(i1,j1,  k1)  *xd
      c11=v14(i1-1,j1-1,k1)  *(1.d0-xd)+v14(i1,j1-1,k1)  *xd
      !
      c0=c00*(1.d0-yd)+c10*yd
      c1=c01*(1.d0-yd)+c11*yd
      !
      v24(i2,j2,k2)=c0*(1.d0-zd)+c1*zd
      !
      c00=v15(i1-1,j1-1,k1-1)*(1.d0-xd)+v15(i1,j1-1,k1-1)*xd
      c01=v15(i1-1,j1,  k1-1)*(1.d0-xd)+v15(i1,j1,  k1-1)*xd
      c10=v15(i1-1,j1,  k1)  *(1.d0-xd)+v15(i1,j1,  k1)  *xd
      c11=v15(i1-1,j1-1,k1)  *(1.d0-xd)+v15(i1,j1-1,k1)  *xd
      !
      c0=c00*(1.d0-yd)+c10*yd
      c1=c01*(1.d0-yd)+c11*yd
      !
      v25(i2,j2,k2)=c0*(1.d0-zd)+c1*zd
      !
      c00=v16(i1-1,j1-1,k1-1)*(1.d0-xd)+v16(i1,j1-1,k1-1)*xd
      c01=v16(i1-1,j1,  k1-1)*(1.d0-xd)+v16(i1,j1,  k1-1)*xd
      c10=v16(i1-1,j1,  k1)  *(1.d0-xd)+v16(i1,j1,  k1)  *xd
      c11=v16(i1-1,j1-1,k1)  *(1.d0-xd)+v16(i1,j1-1,k1)  *xd
      !
      c0=c00*(1.d0-yd)+c10*yd
      c1=c01*(1.d0-yd)+c11*yd
      !
      v26(i2,j2,k2)=c0*(1.d0-zd)+c1*zd
      !
    enddo
    enddo
      !
      nper=nper+1
      !
      rper=100.d0*dble(nper*128)/dble(dkm2)
      !
      write(*,'(1A1,A,F8.4,A,$)')char(13),                            &
                                     '  ** Interpolating ... ',rper,' %'
    enddo
    !$OMP end do
    !$OMP end parallel
    !
    write(*,*)
    !
  end subroutine regularlinearinterp3d6v_cublinear
  subroutine regularlinearinterp3d6v_straight(x1,y1,z1,                         &
                                     v11,v12,v13,v14,v15,v16,          &
                                     x2,y2,z2,                         &
                                     v21,v22,v23,v24,v25,v26 )
    !
    ! arguments
    real(8),intent(in), dimension(0:) :: x1,y1,z1
    real(8),intent(in), dimension(0:,0:,0:) :: v11,v12,v13,v14,v15,v16
    real(8),intent(in), dimension(0:,0:,0:) :: x2,y2,z2
    real(8),intent(out),dimension(0:,0:,0:) :: v21,v22,v23,v24,v25,v26
    !
    ! local data
    integer :: dim1,djm1,dkm1,dim2,djm2,dkm2,i1,j1,k1,i2,j2,k2,nper,ic,jc,kc
    real(8) :: xd,yd,zd,c00,c01,c10,c11,c0,c1,rper
    real(8) :: p(3)
    logical :: lfound
    !
    !$ SAVE p
    !
    !$OMP THREADPRIVATE(p)
    !
    dim1=size(x1)-1
    djm1=size(y1)-1
    dkm1=size(z1)-1
    !
    dim2=size(x2,1)-1
    djm2=size(x2,2)-1
    dkm2=size(x2,3)-1
    !
    p(1)=x2(1,1,1)
    p(2)=y2(1,1,1)
    p(3)=z2(1,1,1)
    !
    print*,' ** interpolating'
    !
    nper=0
    !
    !$OMP parallel default(shared) private(i1,j1,k1,i2,j2,k2,xd,yd,zd,c00,c01,c10,c11,c0,c1,nper,rper)
    !$OMP do
    !
    do k2=0,dkm2
    do j2=0,djm2
    do i2=0,dim2
      !
      loop1: do i1=1,dim1
        !
        if(x2(i2,j2,k2)>=x1(i1-1) .and. x2(i2,j2,k2)<=x1(i1)) then
          !
          do j1=1,djm1
            !
            if(y2(i2,j2,k2)>=y1(j1-1) .and. y2(i2,j2,k2)<=y1(j1)) then
              !
              do k1=1,dkm1
                !
                if(z2(i2,j2,k2)>=z1(k1-1) .and. z2(i2,j2,k2)<=z1(k1)) then
                  !
                  xd=(x2(i2,j2,k2)-x1(i1-1))/(x1(i1)-x1(i1-1))
                  yd=(y2(i2,j2,k2)-y1(j1-1))/(y1(j1)-y1(j1-1))
                  zd=(z2(i2,j2,k2)-z1(k1-1))/(z1(k1)-z1(k1-1))
                  !
                  c00=v11(i1-1,j1-1,k1-1)*(1.d0-xd)+v11(i1,j1-1,k1-1)*xd
                  c01=v11(i1-1,j1,  k1-1)*(1.d0-xd)+v11(i1,j1,  k1-1)*xd
                  c10=v11(i1-1,j1,  k1)  *(1.d0-xd)+v11(i1,j1,  k1)  *xd
                  c11=v11(i1-1,j1-1,k1)  *(1.d0-xd)+v11(i1,j1-1,k1)  *xd
                  !
                  c0=c00*(1.d0-yd)+c10*yd
                  c1=c01*(1.d0-yd)+c11*yd
                  !
                  v21(i2,j2,k2)=c0*(1.d0-zd)+c1*zd
                  !
                  c00=v12(i1-1,j1-1,k1-1)*(1.d0-xd)+v12(i1,j1-1,k1-1)*xd
                  c01=v12(i1-1,j1,  k1-1)*(1.d0-xd)+v12(i1,j1,  k1-1)*xd
                  c10=v12(i1-1,j1,  k1)  *(1.d0-xd)+v12(i1,j1,  k1)  *xd
                  c11=v12(i1-1,j1-1,k1)  *(1.d0-xd)+v12(i1,j1-1,k1)  *xd
                  !
                  c0=c00*(1.d0-yd)+c10*yd
                  c1=c01*(1.d0-yd)+c11*yd
                  !
                  v22(i2,j2,k2)=c0*(1.d0-zd)+c1*zd
                  !
                  c00=v13(i1-1,j1-1,k1-1)*(1.d0-xd)+v13(i1,j1-1,k1-1)*xd
                  c01=v13(i1-1,j1,  k1-1)*(1.d0-xd)+v13(i1,j1,  k1-1)*xd
                  c10=v13(i1-1,j1,  k1)  *(1.d0-xd)+v13(i1,j1,  k1)  *xd
                  c11=v13(i1-1,j1-1,k1)  *(1.d0-xd)+v13(i1,j1-1,k1)  *xd
                  !
                  c0=c00*(1.d0-yd)+c10*yd
                  c1=c01*(1.d0-yd)+c11*yd
                  !
                  v23(i2,j2,k2)=c0*(1.d0-zd)+c1*zd
                  !
                  c00=v14(i1-1,j1-1,k1-1)*(1.d0-xd)+v14(i1,j1-1,k1-1)*xd
                  c01=v14(i1-1,j1,  k1-1)*(1.d0-xd)+v14(i1,j1,  k1-1)*xd
                  c10=v14(i1-1,j1,  k1)  *(1.d0-xd)+v14(i1,j1,  k1)  *xd
                  c11=v14(i1-1,j1-1,k1)  *(1.d0-xd)+v14(i1,j1-1,k1)  *xd
                  !
                  c0=c00*(1.d0-yd)+c10*yd
                  c1=c01*(1.d0-yd)+c11*yd
                  !
                  v24(i2,j2,k2)=c0*(1.d0-zd)+c1*zd
                  !
                  c00=v15(i1-1,j1-1,k1-1)*(1.d0-xd)+v15(i1,j1-1,k1-1)*xd
                  c01=v15(i1-1,j1,  k1-1)*(1.d0-xd)+v15(i1,j1,  k1-1)*xd
                  c10=v15(i1-1,j1,  k1)  *(1.d0-xd)+v15(i1,j1,  k1)  *xd
                  c11=v15(i1-1,j1-1,k1)  *(1.d0-xd)+v15(i1,j1-1,k1)  *xd
                  !
                  c0=c00*(1.d0-yd)+c10*yd
                  c1=c01*(1.d0-yd)+c11*yd
                  !
                  v25(i2,j2,k2)=c0*(1.d0-zd)+c1*zd
                  !
                  c00=v16(i1-1,j1-1,k1-1)*(1.d0-xd)+v16(i1,j1-1,k1-1)*xd
                  c01=v16(i1-1,j1,  k1-1)*(1.d0-xd)+v16(i1,j1,  k1-1)*xd
                  c10=v16(i1-1,j1,  k1)  *(1.d0-xd)+v16(i1,j1,  k1)  *xd
                  c11=v16(i1-1,j1-1,k1)  *(1.d0-xd)+v16(i1,j1-1,k1)  *xd
                  !
                  c0=c00*(1.d0-yd)+c10*yd
                  c1=c01*(1.d0-yd)+c11*yd
                  !
                  v26(i2,j2,k2)=c0*(1.d0-zd)+c1*zd
                  !
                  !
                  exit loop1
                  !
                endif
                !
              enddo
              !
            endif
            !
          enddo
          !
        endif
        !
      enddo loop1
      !
    enddo
    enddo
      !
      nper=nper+1
      !
      rper=100.d0*dble(nper*128)/dble(dkm2)
      !
      write(*,'(1A1,A,F8.4,A,$)')char(13),                            &
                                     '  ** Interpolating ... ',rper,' %'
    enddo
    !$OMP end do
    !$OMP end parallel
    !
    write(*,*)
    !
  end subroutine regularlinearinterp3d6v_straight
  !
  subroutine regularlinearinterp3d6v(x1,y1,z1,                         &
                                     v11,v12,v13,v14,v15,v16,          &
                                     dim1,djm1,dkm1,                   &
                                     x2,y2,z2,                         &
                                     v21,v22,v23,v24,v25,v26,          &
                                     dim2,djm2,dkm2                    )
    !
    integer,intent(in) :: dim1,djm1,dkm1,dim2,djm2,dkm2
    !
    real(8),intent(in) :: x1(0:dim1,0:djm1,0:dkm1),                    &
                          y1(0:dim1,0:djm1,0:dkm1),                    &
                          z1(0:dim1,0:djm1,0:dkm1),                    &
                          v11(0:dim1,0:djm1,0:dkm1),                   &
                          v12(0:dim1,0:djm1,0:dkm1),                   &
                          v13(0:dim1,0:djm1,0:dkm1),                   &
                          v14(0:dim1,0:djm1,0:dkm1),                   &
                          v15(0:dim1,0:djm1,0:dkm1),                   &
                          v16(0:dim1,0:djm1,0:dkm1)
    real(8),intent(in) :: x2(0:dim2,0:djm2,0:dkm2),                    &
                          y2(0:dim2,0:djm2,0:dkm2),                    &
                          z2(0:dim2,0:djm2,0:dkm2)
    real(8),intent(out) :: v21(0:dim2,0:djm2,0:dkm2),                  &
                           v22(0:dim2,0:djm2,0:dkm2),                  &
                           v23(0:dim2,0:djm2,0:dkm2),                  &
                           v24(0:dim2,0:djm2,0:dkm2),                  &
                           v25(0:dim2,0:djm2,0:dkm2),                  &
                           v26(0:dim2,0:djm2,0:dkm2)
    !
    integer :: i,j,k,i1,j1,k1,nper
    !
    real(8) :: vti1(0:dim2,0:djm1,0:dkm1),vtj1(0:dim2,0:djm2,0:dkm1),  &
               vti2(0:dim2,0:djm1,0:dkm1),vtj2(0:dim2,0:djm2,0:dkm1),  &
               vti3(0:dim2,0:djm1,0:dkm1),vtj3(0:dim2,0:djm2,0:dkm1),  &
               vti4(0:dim2,0:djm1,0:dkm1),vtj4(0:dim2,0:djm2,0:dkm1),  &
               vti5(0:dim2,0:djm1,0:dkm1),vtj5(0:dim2,0:djm2,0:dkm1),  &
               vti6(0:dim2,0:djm1,0:dkm1),vtj6(0:dim2,0:djm2,0:dkm1)
    !
    real(8) :: xa,za
    !
    nper=0
    do i=0,dim2
    	!
      xa=x2(i,0,0)
      !
      ! if(xa>x1(dim1,0,0))  xa=xa-x1(dim1,0,0) ! for periodic condition
      !
    	do i1=1,dim1
    		!
    		if( xa>=x1(i1-1,0,0) .and. xa<=x1(i1,0,0) ) then
    		  !
    		  do k1=0,dkm1
    		  do j1=0,djm1
    		    vti1(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v11(i1-1,j1,k1),v11(i1,j1,k1),       &
    		                          xa                            )
    		    vti2(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v12(i1-1,j1,k1),v12(i1,j1,k1),       &
    		                          xa                            )
    		    vti3(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v13(i1-1,j1,k1),v13(i1,j1,k1),       &
    		                          xa                            )
    		    vti4(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v14(i1-1,j1,k1),v14(i1,j1,k1),       &
    		                          xa                            )
    		    vti5(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v15(i1-1,j1,k1),v15(i1,j1,k1),       &
    		                          xa                           )
    		    vti6(i,j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),          &
    		                          v16(i1-1,j1,k1),v16(i1,j1,k1),       &
    		                          xa                            )
    		  end do
    		  end do
    		  !
    		  exit
    		  !
            elseif( xa>=x1(dim1,0,0) ) then
              !
              do k1=0,dkm1
              do j1=0,djm1
                vti1(i,j1,k1)=v11(dim1,j1,k1)
                vti2(i,j1,k1)=v12(dim1,j1,k1)
                vti3(i,j1,k1)=v13(dim1,j1,k1)
                vti4(i,j1,k1)=v14(dim1,j1,k1)
                vti5(i,j1,k1)=v15(dim1,j1,k1)
                vti6(i,j1,k1)=v16(dim1,j1,k1)
                !
                ! vti1(i,j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                !                       v11(dim1-1,j1,k1),v11(dim1,j1,k1),       &
                !                       xa                            )
                ! vti2(i,j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                !                       v12(dim1-1,j1,k1),v12(dim1,j1,k1),       &
                !                       xa                            )
                ! vti3(i,j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                !                       v13(dim1-1,j1,k1),v13(dim1,j1,k1),       &
                !                       xa                            )
                ! vti4(i,j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                !                       v14(dim1-1,j1,k1),v14(dim1,j1,k1),       &
                !                       xa                            )
                ! vti5(i,j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                !                       v15(dim1-1,j1,k1),v15(dim1,j1,k1),       &
                !                       xa                            )
                ! vti6(i,j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                !                       v16(dim1-1,j1,k1),v16(dim1,j1,k1),       &
                !                       xa                            )
              end do
              end do
              !
              exit
              !
            elseif( xa<=x1(0,0,0) ) then
              !
              do k1=0,dkm1
              do j1=0,djm1
                vti1(i,j1,k1)=v11(0,j1,k1)
                vti2(i,j1,k1)=v12(0,j1,k1)
                vti3(i,j1,k1)=v13(0,j1,k1)
                vti4(i,j1,k1)=v14(0,j1,k1)
                vti5(i,j1,k1)=v15(0,j1,k1)
                vti6(i,j1,k1)=v16(0,j1,k1)
                ! vti1(i,j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                !                       v11(0,j1,k1),v11(1,j1,k1),       &
                !                       xa                            )
                ! vti2(i,j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                !                       v12(0,j1,k1),v12(1,j1,k1),       &
                !                       xa                            )
                ! vti3(i,j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                !                       v13(0,j1,k1),v13(1,j1,k1),       &
                !                       xa                            )
                ! vti4(i,j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                !                       v14(0,j1,k1),v14(1,j1,k1),       &
                !                       xa                            )
                ! vti5(i,j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                !                       v15(0,j1,k1),v15(1,j1,k1),       &
                !                       xa                            )
                ! vti6(i,j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                !                       v16(0,j1,k1),v16(1,j1,k1),       &
                !                       xa                            )
              end do
              end do
              !
              exit
              !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    !
    do j=0,djm2
    	!
    	do j1=1,djm1
    		!
    		if( y2(0,j,0)>=y1(0,j1-1,0) .and. y2(0,j,0)<=y1(0,j1,0) ) then
    		  do k1=0,dkm1
    		  do i=0,dim2
    		    vtj1(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti1(i,j1-1,k1),vti1(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj2(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti2(i,j1-1,k1),vti2(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj3(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti3(i,j1-1,k1),vti3(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj4(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti4(i,j1-1,k1),vti4(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj5(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti5(i,j1-1,k1),vti5(i,j1,k1),        &
    		                         y2(0,j,0))
    		    vtj6(i,j,k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
    		                         vti6(i,j1-1,k1),vti6(i,j1,k1),        &
    		                         y2(0,j,0))
    		  end do
    		  end do
    		  !
    		  exit
    		  !
    		elseif( y2(0,j,0)>=y1(0,djm1,0) ) then
    		  do k1=0,dkm1
              do i=0,dim2
                vtj1(i,j,k1)=vti1(i,djm1,k1)
                vtj2(i,j,k1)=vti2(i,djm1,k1)
                vtj3(i,j,k1)=vti3(i,djm1,k1)
                vtj4(i,j,k1)=vti4(i,djm1,k1)
                vtj5(i,j,k1)=vti5(i,djm1,k1)
                vtj6(i,j,k1)=vti6(i,djm1,k1)
                !
                ! vtj1(i,j,k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                !                      vti1(i,djm1-1,k1),vti1(i,djm1,k1),        &
                !                      y2(0,j,0))
                ! vtj2(i,j,k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                !                      vti2(i,djm1-1,k1),vti2(i,djm1,k1),        &
                !                      y2(0,j,0))
                ! vtj3(i,j,k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                !                      vti3(i,djm1-1,k1),vti3(i,djm1,k1),        &
                !                      y2(0,j,0))
                ! vtj4(i,j,k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                !                      vti4(i,djm1-1,k1),vti4(i,djm1,k1),        &
                !                      y2(0,j,0))
                ! vtj5(i,j,k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                !                      vti5(i,djm1-1,k1),vti5(i,djm1,k1),        &
                !                      y2(0,j,0))
                ! vtj6(i,j,k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                !                      vti6(i,djm1-1,k1),vti6(i,djm1,k1),        &
                !                      y2(0,j,0))
              end do
              end do
              !
              exit
              !
            elseif( y2(0,j,0)<=y1(0,0,0) ) then
              do k1=0,dkm1
              do i=0,dim2
                vtj1(i,j,k1)=vti1(i,djm1,k1)
                vtj2(i,j,k1)=vti2(i,djm1,k1)
                vtj3(i,j,k1)=vti3(i,djm1,k1)
                vtj4(i,j,k1)=vti4(i,djm1,k1)
                vtj5(i,j,k1)=vti5(i,djm1,k1)
                vtj6(i,j,k1)=vti6(i,djm1,k1)
                !
                vtj1(i,j,k1)=linear1d(y1(0,0,0),  y1(0,1,0),           &
                                     vti1(i,0,k1),vti1(i,1,k1),        &
                                     y2(0,j,0))
                vtj2(i,j,k1)=linear1d(y1(0,0,0),  y1(0,1,0),           &
                                     vti2(i,0,k1),vti2(i,1,k1),        &
                                     y2(0,j,0))
                vtj3(i,j,k1)=linear1d(y1(0,0,0),  y1(0,1,0),           &
                                     vti3(i,0,k1),vti3(i,1,k1),        &
                                     y2(0,j,0))
                vtj4(i,j,k1)=linear1d(y1(0,0,0),  y1(0,1,0),           &
                                     vti4(i,0,k1),vti4(i,1,k1),        &
                                     y2(0,j,0))
                vtj5(i,j,k1)=linear1d(y1(0,0,0),  y1(0,1,0),           &
                                     vti5(i,0,k1),vti5(i,1,k1),        &
                                     y2(0,j,0))
                vtj6(i,j,k1)=linear1d(y1(0,0,0),  y1(0,1,0),           &
                                     vti6(i,0,k1),vti6(i,1,k1),        &
                                     y2(0,j,0))
              end do
              end do
              !
              exit
              !
    			!
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    !
    do k=0,dkm2
    	!
      za=z2(0,0,k)
      !
      if(za>z1(0,0,dkm1))  za=za-z1(0,0,dkm1) ! for periodic condition
      !
    	do k1=1,dkm1
    		!
    		if( za>=z1(0,0,k1-1) .and. za<=z1(0,0,k1) ) then
    		  do j=0,djm2
    		  do i=0,dim2
    		    v21(i,j,k)=linear1d(z1(0,0,k1-1),     z1(0,0,k1),          &
    		                         vtj1(i,j,k1-1),vtj1(i,j,k1),          &
    		                         za)
    		    v22(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj2(i,j,k1-1),vtj2(i,j,k1),          &
    		                         za)
    		    v23(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj3(i,j,k1-1),vtj3(i,j,k1),          &
    		                         za)
    		    v24(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj4(i,j,k1-1),vtj4(i,j,k1),          &
    		                         za)
    		    v25(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj5(i,j,k1-1),vtj5(i,j,k1),          &
    		                         za)
    		    v26(i,j,k)=linear1d(z1(0,0,k1-1),   z1(0,0,k1),            &
    		                         vtj6(i,j,k1-1),vtj6(i,j,k1),          &
    		                         za)
    		  end do
    		  end do
    		  !
    		  exit
    		  !
            elseif( za>=z1(0,0,dkm1) ) then
              do j=0,djm2
              do i=0,dim2
                v21(i,j,k)=linear1d(z1(0,0,dkm1-1),     z1(0,0,dkm1),          &
                                     vtj1(i,j,dkm1-1),vtj1(i,j,dkm1),          &
                                     za)
                v22(i,j,k)=linear1d(z1(0,0,dkm1-1),   z1(0,0,dkm1),            &
                                     vtj2(i,j,dkm1-1),vtj2(i,j,dkm1),          &
                                     za)
                v23(i,j,k)=linear1d(z1(0,0,dkm1-1),   z1(0,0,dkm1),            &
                                     vtj3(i,j,dkm1-1),vtj3(i,j,dkm1),          &
                                     za)
                v24(i,j,k)=linear1d(z1(0,0,dkm1-1),   z1(0,0,dkm1),            &
                                     vtj4(i,j,dkm1-1),vtj4(i,j,dkm1),          &
                                     za)
                v25(i,j,k)=linear1d(z1(0,0,dkm1-1),   z1(0,0,dkm1),            &
                                     vtj5(i,j,dkm1-1),vtj5(i,j,dkm1),          &
                                     za)
                v26(i,j,k)=linear1d(z1(0,0,dkm1-1),   z1(0,0,dkm1),            &
                                     vtj6(i,j,dkm1-1),vtj6(i,j,dkm1),          &
                                     za)
              end do
              end do
              !
              exit
              !
            elseif( za<=z1(0,0,0) ) then
              do j=0,djm2
              do i=0,dim2
                v21(i,j,k)=linear1d(z1(0,0,0),     z1(0,0,1),          &
                                     vtj1(i,j,0),vtj1(i,j,1),          &
                                     za)
                v22(i,j,k)=linear1d(z1(0,0,0),   z1(0,0,1),            &
                                     vtj2(i,j,0),vtj2(i,j,1),          &
                                     za)
                v23(i,j,k)=linear1d(z1(0,0,0),   z1(0,0,1),            &
                                     vtj3(i,j,0),vtj3(i,j,1),          &
                                     za)
                v24(i,j,k)=linear1d(z1(0,0,0),   z1(0,0,1),            &
                                     vtj4(i,j,0),vtj4(i,j,1),          &
                                     za)
                v25(i,j,k)=linear1d(z1(0,0,0),   z1(0,0,1),            &
                                     vtj5(i,j,0),vtj5(i,j,1),          &
                                     za)
                v26(i,j,k)=linear1d(z1(0,0,0),   z1(0,0,1),            &
                                     vtj6(i,j,0),vtj6(i,j,1),          &
                                     za)
              end do
              end do
              !
              exit
              !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	          '  ** Interpolating ... ',100*nper/(dim2+djm2+dkm2+3),' %'
    end do
    print*
    !
  end subroutine regularlinearinterp3d6v
  !+-------------------------------------------------------------------+
  !| The end of the subroutine regularlinearinterp3d                   |
  !+-------------------------------------------------------------------+
  !!
  subroutine regularlinearinterp3d6v_sp(x1,y1,z1,                      &
                                        v11,v12,v13,v14,v15,v16,       &
                                        dim1,djm1,dkm1,                &
                                        x2,y2,z2,                      &
                                        v21,v22,v23,v24,v25,v26)
    !
    integer,intent(in) :: dim1,djm1,dkm1
    !
    real(8),intent(in) :: x1(0:dim1,0:djm1,0:dkm1),                    &
                          y1(0:dim1,0:djm1,0:dkm1),                    &
                          z1(0:dim1,0:djm1,0:dkm1),                    &
                          v11(0:dim1,0:djm1,0:dkm1),                   &
                          v12(0:dim1,0:djm1,0:dkm1),                   &
                          v13(0:dim1,0:djm1,0:dkm1),                   &
                          v14(0:dim1,0:djm1,0:dkm1),                   &
                          v15(0:dim1,0:djm1,0:dkm1),                   &
                          v16(0:dim1,0:djm1,0:dkm1)
    real(8),intent(in) :: x2,y2,z2
    real(8),intent(out) :: v21,v22,v23,v24,v25,v26
    !
    integer :: i1,j1,k1
    !
    real(8) :: vti1(0:djm1,0:dkm1),vtj1(0:dkm1),  &
               vti2(0:djm1,0:dkm1),vtj2(0:dkm1),  &
               vti3(0:djm1,0:dkm1),vtj3(0:dkm1),  &
               vti4(0:djm1,0:dkm1),vtj4(0:dkm1),  &
               vti5(0:djm1,0:dkm1),vtj5(0:dkm1),  &
               vti6(0:djm1,0:dkm1),vtj6(0:dkm1)
    !
    do i1=1,dim1
      !
      if( x2>=x1(i1-1,0,0) .and. x2<=x1(i1,0,0) ) then
        !
        do k1=0,dkm1
        do j1=0,djm1
          vti1(j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),            &
                               v11(i1-1,j1,k1),v11(i1,j1,k1),x2      )
          vti2(j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),            &
                               v12(i1-1,j1,k1),v12(i1,j1,k1),x2      )
          vti3(j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),            &
                               v13(i1-1,j1,k1),v13(i1,j1,k1),x2      )
          vti4(j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),            &
                               v14(i1-1,j1,k1),v14(i1,j1,k1),x2      )
          vti5(j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),            &
                               v15(i1-1,j1,k1),v15(i1,j1,k1),x2      )
          vti6(j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),            &
                               v16(i1-1,j1,k1),v16(i1,j1,k1),x2      )
        end do
        end do
        !
        exit
        !
      elseif( x2>=x1(dim1,0,0) ) then
        !
        do k1=0,dkm1
        do j1=0,djm1
          vti1(j1,k1)=linear1d(x1(dim1-1,0,0),   x1(dim1,0,0),         &
                               v11(dim1-1,j1,k1),v11(dim1,j1,k1),x2)
          vti2(j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                                v12(dim1-1,j1,k1),v12(dim1,j1,k1),x2)
          vti3(j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                                v13(dim1-1,j1,k1),v13(dim1,j1,k1),x2)
          vti4(j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                                v14(dim1-1,j1,k1),v14(dim1,j1,k1),x2)
          vti5(j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                                v15(dim1-1,j1,k1),v15(dim1,j1,k1),x2)
          vti6(j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                                v16(dim1-1,j1,k1),v16(dim1,j1,k1),x2)
        end do
        end do
        !
        exit
        !
      elseif( x2<=x1(0,0,0) ) then
        !
        do k1=0,dkm1
        do j1=0,djm1
          vti1(j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                               v11(0,j1,k1),v11(1,j1,k1),x2)
          vti2(j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                                v12(0,j1,k1),v12(1,j1,k1),x2)
          vti3(j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                                v13(0,j1,k1),v13(1,j1,k1),x2)
          vti4(j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                                v14(0,j1,k1),v14(1,j1,k1),x2)
          vti5(j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                                v15(0,j1,k1),v15(1,j1,k1),x2)
          vti6(j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                                v16(0,j1,k1),v16(1,j1,k1),x2)
        end do
        end do
        !
        exit
        !
      end if
      !
    end do
    !
    do j1=1,djm1
      !
      if( y2>=y1(0,j1-1,0) .and. y2<=y1(0,j1,0) ) then
        do k1=0,dkm1
          vtj1(k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
                               vti1(j1-1,k1),vti1(j1,k1),y2)
          vtj2(k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
                               vti2(j1-1,k1),vti2(j1,k1),y2)
          vtj3(k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
                               vti3(j1-1,k1),vti3(j1,k1),y2)
          vtj4(k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
                               vti4(j1-1,k1),vti4(j1,k1),y2)
          vtj5(k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
                               vti5(j1-1,k1),vti5(j1,k1),y2)
          vtj6(k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
                               vti6(j1-1,k1),vti6(j1,k1),y2)
        end do
        !
        exit
        !
      elseif( y2>=y1(0,djm1,0) ) then
        do k1=0,dkm1
          vtj1(k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                               vti1(djm1-1,k1),vti1(djm1,k1),y2)
          vtj2(k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                               vti2(djm1-1,k1),vti2(djm1,k1),y2)
          vtj3(k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                               vti3(djm1-1,k1),vti3(djm1,k1),y2)
          vtj4(k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                               vti4(djm1-1,k1),vti4(djm1,k1),y2)
          vtj5(k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                               vti5(djm1-1,k1),vti5(djm1,k1),y2)
          vtj6(k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                               vti6(djm1-1,k1),vti6(djm1,k1),y2)
        end do
        !
        exit
        !
      elseif( y2<=y1(0,0,0) ) then
        do k1=0,dkm1
          vtj1(k1)=linear1d(y1(0,0,0),y1(0,1,0),vti1(0,k1),vti1(1,k1),y2)
          vtj2(k1)=linear1d(y1(0,0,0),y1(0,1,0),vti2(0,k1),vti2(1,k1),y2)
          vtj3(k1)=linear1d(y1(0,0,0),y1(0,1,0),vti3(0,k1),vti3(1,k1),y2)
          vtj4(k1)=linear1d(y1(0,0,0),y1(0,1,0),vti4(0,k1),vti4(1,k1),y2)
          vtj5(k1)=linear1d(y1(0,0,0),y1(0,1,0),vti5(0,k1),vti5(1,k1),y2)
          vtj6(k1)=linear1d(y1(0,0,0),y1(0,1,0),vti6(0,k1),vti6(1,k1),y2)
        end do
        !
        exit
        !
      end if
      !
    end do
    !
    do k1=1,dkm1
      !
      if( z2>=z1(0,0,k1-1) .and. z2<=z1(0,0,k1) ) then
        !
        v21=linear1d(z1(0,0,k1-1),z1(0,0,k1),vtj1(k1-1),vtj1(k1),z2)
        v22=linear1d(z1(0,0,k1-1),z1(0,0,k1),vtj2(k1-1),vtj2(k1),z2)
        v23=linear1d(z1(0,0,k1-1),z1(0,0,k1),vtj3(k1-1),vtj3(k1),z2)
        v24=linear1d(z1(0,0,k1-1),z1(0,0,k1),vtj4(k1-1),vtj4(k1),z2)
        v25=linear1d(z1(0,0,k1-1),z1(0,0,k1),vtj5(k1-1),vtj5(k1),z2)
        v26=linear1d(z1(0,0,k1-1),z1(0,0,k1),vtj6(k1-1),vtj6(k1),z2)
      elseif( z2>=z1(0,0,dkm1) ) then
        v21=linear1d(z1(0,0,dkm1-1),z1(0,0,dkm1),vtj1(dkm1-1),vtj1(dkm1),z2)
        v22=linear1d(z1(0,0,dkm1-1),z1(0,0,dkm1),vtj2(dkm1-1),vtj2(dkm1),z2)
        v23=linear1d(z1(0,0,dkm1-1),z1(0,0,dkm1),vtj3(dkm1-1),vtj3(dkm1),z2)
        v24=linear1d(z1(0,0,dkm1-1),z1(0,0,dkm1),vtj4(dkm1-1),vtj4(dkm1),z2)
        v25=linear1d(z1(0,0,dkm1-1),z1(0,0,dkm1),vtj5(dkm1-1),vtj5(dkm1),z2)
        v26=linear1d(z1(0,0,dkm1-1),z1(0,0,dkm1),vtj6(dkm1-1),vtj6(dkm1),z2)
      elseif( z2<=z1(0,0,0) ) then
        v21=linear1d(z1(0,0,0),z1(0,0,1),vtj1(0),vtj1(1),z2)
        v22=linear1d(z1(0,0,0),z1(0,0,1),vtj2(0),vtj2(1),z2)
        v23=linear1d(z1(0,0,0),z1(0,0,1),vtj3(0),vtj3(1),z2)
        v24=linear1d(z1(0,0,0),z1(0,0,1),vtj4(0),vtj4(1),z2)
        v25=linear1d(z1(0,0,0),z1(0,0,1),vtj5(0),vtj5(1),z2)
        v26=linear1d(z1(0,0,0),z1(0,0,1),vtj6(0),vtj6(1),z2)
      end if
      !
    enddo
    write(*,'(A)')'  ** Interpolating ...  100 %'
    !
  end subroutine regularlinearinterp3d6v_sp

  subroutine regularlinearinterp3d7v_sp(x1,y1,z1,                      &
                                        v11,v12,v13,v14,v15,v16,v17,   &
                                        dim1,djm1,dkm1,                &
                                        x2,y2,z2,                      &
                                        v21,v22,v23,v24,v25,v26,v27)
    !
    integer,intent(in) :: dim1,djm1,dkm1
    !
    real(8),intent(in) :: x1(0:dim1,0:djm1,0:dkm1),                    &
                          y1(0:dim1,0:djm1,0:dkm1),                    &
                          z1(0:dim1,0:djm1,0:dkm1),                    &
                          v11(0:dim1,0:djm1,0:dkm1),                   &
                          v12(0:dim1,0:djm1,0:dkm1),                   &
                          v13(0:dim1,0:djm1,0:dkm1),                   &
                          v14(0:dim1,0:djm1,0:dkm1),                   &
                          v15(0:dim1,0:djm1,0:dkm1),                   &
                          v16(0:dim1,0:djm1,0:dkm1),                   &
                          v17(0:dim1,0:djm1,0:dkm1)
    real(8),intent(in) :: x2,y2,z2
    real(8),intent(out) :: v21,v22,v23,v24,v25,v26,v27
    !
    integer :: i1,j1,k1
    !
    real(8) :: vti1(0:djm1,0:dkm1),vtj1(0:dkm1),  &
               vti2(0:djm1,0:dkm1),vtj2(0:dkm1),  &
               vti3(0:djm1,0:dkm1),vtj3(0:dkm1),  &
               vti4(0:djm1,0:dkm1),vtj4(0:dkm1),  &
               vti5(0:djm1,0:dkm1),vtj5(0:dkm1),  &
               vti6(0:djm1,0:dkm1),vtj6(0:dkm1),  &
               vti7(0:djm1,0:dkm1),vtj7(0:dkm1)
    !
    do i1=1,dim1
      !
      if( x2>=x1(i1-1,0,0) .and. x2<=x1(i1,0,0) ) then
        !
        do k1=0,dkm1
        do j1=0,djm1
          vti1(j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),            &
                               v11(i1-1,j1,k1),v11(i1,j1,k1),x2      )
          vti2(j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),            &
                               v12(i1-1,j1,k1),v12(i1,j1,k1),x2      )
          vti3(j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),            &
                               v13(i1-1,j1,k1),v13(i1,j1,k1),x2      )
          vti4(j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),            &
                               v14(i1-1,j1,k1),v14(i1,j1,k1),x2      )
          vti5(j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),            &
                               v15(i1-1,j1,k1),v15(i1,j1,k1),x2      )
          vti6(j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),            &
                               v16(i1-1,j1,k1),v16(i1,j1,k1),x2      )
          vti7(j1,k1)=linear1d(x1(i1-1,0,0),  x1(i1,0,0),            &
                               v17(i1-1,j1,k1),v17(i1,j1,k1),x2      )
        end do
        end do
        !
        exit
        !
      elseif( x2>=x1(dim1,0,0) ) then
        !
        do k1=0,dkm1
        do j1=0,djm1
          vti1(j1,k1)=linear1d(x1(dim1-1,0,0),   x1(dim1,0,0),         &
                               v11(dim1-1,j1,k1),v11(dim1,j1,k1),x2)
          vti2(j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                                v12(dim1-1,j1,k1),v12(dim1,j1,k1),x2)
          vti3(j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                                v13(dim1-1,j1,k1),v13(dim1,j1,k1),x2)
          vti4(j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                                v14(dim1-1,j1,k1),v14(dim1,j1,k1),x2)
          vti5(j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                                v15(dim1-1,j1,k1),v15(dim1,j1,k1),x2)
          vti6(j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                                v16(dim1-1,j1,k1),v16(dim1,j1,k1),x2)
          vti7(j1,k1)=linear1d(x1(dim1-1,0,0),  x1(dim1,0,0),          &
                                v17(dim1-1,j1,k1),v17(dim1,j1,k1),x2)
        end do
        end do
        !
        exit
        !
      elseif( x2<=x1(0,0,0) ) then
        !
        do k1=0,dkm1
        do j1=0,djm1
          vti1(j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                               v11(0,j1,k1),v11(1,j1,k1),x2)
          vti2(j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                                v12(0,j1,k1),v12(1,j1,k1),x2)
          vti3(j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                                v13(0,j1,k1),v13(1,j1,k1),x2)
          vti4(j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                                v14(0,j1,k1),v14(1,j1,k1),x2)
          vti5(j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                                v15(0,j1,k1),v15(1,j1,k1),x2)
          vti6(j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                                v16(0,j1,k1),v16(1,j1,k1),x2)
          vti7(j1,k1)=linear1d(x1(0,0,0),  x1(1,0,0),          &
                                v17(0,j1,k1),v17(1,j1,k1),x2)
        end do
        end do
        !
        exit
        !
      end if
      !
    end do
    !
    do j1=1,djm1
      !
      if( y2>=y1(0,j1-1,0) .and. y2<=y1(0,j1,0) ) then
        do k1=0,dkm1
          vtj1(k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
                               vti1(j1-1,k1),vti1(j1,k1),y2)
          vtj2(k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
                               vti2(j1-1,k1),vti2(j1,k1),y2)
          vtj3(k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
                               vti3(j1-1,k1),vti3(j1,k1),y2)
          vtj4(k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
                               vti4(j1-1,k1),vti4(j1,k1),y2)
          vtj5(k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
                               vti5(j1-1,k1),vti5(j1,k1),y2)
          vtj6(k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
                               vti6(j1-1,k1),vti6(j1,k1),y2)
          vtj7(k1)=linear1d(y1(0,j1-1,0),  y1(0,j1,0),           &
                               vti7(j1-1,k1),vti7(j1,k1),y2)
        end do
        !
        exit
        !
      elseif( y2>=y1(0,djm1,0) ) then
        do k1=0,dkm1
          vtj1(k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                               vti1(djm1-1,k1),vti1(djm1,k1),y2)
          vtj2(k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                               vti2(djm1-1,k1),vti2(djm1,k1),y2)
          vtj3(k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                               vti3(djm1-1,k1),vti3(djm1,k1),y2)
          vtj4(k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                               vti4(djm1-1,k1),vti4(djm1,k1),y2)
          vtj5(k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                               vti5(djm1-1,k1),vti5(djm1,k1),y2)
          vtj6(k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                               vti6(djm1-1,k1),vti6(djm1,k1),y2)
          vtj7(k1)=linear1d(y1(0,djm1-1,0),  y1(0,djm1,0),           &
                               vti7(djm1-1,k1),vti7(djm1,k1),y2)
        end do
        !
        exit
        !
      elseif( y2<=y1(0,0,0) ) then
        do k1=0,dkm1
          vtj1(k1)=linear1d(y1(0,0,0),y1(0,1,0),vti1(0,k1),vti1(1,k1),y2)
          vtj2(k1)=linear1d(y1(0,0,0),y1(0,1,0),vti2(0,k1),vti2(1,k1),y2)
          vtj3(k1)=linear1d(y1(0,0,0),y1(0,1,0),vti3(0,k1),vti3(1,k1),y2)
          vtj4(k1)=linear1d(y1(0,0,0),y1(0,1,0),vti4(0,k1),vti4(1,k1),y2)
          vtj5(k1)=linear1d(y1(0,0,0),y1(0,1,0),vti5(0,k1),vti5(1,k1),y2)
          vtj6(k1)=linear1d(y1(0,0,0),y1(0,1,0),vti6(0,k1),vti6(1,k1),y2)
          vtj7(k1)=linear1d(y1(0,0,0),y1(0,1,0),vti7(0,k1),vti7(1,k1),y2)
        end do
        !
        exit
        !
      end if
      !
    end do
    !
    do k1=1,dkm1
      !
      if( z2>=z1(0,0,k1-1) .and. z2<=z1(0,0,k1) ) then
        !
        v21=linear1d(z1(0,0,k1-1),z1(0,0,k1),vtj1(k1-1),vtj1(k1),z2)
        v22=linear1d(z1(0,0,k1-1),z1(0,0,k1),vtj2(k1-1),vtj2(k1),z2)
        v23=linear1d(z1(0,0,k1-1),z1(0,0,k1),vtj3(k1-1),vtj3(k1),z2)
        v24=linear1d(z1(0,0,k1-1),z1(0,0,k1),vtj4(k1-1),vtj4(k1),z2)
        v25=linear1d(z1(0,0,k1-1),z1(0,0,k1),vtj5(k1-1),vtj5(k1),z2)
        v26=linear1d(z1(0,0,k1-1),z1(0,0,k1),vtj6(k1-1),vtj6(k1),z2)
        v27=linear1d(z1(0,0,k1-1),z1(0,0,k1),vtj7(k1-1),vtj7(k1),z2)
      elseif( z2>=z1(0,0,dkm1) ) then
        v21=linear1d(z1(0,0,dkm1-1),z1(0,0,dkm1),vtj1(dkm1-1),vtj1(dkm1),z2)
        v22=linear1d(z1(0,0,dkm1-1),z1(0,0,dkm1),vtj2(dkm1-1),vtj2(dkm1),z2)
        v23=linear1d(z1(0,0,dkm1-1),z1(0,0,dkm1),vtj3(dkm1-1),vtj3(dkm1),z2)
        v24=linear1d(z1(0,0,dkm1-1),z1(0,0,dkm1),vtj4(dkm1-1),vtj4(dkm1),z2)
        v25=linear1d(z1(0,0,dkm1-1),z1(0,0,dkm1),vtj5(dkm1-1),vtj5(dkm1),z2)
        v26=linear1d(z1(0,0,dkm1-1),z1(0,0,dkm1),vtj6(dkm1-1),vtj6(dkm1),z2)
        v27=linear1d(z1(0,0,dkm1-1),z1(0,0,dkm1),vtj7(dkm1-1),vtj7(dkm1),z2)
      elseif( z2<=z1(0,0,0) ) then
        v21=linear1d(z1(0,0,0),z1(0,0,1),vtj1(0),vtj1(1),z2)
        v22=linear1d(z1(0,0,0),z1(0,0,1),vtj2(0),vtj2(1),z2)
        v23=linear1d(z1(0,0,0),z1(0,0,1),vtj3(0),vtj3(1),z2)
        v24=linear1d(z1(0,0,0),z1(0,0,1),vtj4(0),vtj4(1),z2)
        v25=linear1d(z1(0,0,0),z1(0,0,1),vtj5(0),vtj5(1),z2)
        v26=linear1d(z1(0,0,0),z1(0,0,1),vtj6(0),vtj6(1),z2)
        v27=linear1d(z1(0,0,0),z1(0,0,1),vtj7(0),vtj7(1),z2)
      end if
      !
    enddo
    write(*,'(A)')'  ** Interpolating ...  100 %'
    !
  end subroutine regularlinearinterp3d7v_sp
  !+-------------------------------------------------------------------+
  !| The end of the subroutine regularlinearinterp3d6v_sp              |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to interpolate 2D flowfield between 2     |
  !| regular meshes with linear approach.                              |
  !+-------------------------------------------------------------------+
  subroutine regularlinearinterp2d1v(x1,y1,v11,dim1,djm1,               &
                                     x2,y2,v21,dim2,djm2)
    !
    integer,intent(in) :: dim1,djm1,dim2,djm2
    real(8),intent(in) :: x1(0:dim1,0:djm1),y1(0:dim1,0:djm1),         &
                          v11(0:dim1,0:djm1)
    real(8),intent(in) :: x2(0:dim2,0:djm2),y2(0:dim2,0:djm2)
    real(8),intent(out) :: v21(0:dim2,0:djm2)
    integer :: i,j,i1,j1,nper
    !
    real(8) :: vti(0:dim2,0:djm1)
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i,0)>=x1(i1-1,0) .and. x2(i,0)<=x1(i1,0) ) then
    		  !
    		  do j1=0,djm1
    		    vti(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                   &
    		                       v11(i1-1,j1),v11(i1,j1),x2(i,0))
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	               '  ** Interpolating ... ',100*nper/(dim2+djm2+2),' %'
    end do
    !
    do j=0,djm2
    	!
    	do j1=1,djm1
    		!
    		if( y2(0,j)>=y1(0,j1-1) .and. y2(0,j)<=y1(0,j1) ) then
    		  do i=0,dim2
    		    v21(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti(i,j1-1),vti(i,j1),y2(0,j))
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	               '  ** Interpolating ... ',100*nper/(dim2+djm2+2),' %'
    end do
    print*
    !
  end subroutine regularlinearinterp2d1v
  !!
  subroutine regularlinearinterp2d2v(x1,y1,v11,v12,dim1,djm1,          &
                                     x2,y2,v21,v22,dim2,djm2)
    !
    integer,intent(in) :: dim1,djm1,dim2,djm2
    real(8),intent(in) :: x1(0:dim1,0:djm1),y1(0:dim1,0:djm1),         &
                          v11(0:dim1,0:djm1),v12(0:dim1,0:djm1)
    real(8),intent(in) :: x2(0:dim2,0:djm2),y2(0:dim2,0:djm2)
    real(8),intent(out) :: v21(0:dim2,0:djm2),v22(0:dim2,0:djm2)
    integer :: i,j,i1,j1,nper
    !
    real(8) :: vti1(0:dim2,0:djm1),vti2(0:dim2,0:djm1)
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i,0)>=x1(i1-1,0) .and. x2(i,0)<=x1(i1,0) ) then
    		  !
    		  do j1=0,djm1
    		    vti1(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v11(i1-1,j1),v11(i1,j1),x2(i,0))
    		    vti2(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v12(i1-1,j1),v12(i1,j1),x2(i,0))
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	               '  ** Interpolating ... ',100*nper/(dim2+djm2+2),' %'
    end do
    !
    do j=0,djm2
    	!
    	do j1=1,djm1
    		!
    		if( y2(0,j)>=y1(0,j1-1) .and. y2(0,j)<=y1(0,j1) ) then
    		  do i=0,dim2
    		    v21(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti1(i,j1-1),vti1(i,j1),y2(0,j))
    		    v22(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti2(i,j1-1),vti2(i,j1),y2(0,j))
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	               '  ** Interpolating ... ',100*nper/(dim2+djm2+2),' %'
    end do
    print*
    !
  end subroutine regularlinearinterp2d2v
  !!
  subroutine regularlinearinterp2d3v(x1,y1,v11,v12,v13,dim1,djm1,      &
                                     x2,y2,v21,v22,v23,dim2,djm2)
    !
    integer,intent(in) :: dim1,djm1,dim2,djm2
    real(8),intent(in) :: x1(0:dim1,0:djm1),y1(0:dim1,0:djm1),         &
                          v11(0:dim1,0:djm1),v12(0:dim1,0:djm1),       &
                          v13(0:dim1,0:djm1)
    real(8),intent(in) :: x2(0:dim2,0:djm2),y2(0:dim2,0:djm2)
    real(8),intent(out) :: v21(0:dim2,0:djm2),v22(0:dim2,0:djm2),      &
                           v23(0:dim2,0:djm2)
    integer :: i,j,i1,j1,nper
    !
    real(8) :: vti1(0:dim2,0:djm1),vti2(0:dim2,0:djm1),                &
               vti3(0:dim2,0:djm1)
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i,0)>=x1(i1-1,0) .and. x2(i,0)<=x1(i1,0) ) then
    		  !
    		  do j1=0,djm1
    		    vti1(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v11(i1-1,j1),v11(i1,j1),x2(i,0))
    		    vti2(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v12(i1-1,j1),v12(i1,j1),x2(i,0))
    		    vti3(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v13(i1-1,j1),v13(i1,j1),x2(i,0))
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	               '  ** Interpolating ... ',100*nper/(dim2+djm2+2),' %'
    end do
    !
    do j=0,djm2
    	!
    	do j1=1,djm1
    		!
    		if( y2(0,j)>=y1(0,j1-1) .and. y2(0,j)<=y1(0,j1) ) then
    		  do i=0,dim2
    		    v21(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti1(i,j1-1),vti1(i,j1),y2(0,j))
    		    v22(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti2(i,j1-1),vti2(i,j1),y2(0,j))
    		    v23(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti3(i,j1-1),vti3(i,j1),y2(0,j))
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	               '  ** Interpolating ... ',100*nper/(dim2+djm2+2),' %'
    end do
    print*
    !
  end subroutine regularlinearinterp2d3v
  !
  subroutine regularlinearinterp2d4v(x1,y1,v11,v12,v13,v14,dim1,djm1,  &
                                     x2,y2,v21,v22,v23,v24,dim2,djm2)
    !
    integer,intent(in) :: dim1,djm1,dim2,djm2
    real(8),intent(in) :: x1(0:dim1,0:djm1),y1(0:dim1,0:djm1),         &
                          v11(0:dim1,0:djm1),v12(0:dim1,0:djm1),       &
                          v13(0:dim1,0:djm1),v14(0:dim1,0:djm1)
    real(8),intent(in) :: x2(0:dim2,0:djm2),y2(0:dim2,0:djm2)
    real(8),intent(out) :: v21(0:dim2,0:djm2),v22(0:dim2,0:djm2),      &
                           v23(0:dim2,0:djm2),v24(0:dim2,0:djm2)
    integer :: i,j,i1,j1,nper
    !
    real(8) :: vti1(0:dim2,0:djm1),vti2(0:dim2,0:djm1),                &
               vti3(0:dim2,0:djm1),vti4(0:dim2,0:djm1)
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i,0)>=x1(i1-1,0) .and. x2(i,0)<=x1(i1,0) ) then
    		  !
    		  do j1=0,djm1
    		    vti1(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v11(i1-1,j1),v11(i1,j1),x2(i,0))
    		    vti2(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v12(i1-1,j1),v12(i1,j1),x2(i,0))
    		    vti3(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v13(i1-1,j1),v13(i1,j1),x2(i,0))
    		    vti4(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v14(i1-1,j1),v14(i1,j1),x2(i,0))
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	               '  ** Interpolating ... ',100*nper/(dim2+djm2+2),' %'
    end do
    !
    do j=0,djm2
    	!
    	do j1=1,djm1
    		!
    		if( y2(0,j)>=y1(0,j1-1) .and. y2(0,j)<=y1(0,j1) ) then
    		  do i=0,dim2
    		    v21(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti1(i,j1-1),vti1(i,j1),y2(0,j))
    		    v22(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti2(i,j1-1),vti2(i,j1),y2(0,j))
    		    v23(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti3(i,j1-1),vti3(i,j1),y2(0,j))
    		    v24(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti4(i,j1-1),vti4(i,j1),y2(0,j))
    		  end do
    		  !
    		  exit
    		  !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	               '  ** Interpolating ... ',100*nper/(dim2+djm2+2),' %'
    end do
    print*
    !
  end subroutine regularlinearinterp2d4v
  !!
  subroutine regularlinearinterp2d5v(x1,y1,v11,v12,v13,v14,v15,        &
                                                     dim1,djm1,        &
                                     x2,y2,v21,v22,v23,v24,v25,        &
                                                     dim2,djm2,        &
                                                   homox,homoy,progress)
    !
    integer,intent(in) :: dim1,djm1,dim2,djm2
    real(8),intent(in) :: x1(0:dim1,0:djm1),y1(0:dim1,0:djm1),         &
                          v11(0:dim1,0:djm1),v12(0:dim1,0:djm1),       &
                          v13(0:dim1,0:djm1),v14(0:dim1,0:djm1),       &
                          v15(0:dim1,0:djm1)
    real(8),intent(in) :: x2(0:dim2,0:djm2),y2(0:dim2,0:djm2)
    logical,intent(in),optional :: homox,homoy,progress
    real(8),intent(out) :: v21(0:dim2,0:djm2),v22(0:dim2,0:djm2),      &
                           v23(0:dim2,0:djm2),v24(0:dim2,0:djm2),      &
                           v25(0:dim2,0:djm2)
    integer :: i,j,i1,j1,nper
    logical :: lihomo,ljhomo,lprog
    !
    real(8) :: vti1(0:dim2,0:djm1),vti2(0:dim2,0:djm1),                &
               vti3(0:dim2,0:djm1),vti4(0:dim2,0:djm1),                &
               vti5(0:dim2,0:djm1)
    real(8) :: dist,damp
    real(8),allocatable,dimension(:) :: xi,yi,xt,yt
    !
    if(present(homox)) then
      lihomo=homox
    else
      lihomo=.false.
    endif
    if(present(homoy)) then
      ljhomo=homoy
    else
      ljhomo=.false.
    endif
    if(present(progress)) then
      lprog=progress
    else
      lprog=.false.
    endif
    !
    allocate(xi(0:dim1),yi(0:djm1),xt(0:dim2),yt(0:djm2))
    !
    xi=x1(:,0)
    yi=y1(0,:)
    xt=x2(:,0)
    yt=y2(0,:)
    !
    if(lihomo) then
      do while(xt(dim2)>xi(dim1))
        do i=0,dim2
          if(xt(i)>xi(dim1)) then
            xt(i)=xt(i)-xi(dim1)
          endif
        enddo
      enddo
    endif
    if(ljhomo) then
      do while(yt(djm2)>yi(djm1))
        do j=0,djm2
          if(yt(j)>yi(djm1)) then
            yt(j)=yt(j)-yi(djm1)
          endif
        enddo
      enddo
    endif
    !
    nper=0
    do i=0,dim2
      !
    	do i1=1,dim1
        !
    		if( xt(i)>=xi(i1-1) .and. xt(i)<=xi(i1) ) then
          !
    		  do j1=0,djm1
            vti1(i,j1)=linear1d(xi(i1-1), xi(i1),                  &
                               v11(i1-1,j1),v11(i1,j1),xt(i))
            vti2(i,j1)=linear1d(xi(i1-1), xi(i1),                  &
                               v12(i1-1,j1),v12(i1,j1),xt(i))
            vti3(i,j1)=linear1d(xi(i1-1), xi(i1),                  &
                               v13(i1-1,j1),v13(i1,j1),xt(i))
            vti4(i,j1)=linear1d(xi(i1-1), xi(i1),                  &
                               v14(i1-1,j1),v14(i1,j1),xt(i))
            vti5(i,j1)=linear1d(xi(i1-1), xi(i1),                  &
                               v15(i1-1,j1),v15(i1,j1),xt(i))
    		  end do
          !
          exit
          !
        elseif(xt(i)>xi(dim1)) then
          !
          do j1=0,djm1
            !
            dist=xt(i)-xi(dim1)
            damp=dexp(-dist**2)
            !
            vti1(i,j1)=v11(dim1,j1)*damp
            vti2(i,j1)=v12(dim1,j1)*damp
            vti3(i,j1)=v13(dim1,j1)*damp
            vti4(i,j1)=v14(dim1,j1)*damp
            vti5(i,j1)=v15(dim1,j1)*damp
            !
          end do
          !
        elseif(xt(i)<=xi(0)) then
          !
          do j1=0,djm1
            !
            vti1(i,j1)=v11(0,j1)
            vti2(i,j1)=v12(0,j1)
            vti3(i,j1)=v13(0,j1)
            vti4(i,j1)=v14(0,j1)
            vti5(i,j1)=v15(0,j1)
            !
          end do
          !
    		end if
        !
    	end do
      !
      nper=nper+1
      !
      if(lprog) write(*,'(1A1,A23,I4,A2,$)')char(13),                   &
                     '  ** Interpolating ... ',100*nper/(dim2+djm2+2),' %'
    end do
    !
    do j=0,djm2
      !
    	do j1=1,djm1
        !
    		if( yt(j)>=yi(j1-1) .and. yt(j)<=yi(j1) ) then
    		  do i=0,dim2
            v21(i,j)=linear1d(yi(j1-1),  yi(j1),                   &
                             vti1(i,j1-1),vti1(i,j1),yt(j))
            v22(i,j)=linear1d(yi(j1-1),  yi(j1),                   &
                             vti2(i,j1-1),vti2(i,j1),yt(j))
            v23(i,j)=linear1d(yi(j1-1),  yi(j1),                   &
                             vti3(i,j1-1),vti3(i,j1),yt(j))
            v24(i,j)=linear1d(yi(j1-1),  yi(j1),                   &
                             vti4(i,j1-1),vti4(i,j1),yt(j))
            v25(i,j)=linear1d(yi(j1-1),  yi(j1),                   &
                             vti5(i,j1-1),vti5(i,j1),yt(j))
    		  end do
          !
          exit
          !
        elseif(yt(j)<=yi(0)) then
          do i=0,dim2
            v21(i,j)=vti1(i,0)
            v22(i,j)=vti2(i,0)
            v23(i,j)=vti3(i,0)
            v24(i,j)=vti4(i,0)
            v25(i,j)=vti5(i,0)
          end do
          !
          exit
          !
        elseif(yt(j)>=yi(djm1)) then
          do i=0,dim2
            v21(i,j)=vti1(i,djm1)
            v22(i,j)=vti2(i,djm1)
            v23(i,j)=vti3(i,djm1)
            v24(i,j)=vti4(i,djm1)
            v25(i,j)=vti5(i,djm1)
          end do
          !
          exit
          !
    		end if
        !
    	end do
      !
      nper=nper+1
      !
      if(lprog) write(*,'(1A1,A23,I4,A2,$)')char(13),                   &
                     '  ** Interpolating ... ',200*nper/(dim2+djm2+2),' %'
    end do
    if(lprog) print*
    !
  end subroutine regularlinearinterp2d5v
  !!
  subroutine regularlinearinterp2d5v_homz(x1,y1,v11,v12,v13,v14,v15,   &
                                                     dim1,djm1,        &
                                     x2,y2,v21,v22,v23,v24,v25,        &
                                                           dim2,djm2,km)
    !
    integer,intent(in) :: dim1,djm1,dim2,djm2,km
    real(8),intent(in) :: x1(0:dim1,0:djm1),y1(0:dim1,0:djm1),         &
                          x2(0:dim2,0:djm2),y2(0:dim2,0:djm2)
    real(8),intent(in) :: v11(0:dim1,0:djm1,0:km),v12(0:dim1,0:djm1,0:km),       &
                          v13(0:dim1,0:djm1,0:km),v14(0:dim1,0:djm1,0:km),       &
                          v15(0:dim1,0:djm1,0:km)
    real(8),intent(out) :: v21(0:dim2,0:djm2,0:km),v22(0:dim2,0:djm2,0:km),      &
                           v23(0:dim2,0:djm2,0:km),v24(0:dim2,0:djm2,0:km),      &
                           v25(0:dim2,0:djm2,0:km)
    integer :: i,j,i1,j1,nper
    !
    real(8) :: vti1(0:dim2,0:djm1,0:km),vti2(0:dim2,0:djm1,0:km),                &
               vti3(0:dim2,0:djm1,0:km),vti4(0:dim2,0:djm1,0:km),                &
               vti5(0:dim2,0:djm1,0:km)
    !
    nper=0
    do i=0,dim2
      !
      do i1=1,dim1
        !
        if( x2(i,0)>=x1(i1-1,0) .and. x2(i,0)<=x1(i1,0) ) then
          !
          do j1=0,djm1
            vti1(i,j1,:)=linear1d(x1(i1-1,0), x1(i1,0),                  &
                               v11(i1-1,j1,:),v11(i1,j1,:),x2(i,0))
            vti2(i,j1,:)=linear1d(x1(i1-1,0), x1(i1,0),                  &
                               v12(i1-1,j1,:),v12(i1,j1,:),x2(i,0))
            vti3(i,j1,:)=linear1d(x1(i1-1,0), x1(i1,0),                  &
                               v13(i1-1,j1,:),v13(i1,j1,:),x2(i,0))
            vti4(i,j1,:)=linear1d(x1(i1-1,0), x1(i1,0),                  &
                               v14(i1-1,j1,:),v14(i1,j1,:),x2(i,0))
            vti5(i,j1,:)=linear1d(x1(i1-1,0), x1(i1,0),                  &
                               v15(i1-1,j1,:),v15(i1,j1,:),x2(i,0))
          end do
          !
          exit
          !
        elseif(x2(i,0)>=x1(dim1,0)) then
          !
          ! do j1=0,djm1
          !   vti1(i,j1,:)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
          !                      v11(dim1-1,j1,:),v11(dim1,j1,:),x2(i,0))
          !   vti2(i,j1,:)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
          !                      v12(dim1-1,j1,:),v12(dim1,j1,:),x2(i,0))
          !   vti3(i,j1,:)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
          !                      v13(dim1-1,j1,:),v13(dim1,j1,:),x2(i,0))
          !   vti4(i,j1,:)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
          !                      v14(dim1-1,j1,:),v14(dim1,j1,:),x2(i,0))
          !   vti5(i,j1,:)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
          !                      v15(dim1-1,j1,:),v15(dim1,j1,:),x2(i,0))
          ! end do
          do j1=0,djm1
            vti1(i,j1,:)=v11(dim1-1,j1,:)
            vti2(i,j1,:)=v12(dim1-1,j1,:)
            vti3(i,j1,:)=v13(dim1-1,j1,:)
            vti4(i,j1,:)=v14(dim1-1,j1,:)
            vti5(i,j1,:)=v15(dim1-1,j1,:)
          end do
          !
          exit
          !
        end if
        !
      end do
      !
      write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
                                '  ** Interpolating ... ',100*i/dim2,' %'
    end do
    !
    do j=0,djm2
      !
      do j1=1,djm1
        !
        if( y2(0,j)>=y1(0,j1-1) .and. y2(0,j)<=y1(0,j1) ) then
          do i=0,dim2
            v21(i,j,:)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
                             vti1(i,j1-1,:),vti1(i,j1,:),y2(0,j))
            v22(i,j,:)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
                             vti2(i,j1-1,:),vti2(i,j1,:),y2(0,j))
            v23(i,j,:)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
                             vti3(i,j1-1,:),vti3(i,j1,:),y2(0,j))
            v24(i,j,:)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
                             vti4(i,j1-1,:),vti4(i,j1,:),y2(0,j))
            v25(i,j,:)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
                             vti5(i,j1-1,:),vti5(i,j1,:),y2(0,j))
          end do
          !
          exit
          !
        elseif( y2(0,j)>=y1(0,djm1))  then
          !
          do i=0,dim2
            v21(i,j,:)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                             vti1(i,djm1-1,:),vti1(i,djm1,:),y2(0,j))
            v22(i,j,:)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                             vti2(i,djm1-1,:),vti2(i,djm1,:),y2(0,j))
            v23(i,j,:)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                             vti3(i,djm1-1,:),vti3(i,djm1,:),y2(0,j))
            v24(i,j,:)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                             vti4(i,djm1-1,:),vti4(i,djm1,:),y2(0,j))
            v25(i,j,:)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                             vti5(i,djm1-1,:),vti5(i,djm1,:),y2(0,j))
          end do
          !
          exit
          !
        end if
        !
      end do
      !
      nper=nper+1
      !
      write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
                                '  ** Interpolating ... ',200*j/djm2,' %'
    end do
    ! print*
    !
  end subroutine regularlinearinterp2d5v_homz
  !
  subroutine regularlinearinterp2d6v(x1,y1,v11,v12,v13,v14,v15,v16,    &
                                                     dim1,djm1,        &
                                     x2,y2,v21,v22,v23,v24,v25,v26,    &
                                                              dim2,djm2)
    !
    integer,intent(in) :: dim1,djm1,dim2,djm2
    real(8),intent(in) :: x1(0:dim1,0:djm1),y1(0:dim1,0:djm1),         &
                          v11(0:dim1,0:djm1),v12(0:dim1,0:djm1),       &
                          v13(0:dim1,0:djm1),v14(0:dim1,0:djm1),       &
                          v15(0:dim1,0:djm1),v16(0:dim1,0:djm1)
    real(8),intent(in) :: x2(0:dim2,0:djm2),y2(0:dim2,0:djm2)
    real(8),intent(out) :: v21(0:dim2,0:djm2),v22(0:dim2,0:djm2),      &
                           v23(0:dim2,0:djm2),v24(0:dim2,0:djm2),      &
                           v25(0:dim2,0:djm2),v26(0:dim2,0:djm2)
    integer :: i,j,i1,j1,nper
    !
    real(8) :: vti1(0:dim2,0:djm1),vti2(0:dim2,0:djm1),                &
               vti3(0:dim2,0:djm1),vti4(0:dim2,0:djm1),                &
               vti5(0:dim2,0:djm1),vti6(0:dim2,0:djm1)
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i,0)>=x1(i1-1,0) .and. x2(i,0)<=x1(i1,0) ) then
    		  !
    		  do j1=0,djm1
    		    vti1(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v11(i1-1,j1),v11(i1,j1),x2(i,0))
    		    vti2(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v12(i1-1,j1),v12(i1,j1),x2(i,0))
    		    vti3(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v13(i1-1,j1),v13(i1,j1),x2(i,0))
    		    vti4(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v14(i1-1,j1),v14(i1,j1),x2(i,0))
    		    vti5(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v15(i1-1,j1),v15(i1,j1),x2(i,0))
    		    vti6(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v16(i1-1,j1),v16(i1,j1),x2(i,0))
    		  end do
    		  !
    		  exit
    		  !
    		elseif(x2(i,0)>=x1(dim1,0)) then
              !
              do j1=0,djm1
                vti1(i,j1)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v11(dim1-1,j1),v11(dim1,j1),x2(i,0))
                vti2(i,j1)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v12(dim1-1,j1),v12(dim1,j1),x2(i,0))
                vti3(i,j1)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v13(dim1-1,j1),v13(dim1,j1),x2(i,0))
                vti4(i,j1)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v14(dim1-1,j1),v14(dim1,j1),x2(i,0))
                vti5(i,j1)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v15(dim1-1,j1),v15(dim1,j1),x2(i,0))
                vti6(i,j1)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v16(dim1-1,j1),v16(dim1,j1),x2(i,0))
              end do
              !
              exit
              !
            end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	               '  ** Interpolating ... ',100*nper/(dim2+djm2+2),' %'
    end do
    !
    do j=0,djm2
    	!
    	do j1=1,djm1
    		!
    		if( y2(0,j)>=y1(0,j1-1) .and. y2(0,j)<=y1(0,j1) ) then
    		  do i=0,dim2
    		    v21(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti1(i,j1-1),vti1(i,j1),y2(0,j))
    		    v22(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti2(i,j1-1),vti2(i,j1),y2(0,j))
    		    v23(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti3(i,j1-1),vti3(i,j1),y2(0,j))
    		    v24(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti4(i,j1-1),vti4(i,j1),y2(0,j))
    		    v25(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti5(i,j1-1),vti5(i,j1),y2(0,j))
    		    v26(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti6(i,j1-1),vti6(i,j1),y2(0,j))
    		  end do
    		  !
    		  exit
    		  !
            elseif( y2(0,j)>=y1(0,djm1))  then
              !
              do i=0,dim2
                v21(i,j)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                                 vti1(i,djm1-1),vti1(i,djm1),y2(0,j))
                v22(i,j)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                                 vti2(i,djm1-1),vti2(i,djm1),y2(0,j))
                v23(i,j)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                                 vti3(i,djm1-1),vti3(i,djm1),y2(0,j))
                v24(i,j)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                                 vti4(i,djm1-1),vti4(i,djm1),y2(0,j))
                v25(i,j)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                                 vti5(i,djm1-1),vti5(i,djm1),y2(0,j))
                v26(i,j)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                                 vti6(i,djm1-1),vti6(i,djm1),y2(0,j))
              end do
              !
              exit
              !
            end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	               '  ** Interpolating ... ',200*nper/(dim2+djm2+2),' %'
    end do
    print*
    !
  end subroutine regularlinearinterp2d6v
  !
  subroutine regularlinearinterp2d6v_homz(x1,y1,v11,v12,v13,v14,v15,v16,    &
                                                     dim1,djm1,        &
                                     x2,y2,v21,v22,v23,v24,v25,v26,    &
                                                           dim2,djm2,km)
    !
    integer,intent(in) :: dim1,djm1,dim2,djm2,km
    real(8),intent(in) :: x1(0:dim1,0:djm1),y1(0:dim1,0:djm1),         &
                          x2(0:dim2,0:djm2),y2(0:dim2,0:djm2)
    real(8),intent(in) :: v11(0:dim1,0:djm1,0:km),v12(0:dim1,0:djm1,0:km),       &
                          v13(0:dim1,0:djm1,0:km),v14(0:dim1,0:djm1,0:km),       &
                          v15(0:dim1,0:djm1,0:km),v16(0:dim1,0:djm1,0:km)
    real(8),intent(out) :: v21(0:dim2,0:djm2,0:km),v22(0:dim2,0:djm2,0:km),      &
                           v23(0:dim2,0:djm2,0:km),v24(0:dim2,0:djm2,0:km),      &
                           v25(0:dim2,0:djm2,0:km),v26(0:dim2,0:djm2,0:km)
    integer :: i,j,i1,j1,nper
    !
    real(8) :: vti1(0:dim2,0:djm1,0:km),vti2(0:dim2,0:djm1,0:km),                &
               vti3(0:dim2,0:djm1,0:km),vti4(0:dim2,0:djm1,0:km),                &
               vti5(0:dim2,0:djm1,0:km),vti6(0:dim2,0:djm1,0:km)
    !
    nper=0
    do i=0,dim2
      !
      do i1=1,dim1
        !
        if( x2(i,0)>=x1(i1-1,0) .and. x2(i,0)<=x1(i1,0) ) then
          !
          do j1=0,djm1
            vti1(i,j1,:)=linear1d(x1(i1-1,0), x1(i1,0),                  &
                               v11(i1-1,j1,:),v11(i1,j1,:),x2(i,0))
            vti2(i,j1,:)=linear1d(x1(i1-1,0), x1(i1,0),                  &
                               v12(i1-1,j1,:),v12(i1,j1,:),x2(i,0))
            vti3(i,j1,:)=linear1d(x1(i1-1,0), x1(i1,0),                  &
                               v13(i1-1,j1,:),v13(i1,j1,:),x2(i,0))
            vti4(i,j1,:)=linear1d(x1(i1-1,0), x1(i1,0),                  &
                               v14(i1-1,j1,:),v14(i1,j1,:),x2(i,0))
            vti5(i,j1,:)=linear1d(x1(i1-1,0), x1(i1,0),                  &
                               v15(i1-1,j1,:),v15(i1,j1,:),x2(i,0))
            vti6(i,j1,:)=linear1d(x1(i1-1,0), x1(i1,0),                  &
                               v16(i1-1,j1,:),v16(i1,j1,:),x2(i,0))
          end do
          !
          exit
          !
        elseif(x2(i,0)>=x1(dim1,0)) then
              !
              do j1=0,djm1
                vti1(i,j1,:)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v11(dim1-1,j1,:),v11(dim1,j1,:),x2(i,0))
                vti2(i,j1,:)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v12(dim1-1,j1,:),v12(dim1,j1,:),x2(i,0))
                vti3(i,j1,:)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v13(dim1-1,j1,:),v13(dim1,j1,:),x2(i,0))
                vti4(i,j1,:)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v14(dim1-1,j1,:),v14(dim1,j1,:),x2(i,0))
                vti5(i,j1,:)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v15(dim1-1,j1,:),v15(dim1,j1,:),x2(i,0))
                vti6(i,j1,:)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v16(dim1-1,j1,:),v16(dim1,j1,:),x2(i,0))
              end do
              !
              exit
              !
            end if
        !
      end do
      !
      write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
                                '  ** Interpolating ... ',100*i/dim2,' %'
    end do
    !
    do j=0,djm2
      !
      do j1=1,djm1
        !
        if( y2(0,j)>=y1(0,j1-1) .and. y2(0,j)<=y1(0,j1) ) then
          do i=0,dim2
            v21(i,j,:)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
                             vti1(i,j1-1,:),vti1(i,j1,:),y2(0,j))
            v22(i,j,:)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
                             vti2(i,j1-1,:),vti2(i,j1,:),y2(0,j))
            v23(i,j,:)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
                             vti3(i,j1-1,:),vti3(i,j1,:),y2(0,j))
            v24(i,j,:)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
                             vti4(i,j1-1,:),vti4(i,j1,:),y2(0,j))
            v25(i,j,:)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
                             vti5(i,j1-1,:),vti5(i,j1,:),y2(0,j))
            v26(i,j,:)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
                             vti6(i,j1-1,:),vti6(i,j1,:),y2(0,j))
          end do
          !
          exit
          !
        elseif( y2(0,j)>=y1(0,djm1))  then
          !
          do i=0,dim2
            v21(i,j,:)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                             vti1(i,djm1-1,:),vti1(i,djm1,:),y2(0,j))
            v22(i,j,:)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                             vti2(i,djm1-1,:),vti2(i,djm1,:),y2(0,j))
            v23(i,j,:)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                             vti3(i,djm1-1,:),vti3(i,djm1,:),y2(0,j))
            v24(i,j,:)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                             vti4(i,djm1-1,:),vti4(i,djm1,:),y2(0,j))
            v25(i,j,:)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                             vti5(i,djm1-1,:),vti5(i,djm1,:),y2(0,j))
            v26(i,j,:)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                             vti6(i,djm1-1,:),vti6(i,djm1,:),y2(0,j))
          end do
          !
          exit
          !
        end if
        !
      end do
      !
      nper=nper+1
      !
      write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
                                '  ** Interpolating ... ',200*j/djm2,' %'
    end do
    ! print*
    !
  end subroutine regularlinearinterp2d6v_homz
  !
  subroutine regularlinearinterp2d7v(x1,y1,v11,v12,v13,v14,v15,v16,v17,&
                                                     dim1,djm1,        &
                                     x2,y2,v21,v22,v23,v24,v25,v26,v27,&
                                                              dim2,djm2)
    !
    integer,intent(in) :: dim1,djm1,dim2,djm2
    real(8),intent(in) :: x1(0:dim1,0:djm1),y1(0:dim1,0:djm1),         &
                          v11(0:dim1,0:djm1),v12(0:dim1,0:djm1),       &
                          v13(0:dim1,0:djm1),v14(0:dim1,0:djm1),       &
                          v15(0:dim1,0:djm1),v16(0:dim1,0:djm1),       &
                          v17(0:dim1,0:djm1)
    real(8),intent(in) :: x2(0:dim2,0:djm2),y2(0:dim2,0:djm2)
    real(8),intent(out) :: v21(0:dim2,0:djm2),v22(0:dim2,0:djm2),      &
                           v23(0:dim2,0:djm2),v24(0:dim2,0:djm2),      &
                           v25(0:dim2,0:djm2),v26(0:dim2,0:djm2),      &
                           v27(0:dim2,0:djm2)
    integer :: i,j,i1,j1,nper
    !
    real(8) :: vti1(0:dim2,0:djm1),vti2(0:dim2,0:djm1),                &
               vti3(0:dim2,0:djm1),vti4(0:dim2,0:djm1),                &
               vti5(0:dim2,0:djm1),vti6(0:dim2,0:djm1),                &
               vti7(0:dim2,0:djm1)
    !
    nper=0
    do i=0,dim2
    	!
    	do i1=1,dim1
    		!
    		if( x2(i,0)>=x1(i1-1,0) .and. x2(i,0)<=x1(i1,0) ) then
    		  !
    		  do j1=0,djm1
    		    vti1(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v11(i1-1,j1),v11(i1,j1),x2(i,0))
    		    vti2(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v12(i1-1,j1),v12(i1,j1),x2(i,0))
    		    vti3(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v13(i1-1,j1),v13(i1,j1),x2(i,0))
    		    vti4(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v14(i1-1,j1),v14(i1,j1),x2(i,0))
    		    vti5(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v15(i1-1,j1),v15(i1,j1),x2(i,0))
    		    vti6(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v16(i1-1,j1),v16(i1,j1),x2(i,0))
    		    vti7(i,j1)=linear1d(x1(i1-1,0), x1(i1,0),                  &
    		                       v17(i1-1,j1),v17(i1,j1),x2(i,0))
    		  end do
    		  !
    		  exit
    		  !
            elseif(x2(i,0)>=x1(dim1,0)) then
              !
              do j1=0,djm1
                vti1(i,j1)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v11(dim1-1,j1),v11(dim1,j1),x2(i,0))
                vti2(i,j1)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v12(dim1-1,j1),v12(dim1,j1),x2(i,0))
                vti3(i,j1)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v13(dim1-1,j1),v13(dim1,j1),x2(i,0))
                vti4(i,j1)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v14(dim1-1,j1),v14(dim1,j1),x2(i,0))
                vti5(i,j1)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v15(dim1-1,j1),v15(dim1,j1),x2(i,0))
                vti6(i,j1)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v16(dim1-1,j1),v16(dim1,j1),x2(i,0))
                vti7(i,j1)=linear1d(x1(dim1-1,0), x1(dim1,0),                  &
                                   v17(dim1-1,j1),v17(dim1,j1),x2(i,0))
              end do
              !
              exit
              !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	               '  ** Interpolating ... ',100*nper/(dim2+djm2+2),' %'
    end do
    !
    do j=0,djm2
    	!
    	do j1=1,djm1
    		!
    		if( y2(0,j)>=y1(0,j1-1) .and. y2(0,j)<=y1(0,j1) ) then
    		  do i=0,dim2
    		    v21(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti1(i,j1-1),vti1(i,j1),y2(0,j))
    		    v22(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti2(i,j1-1),vti2(i,j1),y2(0,j))
    		    v23(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti3(i,j1-1),vti3(i,j1),y2(0,j))
    		    v24(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti4(i,j1-1),vti4(i,j1),y2(0,j))
    		    v25(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti5(i,j1-1),vti5(i,j1),y2(0,j))
    		    v26(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti6(i,j1-1),vti6(i,j1),y2(0,j))
    		    v27(i,j)=linear1d(y1(0,j1-1),  y1(0,j1),                   &
    		                     vti7(i,j1-1),vti7(i,j1),y2(0,j))
    		  end do
    		  !
    		  exit
    		  !
            elseif( y2(0,j)>=y1(0,djm1))  then
              !
              do i=0,dim2
                v21(i,j)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                                 vti1(i,djm1-1),vti1(i,djm1),y2(0,j))
                v22(i,j)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                                 vti2(i,djm1-1),vti2(i,djm1),y2(0,j))
                v23(i,j)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                                 vti3(i,djm1-1),vti3(i,djm1),y2(0,j))
                v24(i,j)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                                 vti4(i,djm1-1),vti4(i,djm1),y2(0,j))
                v25(i,j)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                                 vti5(i,djm1-1),vti5(i,djm1),y2(0,j))
                v26(i,j)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                                 vti6(i,djm1-1),vti6(i,djm1),y2(0,j))
                v27(i,j)=linear1d(y1(0,djm1-1),  y1(0,djm1),                   &
                                 vti7(i,djm1-1),vti7(i,djm1),y2(0,j))
              end do
              !
              exit
              !
    		end if
    		!
    	end do
    	!
    	nper=nper+1
    	!
    	write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
    	               '  ** Interpolating ... ',100*nper/(dim2+djm2+2),' %'
    end do
    print*
    !
  end subroutine regularlinearinterp2d7v
  !+-------------------------------------------------------------------+
  !| The end of the subroutine regularlinearinterp2d                   |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is a linear interpolation function.                 |
  !+-------------------------------------------------------------------+
  function linear1d_s(xx1,xx2,yy1,yy2,xx) result(yy)
    !
    real(8),intent(in) :: xx1,xx2,yy1,yy2,xx
    real(8) :: yy
    !
    real(8) :: var1
    !
    var1=(yy2-yy1)/(xx2-xx1)
    yy=var1*(xx-xx1)+yy1
    !
    return
    !
  end function linear1d_s
  !
  function linear1d_a1(xx1,xx2,yy1,yy2,xx) result(yy)
    !
    real(8),intent(in) :: xx1,xx2,xx
    real(8),intent(in) ::  yy1(:),yy2(:)
    real(8) :: yy(1:size(yy1))
    !
    real(8) :: var1
    !
    var1=(xx-xx1)/(xx2-xx1)
    yy=(yy2-yy1)*var1+yy1
    !
    return
    !
  end function linear1d_a1
  !
  function linear1d_a2(xx1,xx2,yy1,yy2,xx) result(yy)
    !
    real(8),intent(in) :: xx1,xx2,xx
    real(8),intent(in) ::  yy1(:,:),yy2(:,:)
    real(8) :: yy(1:size(yy1,1),1:size(yy1,2))
    !
    real(8) :: var1
    !
    var1=(xx-xx1)/(xx2-xx1)
    yy=(yy2-yy1)*var1+yy1
    !
    return
    !
  end function linear1d_a2
  !+-------------------------------------------------------------------+
  !| The end of the function linear1d                                  |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to interpolate flow field between two grids.
  ! The Inverse-Distance Algorithm is used for interpolation in the 
  ! x-y plane
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2010-09-21.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pointinterp(x1,y1,va1im,va2im,dim1,dim2,              &
                         x2,y2,va1in,va2in              )
    !
    integer,intent(in) :: dim1,dim2
    real(8),intent(in) :: x1(0:dim1,0:dim2),y1(0:dim1,0:dim2),         &
                          x2,y2
    real(8),intent(in) :: va1im(0:dim1,0:dim2),va2im(0:dim1,0:dim2)
    real(8),intent(out) :: va1in,va2in
    !
    integer :: i,j,k1,is,js,n1,n2
    real(8) :: var1,dismin
    integer :: perc
    real(8) :: wei(-1:1,-1:1),dis(-1:1,-1:1),temp(2)
    !
    !!$ SAVE wei,dis
    !
    !!$OMP THREADPRIVATE(wei,dis)
    !
    !!$OMP parallel default(shared) private(i,j,i1,j1,k1,is,js,n1,n2,   &
    !!$OMP                                              var1,dismin,perc)
    !!$OMP do
      dismin=1d10
      do j=0,dim2
      do i=0,dim1
        !
        var1=(x1(i,j)-x2)**2+(y1(i,j)-y2)**2
        if(var1<dismin) then
          dismin=var1
          is=i
          js=j
        end if
        !
      end do
      end do
      !
      if(dismin<=1.d-15) then
        !
        va1in=va1im(is,js)
        va2in=va2im(is,js)
        !
      else
        !
        var1=0.d0
        do n1=-1,1
        do n2=-1,1
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            wei(n1,n2)=0.d0
          else
            dis(n1,n2)=(x2-x1(is+n1,js+n2))**2+                 &
                       (y2-y1(is+n1,js+n2))**2
            wei(n1,n2)=1.d0/dis(n1,n2)
          end if
          !
          var1=var1+wei(n1,n2)
        end do
        end do
        !
        do n1=-1,1
        do n2=-1,1
          wei(n1,n2)=wei(n1,n2)/var1
        end do
        end do
        !
        temp=0.d0
        do n1=-1,1
        do n2=-1,1
          !
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            continue
          else
            temp(1)=temp(1)+wei(n1,n2)*va1im(is+n1,js+n2)
            temp(2)=temp(2)+wei(n1,n2)*va2im(is+n1,js+n2)
          end if
          !
        end do
        end do
        !
        va1in=temp(1)
        va2in=temp(2)
        !
      end if
      !
    !!$OMP end do
    !!$OMP end parallel
    !
  end subroutine pointinterp
  !
  subroutine InvdisInterp2D1(x1,y1,va1im,dim1,dim2,                    &
                             x2,y2,va1in,din1,din2                    )
    !
    integer,intent(in) :: dim1,dim2,din1,din2
    real(8),intent(in) :: x1(0:dim1,0:dim2),y1(0:dim1,0:dim2),         &
                          x2(0:din1,0:din2),y2(0:din1,0:din2)
    real(8),intent(in) :: va1im(0:dim1,0:dim2)
    real(8),intent(inout) :: va1in(0:din1,0:din2)
    !
    integer :: i,j,i1,j1,k1,is,js,n1,n2
    real(8) :: var1,dismin
    integer :: perc
    real(8) :: wei(-1:1,-1:1),dis(-1:1,-1:1),temp(1)
    !
    !!$ SAVE wei,dis
    !
    !!$OMP THREADPRIVATE(wei,dis)
    !
    !!$OMP parallel default(shared) private(i,j,i1,j1,k1,is,js,n1,n2,   &
    !!$OMP                                              var1,dismin,perc)
    !!$OMP do
    do j1=0,din2
    do i1=0,din1
      !
      dismin=1d10
      do j=0,dim2
      do i=0,dim1
        !
        var1=(x1(i,j)-x2(i1,j1))**2+(y1(i,j)-y2(i1,j1))**2
        if(var1<dismin) then
          dismin=var1
          is=i
          js=j
        end if
        !
      end do
      end do
      !
      if(dismin<=1.d-15) then
        !
        va1in(i1,j1)=va1im(is,js)
        !
      else
        !
        var1=0.d0
        do n1=-1,1
        do n2=-1,1
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            wei(n1,n2)=0.d0
          else
            dis(n1,n2)=(x2(i1,j1)-x1(is+n1,js+n2))**2+                 &
                       (y2(i1,j1)-y1(is+n1,js+n2))**2
            wei(n1,n2)=1.d0/dis(n1,n2)
          end if
          !
          var1=var1+wei(n1,n2)
        end do
        end do
        !
        do n1=-1,1
        do n2=-1,1
          wei(n1,n2)=wei(n1,n2)/var1
        end do
        end do
        !
        temp=0.d0
        do n1=-1,1
        do n2=-1,1
          !
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            continue
          else
            temp(1)=temp(1)+wei(n1,n2)*va1im(is+n1,js+n2)
          end if
          !
        end do
        end do
        !
        va1in(i1,j1)=temp(1)
        !
      end if
      !
    end do
      !
      perc=100*(j1)/(din2)
      write(*,'(1A1,A23,I3,A2,$)')char(13),'  ** Interpolating ... ',perc,' %'
      !
    end do
    print*
    !!$OMP end do
    !!$OMP end parallel
    !
  end subroutine InvdisInterp2D1
  !
  subroutine InvdisInterp2D2(x1,y1,va1im,va2im,dim1,dim2,              &
                             x2,y2,va1in,va2in,din1,din2              )
    !
    integer,intent(in) :: dim1,dim2,din1,din2
    real(8),intent(in) :: x1(0:dim1,0:dim2),y1(0:dim1,0:dim2),         &
                          x2(0:din1,0:din2),y2(0:din1,0:din2)
    real(8),intent(in) :: va1im(0:dim1,0:dim2),va2im(0:dim1,0:dim2)
    real(8),intent(inout) :: va1in(0:din1,0:din2),va2in(0:din1,0:din2)
    !
    integer :: i,j,i1,j1,k1,is,js,n1,n2
    real(8) :: var1,dismin
    integer :: perc
    real(8) :: wei(-1:1,-1:1),dis(-1:1,-1:1),temp(2)
    !
    !!$ SAVE wei,dis
    !
    !!$OMP THREADPRIVATE(wei,dis)
    !
    !!$OMP parallel default(shared) private(i,j,i1,j1,k1,is,js,n1,n2,   &
    !!$OMP                                              var1,dismin,perc)
    !!$OMP do
    do j1=0,din2
    do i1=0,din1
      !
      dismin=1d10
      do j=0,dim2
      do i=0,dim1
        !
        var1=(x1(i,j)-x2(i1,j1))**2+(y1(i,j)-y2(i1,j1))**2
        if(var1<dismin) then
          dismin=var1
          is=i
          js=j
        end if
        !
      end do
      end do
      !
      if(dismin<=1.d-15) then
        !
        va1in(i1,j1)=va1im(is,js)
        va2in(i1,j1)=va2im(is,js)
        !
      else
        !
        var1=0.d0
        do n1=-1,1
        do n2=-1,1
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            wei(n1,n2)=0.d0
          else
            dis(n1,n2)=(x2(i1,j1)-x1(is+n1,js+n2))**2+                 &
                       (y2(i1,j1)-y1(is+n1,js+n2))**2
            wei(n1,n2)=1.d0/dis(n1,n2)
          end if
          !
          var1=var1+wei(n1,n2)
        end do
        end do
        !
        do n1=-1,1
        do n2=-1,1
          wei(n1,n2)=wei(n1,n2)/var1
        end do
        end do
        !
        temp=0.d0
        do n1=-1,1
        do n2=-1,1
          !
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            continue
          else
            temp(1)=temp(1)+wei(n1,n2)*va1im(is+n1,js+n2)
            temp(2)=temp(2)+wei(n1,n2)*va2im(is+n1,js+n2)
          end if
          !
        end do
        end do
        !
        va1in(i1,j1)=temp(1)
        va2in(i1,j1)=temp(2)
        !
      end if
      !
    end do
      !
      perc=100*(j1)/(din2)
      write(*,'(1A1,A23,I3,A2,$)')char(13),'  ** Interpolating ... ',perc,' %'
      !
    end do
    !!$OMP end do
    !!$OMP end parallel
    !
  end subroutine InvdisInterp2D2
  !!
  subroutine InvdisInterp2D3(x1,y1,va1im,va2im,va3im,dim1,dim2,        &
                             x2,y2,va1in,va2in,va3in,din1,din2        )
    !
    integer,intent(in) :: dim1,dim2,din1,din2
    real(8),intent(in) :: x1(0:dim1,0:dim2),y1(0:dim1,0:dim2),         &
                          x2(0:din1,0:din2),y2(0:din1,0:din2)
    real(8),intent(in) :: va1im(0:dim1,0:dim2),va2im(0:dim1,0:dim2),   &
                          va3im(0:dim1,0:dim2)
    real(8),intent(inout) :: va1in(0:din1,0:din2),va2in(0:din1,0:din2),&
                             va3in(0:din1,0:din2)
    !
    integer :: i,j,i1,j1,k1,is,js,n1,n2
    real(8) :: var1,dismin
    integer :: perc
    real(8) :: wei(-1:1,-1:1),dis(-1:1,-1:1),temp(3)
    !
    !!$ SAVE wei,dis
    !
    !!$OMP THREADPRIVATE(wei,dis)
    !
    !!$OMP parallel default(shared) private(i,j,i1,j1,k1,is,js,n1,n2,   &
    !!$OMP                                              var1,dismin,perc)
    !!$OMP do
    do j1=0,din2
    do i1=0,din1
      !
      dismin=1d10
      do j=0,dim2
      do i=0,dim1
        !
        var1=(x1(i,j)-x2(i1,j1))**2+(y1(i,j)-y2(i1,j1))**2
        if(var1<dismin) then
          dismin=var1
          is=i
          js=j
        end if
        !
      end do
      end do
      !
      if(dismin<=1.d-15) then
        !
        va1in(i1,j1)=va1im(is,js)
        va2in(i1,j1)=va2im(is,js)
        va3in(i1,j1)=va3im(is,js)
        !
      else
        !
        var1=0.d0
        do n1=-1,1
        do n2=-1,1
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            wei(n1,n2)=0.d0
          else
            dis(n1,n2)=(x2(i1,j1)-x1(is+n1,js+n2))**2+                 &
                       (y2(i1,j1)-y1(is+n1,js+n2))**2
            wei(n1,n2)=1.d0/dis(n1,n2)
          end if
          !
          var1=var1+wei(n1,n2)
        end do
        end do
        !
        do n1=-1,1
        do n2=-1,1
          wei(n1,n2)=wei(n1,n2)/var1
        end do
        end do
        !
        temp=0.d0
        do n1=-1,1
        do n2=-1,1
          !
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            continue
          else
            temp(1)=temp(1)+wei(n1,n2)*va1im(is+n1,js+n2)
            temp(2)=temp(2)+wei(n1,n2)*va2im(is+n1,js+n2)
            temp(3)=temp(3)+wei(n1,n2)*va3im(is+n1,js+n2)
          end if
          !
        end do
        end do
        !
        va1in(i1,j1)=temp(1)
        va2in(i1,j1)=temp(2)
        va3in(i1,j1)=temp(3)
        !
      end if
      !
    end do
      !
      perc=100*(j1)/(din2)
      write(*,'(1A1,A23,I3,A2,$)')char(13),'  ** Interpolating ... ',perc,' %'
      !
    end do
    !!$OMP end do
    !!$OMP end parallel
    !
  end subroutine InvdisInterp2D3
  !
  subroutine InvdisInterp2D4(x1,y1,va1im,va2im,va3im,va4im,dim1,dim2,  &
                             x2,y2,va1in,va2in,va3in,va4in,din1,din2  )
    !
    integer,intent(in) :: dim1,dim2,din1,din2
    real(8),intent(in) :: x1(0:dim1,0:dim2),y1(0:dim1,0:dim2),         &
                          x2(0:din1,0:din2),y2(0:din1,0:din2)
    real(8),intent(in) :: va1im(0:dim1,0:dim2),va2im(0:dim1,0:dim2),   &
                          va3im(0:dim1,0:dim2),va4im(0:dim1,0:dim2)
    real(8),intent(inout) :: va1in(0:din1,0:din2),va2in(0:din1,0:din2),&
                             va3in(0:din1,0:din2),va4in(0:din1,0:din2)
    !
    integer :: i,j,i1,j1,k1,is,js,n1,n2
    real(8) :: var1,dismin
    integer :: perc
    real(8) :: wei(-1:1,-1:1),dis(-1:1,-1:1),temp(4)
    !
    !!$ SAVE wei,dis
    !
    !!$OMP THREADPRIVATE(wei,dis)
    !
    !!$OMP parallel default(shared) private(i,j,i1,j1,k1,is,js,n1,n2,   &
    !!$OMP                                              var1,dismin,perc)
    !!$OMP do
    open(18,file='interprofile.dat')
    do j1=0,din2
    do i1=0,din1
      !
      dismin=1d10
      do j=0,dim2
      do i=0,dim1
        !
        var1=(x1(i,j)-x2(i1,j1))**2+(y1(i,j)-y2(i1,j1))**2
        if(var1<dismin) then
          dismin=var1
          is=i
          js=j
        end if
        !
      end do
      end do
      !
      if(dismin<=1.d-15) then
        !
        va1in(i1,j1)=va1im(is,js)
        va2in(i1,j1)=va2im(is,js)
        va3in(i1,j1)=va3im(is,js)
        va4in(i1,j1)=va4im(is,js)
        !
      else
        !
        var1=0.d0
        do n1=-1,1
        do n2=-1,1
        	!
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            wei(n1,n2)=0.d0
          else
            dis(n1,n2)=(x2(i1,j1)-x1(is+n1,js+n2))**2+                 &
                       (y2(i1,j1)-y1(is+n1,js+n2))**2
            wei(n1,n2)=1.d0/dis(n1,n2)
          end if
          !
          var1=var1+wei(n1,n2)
        end do
        end do
        !
        do n1=-1,1
        do n2=-1,1
          wei(n1,n2)=wei(n1,n2)/var1
        end do
        end do
        !
        temp=0.d0
        do n1=-1,1
        do n2=-1,1
          !
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            continue
          else
            temp(1)=temp(1)+wei(n1,n2)*va1im(is+n1,js+n2)
            temp(2)=temp(2)+wei(n1,n2)*va2im(is+n1,js+n2)
            temp(3)=temp(3)+wei(n1,n2)*va3im(is+n1,js+n2)
            temp(4)=temp(4)+wei(n1,n2)*va4im(is+n1,js+n2)
          end if
          !
        end do
        end do
        !
        va1in(i1,j1)=temp(1)
        va2in(i1,j1)=temp(2)
        va3in(i1,j1)=temp(3)
        va4in(i1,j1)=temp(4)
        !
      end if
      !
      if(i1==0) then
        write(18,*)x2(i1,j1),y2(i1,j1),va1in(i1,j1)
      end if
      !
    end do
      !
      if( mod(j,10)==0 ) then
        perc=100*(j1)/(din2)
        write(*,'(1A1,A23,I3,A2,$)')char(13),'  ** Interpolating ... ',perc,' %'
      end if
      !
    end do
    print*,'  ** Interpolating ------- done.'
    !
    close(18)
    !
    !!$OMP end do
    !!$OMP end parallel
    !
  end subroutine InvdisInterp2D4
  !
  subroutine InvdisInterp2D5(x1,y1,va1im,va2im,va3im,va4im,va5im,      &
                             dim1,dim2,  &
                             x2,y2,va1in,va2in,va3in,va4in,va5in,      &
                             din1,din2   )
    !
    integer,intent(in) :: dim1,dim2,din1,din2
    real(8),intent(in) :: x1(0:dim1,0:dim2),y1(0:dim1,0:dim2),         &
                          x2(0:din1,0:din2),y2(0:din1,0:din2)
    real(8),intent(in) :: va1im(0:dim1,0:dim2),va2im(0:dim1,0:dim2),   &
                          va3im(0:dim1,0:dim2),va4im(0:dim1,0:dim2),   &
                          va5im(0:dim1,0:dim2)
    real(8),intent(inout) :: va1in(0:din1,0:din2),va2in(0:din1,0:din2),&
                             va3in(0:din1,0:din2),va4in(0:din1,0:din2),&
                             va5in(0:din1,0:din2)
    !
    integer :: i,j,i1,j1,k1,is,js,n1,n2
    real(8) :: var1,dismin
    integer :: perc
    real(8) :: wei(-1:1,-1:1),dis(-1:1,-1:1),temp(5)
    !
    !$ SAVE wei,dis,temp
    
    !$OMP THREADPRIVATE(wei,dis,temp)
    !
    !$OMP parallel default(shared) private(i,j,i1,j1,k1,is,js,n1,n2,   &
    !$OMP                                              var1,dismin,perc)
    !$OMP do
    do j1=0,din2
    do i1=0,din1
      !
      dismin=1d10
      do j=0,dim2
      do i=0,dim1
        !
        var1=(x1(i,j)-x2(i1,j1))**2+(y1(i,j)-y2(i1,j1))**2
        if(var1<dismin) then
          dismin=var1
          is=i
          js=j
        end if
        !
      end do
      end do
      !
      if(dismin<=1.d-15) then
        !
        va1in(i1,j1)=va1im(is,js)
        va2in(i1,j1)=va2im(is,js)
        va3in(i1,j1)=va3im(is,js)
        va4in(i1,j1)=va4im(is,js)
        va5in(i1,j1)=va5im(is,js)
        !
      else
        !
        var1=0.d0
        do n1=-1,1
        do n2=-1,1
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            wei(n1,n2)=0.d0
          else
            dis(n1,n2)=(x2(i1,j1)-x1(is+n1,js+n2))**2+                 &
                       (y2(i1,j1)-y1(is+n1,js+n2))**2
            wei(n1,n2)=1.d0/dis(n1,n2)
          end if
          !
          var1=var1+wei(n1,n2)
        end do
        end do
        !
        do n1=-1,1
        do n2=-1,1
          wei(n1,n2)=wei(n1,n2)/var1
        end do
        end do
        !
        temp=0.d0
        do n1=-1,1
        do n2=-1,1
          !
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            continue
          else
            temp(1)=temp(1)+wei(n1,n2)*va1im(is+n1,js+n2)
            temp(2)=temp(2)+wei(n1,n2)*va2im(is+n1,js+n2)
            temp(3)=temp(3)+wei(n1,n2)*va3im(is+n1,js+n2)
            temp(4)=temp(4)+wei(n1,n2)*va4im(is+n1,js+n2)
            temp(5)=temp(5)+wei(n1,n2)*va5im(is+n1,js+n2)
          end if
          !
        end do
        end do
        !
        va1in(i1,j1)=temp(1)
        va2in(i1,j1)=temp(2)
        va3in(i1,j1)=temp(3)
        va4in(i1,j1)=temp(4)
        va5in(i1,j1)=temp(5)
        !
      end if
      !
    end do
      !
      perc=100*(j1)/(din2)
      write(*,'(1A1,A23,I3,A2,$)')char(13),'  ** Interpolating ... ',perc,' %'
      !
    end do
    !$OMP end do
    !$OMP end parallel
    !
  end subroutine InvdisInterp2D5
  !
  subroutine InvdisInterp2D6(x1,y1,va1im,va2im,va3im,va4im,va5im,      &
                             va6im,dim1,dim2,  &
                             x2,y2,va1in,va2in,va3in,va4in,va5in,      &
                             va6in,din1,din2   )
    !
    integer,intent(in) :: dim1,dim2,din1,din2
    real(8),intent(in) :: x1(0:dim1,0:dim2),y1(0:dim1,0:dim2),         &
                          x2(0:din1,0:din2),y2(0:din1,0:din2)
    real(8),intent(in) :: va1im(0:dim1,0:dim2),va2im(0:dim1,0:dim2),   &
                          va3im(0:dim1,0:dim2),va4im(0:dim1,0:dim2),   &
                          va5im(0:dim1,0:dim2),va6im(0:dim1,0:dim2)
    real(8),intent(inout) :: va1in(0:din1,0:din2),va2in(0:din1,0:din2),&
                             va3in(0:din1,0:din2),va4in(0:din1,0:din2),&
                             va5in(0:din1,0:din2),va6in(0:din1,0:din2)
    !
    integer :: i,j,i1,j1,k1,is,js,n1,n2
    real(8) :: var1,dismin
    real(4) :: perc
    real(8) :: wei(-1:1,-1:1),dis(-1:1,-1:1),temp(6)
    !
    !$ SAVE wei,dis,temp
    !
    !$OMP THREADPRIVATE(wei,dis,temp)
    !
    !$OMP parallel default(shared) private(i,j,i1,j1,k1,is,js,n1,n2,   &
    !$OMP                                              var1,dismin,perc)
    !$OMP do
    do j1=0,din2
    do i1=0,din1
      !
      dismin=1d10
      do j=0,dim2
      do i=0,dim1
        !
        var1=(x1(i,j)-x2(i1,j1))**2+(y1(i,j)-y2(i1,j1))**2
        if(var1<dismin) then
          dismin=var1
          is=i
          js=j
        end if
        !
      end do
      end do
      !
      if(dismin<=1.d-15) then
        !
        va1in(i1,j1)=va1im(is,js)
        va2in(i1,j1)=va2im(is,js)
        va3in(i1,j1)=va3im(is,js)
        va4in(i1,j1)=va4im(is,js)
        va5in(i1,j1)=va5im(is,js)
        va6in(i1,j1)=va6im(is,js)
        !
      else
        !
        var1=0.d0
        do n1=-1,1
        do n2=-1,1
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            wei(n1,n2)=0.d0
          else
            dis(n1,n2)=(x2(i1,j1)-x1(is+n1,js+n2))**2+                 &
                       (y2(i1,j1)-y1(is+n1,js+n2))**2
            wei(n1,n2)=1.d0/dis(n1,n2)
          end if
          !
          var1=var1+wei(n1,n2)
        end do
        end do
        !
        do n1=-1,1
        do n2=-1,1
          wei(n1,n2)=wei(n1,n2)/var1
        end do
        end do
        !
        temp=0.d0
        do n1=-1,1
        do n2=-1,1
          !
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            continue
          else
            temp(1)=temp(1)+wei(n1,n2)*va1im(is+n1,js+n2)
            temp(2)=temp(2)+wei(n1,n2)*va2im(is+n1,js+n2)
            temp(3)=temp(3)+wei(n1,n2)*va3im(is+n1,js+n2)
            temp(4)=temp(4)+wei(n1,n2)*va4im(is+n1,js+n2)
            temp(5)=temp(5)+wei(n1,n2)*va5im(is+n1,js+n2)
            temp(6)=temp(6)+wei(n1,n2)*va6im(is+n1,js+n2)
          end if
          !
        end do
        end do
        !
        va1in(i1,j1)=temp(1)
        va2in(i1,j1)=temp(2)
        va3in(i1,j1)=temp(3)
        va4in(i1,j1)=temp(4)
        va5in(i1,j1)=temp(5)
        va6in(i1,j1)=temp(6)
        !
      end if
      !
    end do
      !
      perc=100.*(j1)/(din2)
      !
      write(*,'(1A1,A23,F5.1,A2,$)')char(13),'  ** Interpolating ... ',perc,' %'
    end do
    !$OMP end do
    !$OMP end parallel
    !
  end subroutine InvdisInterp2D6
  !
  subroutine InvdisInterp2D6_km(x1,y1,va1im,va2im,va3im,va4im,va5im,va6im,   &
                                x2,y2,va1in,va2in,va3in,va4in,va5in,va6in )
    !
    real(8),intent(in),dimension(0:,0:)  :: x1,y1,x2,y2
    real(8),intent(in),dimension(0:,0:,0:) ::    va1im,va2im,va3im,va4im,va5im,va6im
    real(8),intent(inout),dimension(0:,0:,0:) :: va1in,va2in,va3in,va4in,va5in,va6in
    !
    integer :: dim1,dim2,dim3,din1,din2,din3,dimz
    integer :: i,j,k,i1,j1,k1,is,js,n1,n2
    real(8) :: var1,dismin
    integer :: perc
    real(8) :: wei(-1:1,-1:1),dis(-1:1,-1:1)
    real(8),allocatable :: temp(:,:)
    !
    !$ SAVE wei,dis,temp
    !$OMP THREADPRIVATE(wei,dis,temp)
    !
    dim1=size(va1im,1)-1
    dim2=size(va1im,2)-1
    dim3=size(va1im,3)-1
    !
    din1=size(va1in,1)-1
    din2=size(va1in,2)-1
    din3=size(va1in,3)-1
    !
    if(din3 .ne. dim3) then
      print*,' !! the third dimension not match',dim3,din3
      stop 
    else
      dimz=dim3
    endif
    !
    print*,dim1,dim2,dim3
    print*,din1,din2,din3
    !
    !$OMP parallel default(shared) private(i,j,i1,j1,k1,is,js,n1,n2,   &
    !$OMP                                              var1,dismin,perc)
    allocate(temp(6,0:dimz))
    !
    !$OMP do
    do j1=0,din2
    do i1=0,din1
      !
      dismin=1d10
      do j=0,dim2
      do i=0,dim1
        !
        var1=(x1(i,j)-x2(i1,j1))**2+(y1(i,j)-y2(i1,j1))**2
        if(var1<dismin) then
          dismin=var1
          is=i
          js=j
        end if
        !
      end do
      end do
      !
      if(dismin<=1.d-15) then
        !
        va1in(i1,j1,:)=va1im(is,js,:)
        va2in(i1,j1,:)=va2im(is,js,:)
        va3in(i1,j1,:)=va3im(is,js,:)
        va4in(i1,j1,:)=va4im(is,js,:)
        va5in(i1,j1,:)=va5im(is,js,:)
        va6in(i1,j1,:)=va6im(is,js,:)
        !
      else
        !
        var1=0.d0
        do n1=-1,1
        do n2=-1,1
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            wei(n1,n2)=0.d0
          else
            dis(n1,n2)=(x2(i1,j1)-x1(is+n1,js+n2))**2+                 &
                       (y2(i1,j1)-y1(is+n1,js+n2))**2
            wei(n1,n2)=1.d0/dis(n1,n2)
          end if
          !
          var1=var1+wei(n1,n2)
        end do
        end do
        !
        do n1=-1,1
        do n2=-1,1
          wei(n1,n2)=wei(n1,n2)/var1
        end do
        end do
        !
        temp=0.d0
        do n1=-1,1
        do n2=-1,1
          !
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            continue
          else
            do k=0,dimz
              temp(1,k)=temp(1,k)+wei(n1,n2)*va1im(is+n1,js+n2,k)
              temp(2,k)=temp(2,k)+wei(n1,n2)*va2im(is+n1,js+n2,k)
              temp(3,k)=temp(3,k)+wei(n1,n2)*va3im(is+n1,js+n2,k)
              temp(4,k)=temp(4,k)+wei(n1,n2)*va4im(is+n1,js+n2,k)
              temp(5,k)=temp(5,k)+wei(n1,n2)*va5im(is+n1,js+n2,k)
              temp(6,k)=temp(6,k)+wei(n1,n2)*va6im(is+n1,js+n2,k)
            enddo
          end if
          !
        end do
        end do
        !
        va1in(i1,j1,:)=temp(1,:)
        va2in(i1,j1,:)=temp(2,:)
        va3in(i1,j1,:)=temp(3,:)
        va4in(i1,j1,:)=temp(4,:)
        va5in(i1,j1,:)=temp(5,:)
        va6in(i1,j1,:)=temp(6,:)
        !
      end if
      !
    end do
      !
      perc=100*(j1)/(din2)
      write(*,'(1A1,A23,I3,A2,$)')char(13),'  ** Interpolating ... ',perc,' %'
      !
    end do
    !$OMP end do
    !$OMP end parallel
    !
  end subroutine InvdisInterp2D6_km
  !
  subroutine InvdisInterp2D4_km(x1,y1,va1im,va2im,va3im,va4im,   &
                                x2,y2,va1in,va2in,va3in,va4in)
    !
    real(8),intent(in),dimension(0:,0:)  :: x1,y1,x2,y2
    real(8),intent(in),dimension(0:,0:,0:) ::    va1im,va2im,va3im,va4im
    real(8),intent(inout),dimension(0:,0:,0:) :: va1in,va2in,va3in,va4in
    !
    integer :: dim1,dim2,dim3,din1,din2,din3,dimz
    integer :: i,j,k,i1,j1,k1,is,js,n1,n2
    real(8) :: var1,dismin
    integer :: perc
    real(8) :: wei(-1:1,-1:1),dis(-1:1,-1:1)
    real(8),allocatable :: temp(:,:)
    !
    !$ SAVE wei,dis,temp
    !$OMP THREADPRIVATE(wei,dis,temp)
    !
    dim1=size(va1im,1)-1
    dim2=size(va1im,2)-1
    dim3=size(va1im,3)-1
    !
    din1=size(va1in,1)-1
    din2=size(va1in,2)-1
    din3=size(va1in,3)-1
    !
    if(din3 .ne. dim3) then
      print*,' !! the third dimension not match',dim3,din3
      stop 
    else
      dimz=dim3
    endif
    !
    print*,dim1,dim2,dim3
    print*,din1,din2,din3
    !
    !$OMP parallel default(shared) private(i,j,i1,j1,k1,is,js,n1,n2,   &
    !$OMP                                              var1,dismin,perc)
    allocate(temp(6,0:dimz))
    !
    !$OMP do
    do j1=0,din2
    do i1=0,din1
      !
      dismin=1d10
      do j=0,dim2
      do i=0,dim1
        !
        var1=(x1(i,j)-x2(i1,j1))**2+(y1(i,j)-y2(i1,j1))**2
        if(var1<dismin) then
          dismin=var1
          is=i
          js=j
        end if
        !
      end do
      end do
      !
      if(dismin<=1.d-15) then
        !
        va1in(i1,j1,:)=va1im(is,js,:)
        va2in(i1,j1,:)=va2im(is,js,:)
        va3in(i1,j1,:)=va3im(is,js,:)
        va4in(i1,j1,:)=va4im(is,js,:)
        !
      else
        !
        var1=0.d0
        do n1=-1,1
        do n2=-1,1
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            wei(n1,n2)=0.d0
          else
            dis(n1,n2)=(x2(i1,j1)-x1(is+n1,js+n2))**2+                 &
                       (y2(i1,j1)-y1(is+n1,js+n2))**2
            wei(n1,n2)=1.d0/dis(n1,n2)
          end if
          !
          var1=var1+wei(n1,n2)
        end do
        end do
        !
        do n1=-1,1
        do n2=-1,1
          wei(n1,n2)=wei(n1,n2)/var1
        end do
        end do
        !
        temp=0.d0
        do n1=-1,1
        do n2=-1,1
          !
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            continue
          else
            do k=0,dimz
              temp(1,k)=temp(1,k)+wei(n1,n2)*va1im(is+n1,js+n2,k)
              temp(2,k)=temp(2,k)+wei(n1,n2)*va2im(is+n1,js+n2,k)
              temp(3,k)=temp(3,k)+wei(n1,n2)*va3im(is+n1,js+n2,k)
              temp(4,k)=temp(4,k)+wei(n1,n2)*va4im(is+n1,js+n2,k)
            enddo
          end if
          !
        end do
        end do
        !
        va1in(i1,j1,:)=temp(1,:)
        va2in(i1,j1,:)=temp(2,:)
        va3in(i1,j1,:)=temp(3,:)
        va4in(i1,j1,:)=temp(4,:)
        !
      end if
      !
    end do
      !
      perc=100*(j1)/(din2)
      write(*,'(1A1,A23,I3,A2,$)')char(13),'  ** Interpolating ... ',perc,' %'
      !
    end do
    !$OMP end do
    !$OMP end parallel
    !
  end subroutine InvdisInterp2D4_km
  !
  subroutine InvdisInterp2D4_1_km(x1,y1,va1im,va2im,va3im,va4im,   &
                                  x2,y2,va1in,va2in,va3in,va4in)
    !
    real(8),intent(in),dimension(0:,0:)  :: x1,y1
    real(8),intent(in),dimension(0:)  :: x2,y2
    real(8),intent(in),dimension(0:,0:,0:) ::    va1im,va2im,va3im,va4im
    real(8),intent(inout),dimension(0:,0:) :: va1in,va2in,va3in,va4in
    !
    integer :: dim1,dim2,dim3,din1,din2,din3,dimz
    integer :: i,j,k,i1,j1,k1,is,js,n1,n2
    real(8) :: var1,dismin
    integer :: perc
    real(8) :: wei(-1:1,-1:1),dis(-1:1,-1:1)
    real(8),allocatable :: temp(:,:)
    !
    !!$ SAVE wei,dis,temp
    !!$OMP THREADPRIVATE(wei,dis,temp)
    !
    dim1=size(va1im,1)-1
    dim2=size(va1im,2)-1
    dim3=size(va1im,3)-1
    !
    din1=size(va1in,1)-1
    din3=size(va1in,2)-1
    !
    if(din3 .ne. dim3) then
      print*,' !! the third dimension not match',dim3,din3
      stop 
    else
      dimz=dim3
    endif
    !
    din2=0
    !
    print*,dim1,dim2,dim3
    print*,din1,din2,din3
    !
    !!$OMP parallel default(shared) private(i,j,i1,j1,k1,is,js,n1,n2,   &
    !!$OMP                                              var1,dismin,perc)
    allocate(temp(4,0:dimz))
    !
    !!$OMP do
    do i1=0,din1
      !
      dismin=1d10
      do j=0,dim2
      do i=0,dim1
        !
        var1=(x1(i,j)-x2(i1))**2+(y1(i,j)-y2(i1))**2
        if(var1<dismin) then
          dismin=var1
          is=i
          js=j
        end if
        !
      end do
      end do
      !
      if(dismin<=1.d-15) then
        !
        va1in(i1,:)=va1im(is,js,:)
        va2in(i1,:)=va2im(is,js,:)
        va3in(i1,:)=va3im(is,js,:)
        va4in(i1,:)=va4im(is,js,:)
        !
      else
        !
        var1=0.d0
        do n1=-1,1
        do n2=-1,1
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            wei(n1,n2)=0.d0
          else
            dis(n1,n2)=(x2(i1)-x1(is+n1,js+n2))**2+                 &
                       (y2(i1)-y1(is+n1,js+n2))**2
            wei(n1,n2)=1.d0/dis(n1,n2)
          end if
          !
          var1=var1+wei(n1,n2)
        end do
        end do
        !
        wei=wei/var1
        !
        temp=0.d0
        do n1=-1,1
        do n2=-1,1
          !
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            continue
          else
            do k=0,dimz
              temp(1,k)=temp(1,k)+wei(n1,n2)*va1im(is+n1,js+n2,k)
              temp(2,k)=temp(2,k)+wei(n1,n2)*va2im(is+n1,js+n2,k)
              temp(3,k)=temp(3,k)+wei(n1,n2)*va3im(is+n1,js+n2,k)
              temp(4,k)=temp(4,k)+wei(n1,n2)*va4im(is+n1,js+n2,k)
            enddo
          end if
          !
        end do
        end do
        !
        va1in(i1,:)=temp(1,:)
        va2in(i1,:)=temp(2,:)
        va3in(i1,:)=temp(3,:)
        va4in(i1,:)=temp(4,:)
        !
      end if
      !
      perc=100*(i1)/(din1)
      write(*,'(1A1,A23,I3,A2,$)')char(13),'  ** Interpolating ... ',perc,' %'
      !
    end do
    !!$OMP end do
    !
    deallocate(temp)
    !
    !!$OMP end parallel
    !
  end subroutine InvdisInterp2D4_1_km
  !
  subroutine InvdisInterp2D5_km(x1,y1,va1im,va2im,va3im,va4im,va5im,   &
                                x2,y2,va1in,va2in,va3in,va4in,va5in)
    !
    real(8),intent(in),dimension(0:,0:)  :: x1,y1,x2,y2
    real(8),intent(in),dimension(0:,0:,0:) ::    va1im,va2im,va3im,va4im,va5im
    real(8),intent(inout),dimension(0:,0:,0:) :: va1in,va2in,va3in,va4in,va5in
    !
    integer :: dim1,dim2,dim3,din1,din2,din3,dimz
    integer :: i,j,k,i1,j1,k1,is,js,n1,n2
    real(8) :: var1,dismin
    integer :: perc
    real(8) :: wei(-1:1,-1:1),dis(-1:1,-1:1)
    real(8),allocatable :: temp(:,:)
    !
    !$ SAVE wei,dis,temp
    !$OMP THREADPRIVATE(wei,dis,temp)
    !
    dim1=size(va1im,1)-1
    dim2=size(va1im,2)-1
    dim3=size(va1im,3)-1
    !
    din1=size(va1in,1)-1
    din2=size(va1in,2)-1
    din3=size(va1in,3)-1
    !
    if(din3 .ne. dim3) then
      print*,' !! the third dimension not match',dim3,din3
      stop 
    else
      dimz=dim3
    endif
    !
    print*,dim1,dim2,dim3
    print*,din1,din2,din3
    !
    !$OMP parallel default(shared) private(i,j,i1,j1,k1,is,js,n1,n2,   &
    !$OMP                                              var1,dismin,perc)
    allocate(temp(5,0:dimz))
    !
    !$OMP do
    do j1=0,din2
    do i1=0,din1
      !
      dismin=1d10
      do j=0,dim2
      do i=0,dim1
        !
        var1=(x1(i,j)-x2(i1,j1))**2+(y1(i,j)-y2(i1,j1))**2
        if(var1<dismin) then
          dismin=var1
          is=i
          js=j
        end if
        !
      end do
      end do
      !
      if(dismin<=1.d-15) then
        !
        va1in(i1,j1,:)=va1im(is,js,:)
        va2in(i1,j1,:)=va2im(is,js,:)
        va3in(i1,j1,:)=va3im(is,js,:)
        va4in(i1,j1,:)=va4im(is,js,:)
        va5in(i1,j1,:)=va5im(is,js,:)
        !
      else
        !
        var1=0.d0
        do n1=-1,1
        do n2=-1,1
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            wei(n1,n2)=0.d0
          else
            dis(n1,n2)=(x2(i1,j1)-x1(is+n1,js+n2))**2+                 &
                       (y2(i1,j1)-y1(is+n1,js+n2))**2
            wei(n1,n2)=1.d0/dis(n1,n2)
          end if
          !
          var1=var1+wei(n1,n2)
        end do
        end do
        !
        do n1=-1,1
        do n2=-1,1
          wei(n1,n2)=wei(n1,n2)/var1
        end do
        end do
        !
        temp=0.d0
        do n1=-1,1
        do n2=-1,1
          !
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2) then
            continue
          else
            do k=0,dimz
              temp(1,k)=temp(1,k)+wei(n1,n2)*va1im(is+n1,js+n2,k)
              temp(2,k)=temp(2,k)+wei(n1,n2)*va2im(is+n1,js+n2,k)
              temp(3,k)=temp(3,k)+wei(n1,n2)*va3im(is+n1,js+n2,k)
              temp(4,k)=temp(4,k)+wei(n1,n2)*va4im(is+n1,js+n2,k)
              temp(5,k)=temp(5,k)+wei(n1,n2)*va5im(is+n1,js+n2,k)
            enddo
          end if
          !
        end do
        end do
        !
        va1in(i1,j1,:)=temp(1,:)
        va2in(i1,j1,:)=temp(2,:)
        va3in(i1,j1,:)=temp(3,:)
        va4in(i1,j1,:)=temp(4,:)
        va5in(i1,j1,:)=temp(5,:)
        !
      end if
      !
    end do
      !
      perc=100*(j1)/(din2)
      write(*,'(1A1,A23,I3,A2,$)')char(13),'  ** Interpolating ... ',perc,' %'
      !
    end do
    !$OMP end do
    !$OMP end parallel
    !
  end subroutine InvdisInterp2D5_km
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The end of the subroutine InvdisInterp2D.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine InvdisInterp3D6(x1,y1,z1,va1im,va2im,va3im,va4im,va5im,      &
                             va6im,dim1,dim2,dim3,  &
                             x2,y2,z2,va1in,va2in,va3in,va4in,va5in,      &
                             va6in,din1,din2,din3   )
    !
    integer,intent(in) :: dim1,dim2,dim3,din1,din2,din3
    real(8),intent(in) :: x1(0:dim1,0:dim2,0:dim3),y1(0:dim1,0:dim2,0:dim3),z1(0:dim1,0:dim2,0:dim3), &
                          x2(0:din1,0:din2,0:din3),y2(0:din1,0:din2,0:din3),z2(0:din1,0:din2,0:din3)
    real(8),intent(in) :: va1im(0:dim1,0:dim2,0:dim3),va2im(0:dim1,0:dim2,0:dim3),   &
                          va3im(0:dim1,0:dim2,0:dim3),va4im(0:dim1,0:dim2,0:dim3),   &
                          va5im(0:dim1,0:dim2,0:dim3),va6im(0:dim1,0:dim2,0:dim3)
    real(8),intent(inout) :: va1in(0:din1,0:din2,0:din3),va2in(0:din1,0:din2,0:din3),&
                             va3in(0:din1,0:din2,0:din3),va4in(0:din1,0:din2,0:din3),&
                             va5in(0:din1,0:din2,0:din3),va6in(0:din1,0:din2,0:din3)
    !
    integer :: i,j,k,i1,j1,k1,is,js,ks,n1,n2,n3
    real(8) :: var1,dismin
    real(4) :: perc
    real(8) :: wei(-1:1,-1:1,-1:1),dis(-1:1,-1:1,-1:1),temp(6)
    !
    !!!$ SAVE wei,dis,temp
    !
    open(13,file='interp.log')
    !
    !!!$OMP THREADPRIVATE(wei,dis,temp)
    !
    !$OMP parallel default(shared) private(i,j,k,i1,j1,k1,is,js,ks,n1,n2,n3,   &
    !$OMP                                          var1,dismin,perc,wei,dis,temp)
    !$OMP do
    do k1=0,din3
    do j1=0,din2
    do i1=0,din1
      !
      dismin=1d10
      do k=0,dim3
      do j=0,dim2
      do i=0,dim1
        !
        var1=(x1(i,j,k)-x2(i1,j1,k1))**2+(y1(i,j,k)-y2(i1,j1,k1))**2+(z1(i,j,k)-z2(i1,j1,k1))**2
        if(var1<dismin) then
          dismin=var1
          is=i
          js=j
          ks=k
        end if
        !
      end do
      end do
      end do
      !
      if(dismin<=1.d-15) then
        !
        va1in(i1,j1,k1)=va1im(is,js,ks)
        va2in(i1,j1,k1)=va2im(is,js,ks)
        va3in(i1,j1,k1)=va3im(is,js,ks)
        va4in(i1,j1,k1)=va4im(is,js,ks)
        va5in(i1,j1,k1)=va5im(is,js,ks)
        va6in(i1,j1,k1)=va6im(is,js,ks)
        !
      else
        !
        var1=0.d0
        do n3=-1,1
        do n2=-1,1
        do n1=-1,1
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2 .or. ks+n3<0 .or. ks+n3>dim3) then
            wei(n1,n2,n3)=0.d0
          else
            dis(n1,n2,n3)=(x2(i1,j1,k1)-x1(is+n1,js+n2,ks+n3))**2+                 &
                          (y2(i1,j1,k1)-y1(is+n1,js+n2,ks+n3))**2+                 &
                          (z2(i1,j1,k1)-z1(is+n1,js+n2,ks+n3))**2
            wei(n1,n2,n3)=1.d0/sqrt(dis(n1,n2,n3))
          end if
          !
          var1=var1+wei(n1,n2,n3)
        end do
        end do
        end do
        !
        do n3=-1,1
        do n2=-1,1
        do n1=-1,1
          wei(n1,n2,n3)=wei(n1,n2,n3)/var1
        end do
        end do
        end do
        !
        temp=0.d0
        do n3=-1,1
        do n2=-1,1
        do n1=-1,1
          !
          if(is+n1<0 .or. is+n1>dim1 .or. js+n2<0 .or. js+n2>dim2 .or. ks+n3<0 .or. ks+n3>dim3) then
            continue
          else
            temp(1)=temp(1)+wei(n1,n2,n3)*va1im(is+n1,js+n2,ks+n3)
            temp(2)=temp(2)+wei(n1,n2,n3)*va2im(is+n1,js+n2,ks+n3)
            temp(3)=temp(3)+wei(n1,n2,n3)*va3im(is+n1,js+n2,ks+n3)
            temp(4)=temp(4)+wei(n1,n2,n3)*va4im(is+n1,js+n2,ks+n3)
            temp(5)=temp(5)+wei(n1,n2,n3)*va5im(is+n1,js+n2,ks+n3)
            temp(6)=temp(6)+wei(n1,n2,n3)*va6im(is+n1,js+n2,ks+n3)
          end if
          !
        end do
        end do
        end do
        !
        va1in(i1,j1,k1)=temp(1)
        va2in(i1,j1,k1)=temp(2)
        va3in(i1,j1,k1)=temp(3)
        va4in(i1,j1,k1)=temp(4)
        va5in(i1,j1,k1)=temp(5)
        va6in(i1,j1,k1)=temp(6)
        !
      end if
      !
    end do
      !
      perc=100.*(j1)/(din2)
      ! write(*,'(1A1,A23,F5.1,A2,$)')char(13),'  ** Interpolating ... ',perc,' %'
      write(13,*)'  ** Interpolating ... ',perc,' %'
      flush(13)
      !
    end do
    end do
    !$OMP end do
    !$OMP end parallel
    !
    close(13)
    !
  end subroutine InvdisInterp3D6
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to search a point in a structure mesh.         |
  !| return the number of i,j,k where the point is belong.             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-Aug-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine ijksearch(p,ip,jp,kp,lfound,initial,xin,yin,zin)
    !
    ! arguments
    real(8),intent(in),dimension(0:,0:,0:),optional :: xin,yin,zin
    real(8),intent(in) :: p(3)
    logical,intent(in) :: initial
    integer,intent(out) :: ip,jp,kp
    logical,intent(out) :: lfound
    !
    ! local data
    type :: subdomaim
      integer :: level,i0,j0,k0,im,jm,km
      integer :: next(8)
      real(8) :: xmin(3),xmax(3)
    end type subdomaim
    !
    real(8),allocatable,save :: x(:,:,:),y(:,:,:),z(:,:,:)
    type(subdomaim),allocatable,save :: domain(:)
    integer,allocatable,save :: ndomai(:)
    integer,save :: nlevel,nalldom
    !
    integer :: ilevel,im,jm,km,ic,id,i,j,k,in,jn,kn,                   &
               ncou,nco2,imm,jmm,kmm,is,js,ks,ie,je,ke,isub,n
    integer,allocatable :: idsaver(:),idstemp(:)
    !
    if(initial) then
      ! divide the domain into many subdomains using a tree structure
      !
      if(present(xin) .and. present(yin) .and. present(zin)) then
        im=size(xin,1)-1
        jm=size(xin,2)-1
        km=size(xin,3)-1
        !
        allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
        x=xin
        y=yin
        z=zin
      else
        stop 'xin,yin,zin missing'
      endif
      !
      nlevel=1
      do while(2**nlevel<max(im,jm,km))
        nlevel=nlevel+1
      enddo
      !
      nlevel=nlevel-2
      !
      allocate(ndomai(1:nlevel))
      ndomai(1)=1
      nalldom=1
      do ilevel=2,nlevel
        ndomai(ilevel)=ndomai(ilevel-1)*8
        nalldom=nalldom+ndomai(ilevel)
      enddo
      !
      print*,' ** number of levels  :',nlevel
      print*,' ** number of domains :',nalldom
      !
      allocate(domain(nalldom))
      !
      domain(1)%level=1
      domain(1)%i0=0; domain(1)%im=im
      domain(1)%j0=0; domain(1)%jm=jm
      domain(1)%k0=0; domain(1)%km=km
      !
      ic=1
      id=1
      ncou=1
      !
      do while(ncou<nalldom)
        !
        ! go through allocated domains
        do i=ic,id
          !
          nco2=0
          do kn=1,2
          do jn=1,2
          do in=1,2
            !
            ncou=ncou+1
            !
            nco2=nco2+1
            !
            ! assgin the children domains
            domain(i)%next(nco2)=ncou
            !
            imm=(domain(i)%im-domain(i)%i0)/2
            jmm=(domain(i)%jm-domain(i)%j0)/2
            kmm=(domain(i)%km-domain(i)%k0)/2
            !
            ! set the children domains
            domain(ncou)%level=domain(i)%level+1
            domain(ncou)%next=-1
            !
            if(in==1) then
              domain(ncou)%i0=domain(i)%i0
              domain(ncou)%im=domain(i)%i0 + imm
            elseif(in==2) then
              domain(ncou)%i0=domain(i)%i0 + imm
              domain(ncou)%im=domain(i)%im
            endif
            !
            if(jn==1) then
              domain(ncou)%j0=domain(i)%j0
              domain(ncou)%jm=domain(i)%j0 + jmm
            elseif(jn==2) then
              domain(ncou)%j0=domain(i)%j0 + jmm
              domain(ncou)%jm=domain(i)%jm
            endif
            !
            if(kn==1) then
              domain(ncou)%k0=domain(i)%k0
              domain(ncou)%km=domain(i)%k0 + kmm
            elseif(kn==2) then
              domain(ncou)%k0=domain(i)%k0 + kmm
              domain(ncou)%km=domain(i)%km
            endif
            !
            !
          enddo
          enddo
          enddo
          !
        enddo
        !
        ic=id+1
        id=ncou
        !
        print*,' ** search domain: ',ic,id
        !
      enddo
      !
      ! set search x,y,z extent
      do id=1,nalldom
        !
        ! only check the outer layer of the domain
        !
        is=domain(id)%i0
        ie=domain(id)%im
        js=domain(id)%j0
        je=domain(id)%jm
        ks=domain(id)%k0
        ke=domain(id)%km
        !
        domain(id)%xmin= 1.d10
        domain(id)%xmax=-1.d10
        do k=ks,ke
        do j=js,je
        do i=is,ie
          !
          if(i==is .or. j==js .or. k==ks .or. i==ie .or. j==je .or. k==ke) then
            ! only search the outerlayer
            domain(id)%xmin(1)=min(domain(id)%xmin(1),x(i,j,k))
            domain(id)%xmin(2)=min(domain(id)%xmin(2),y(i,j,k))
            domain(id)%xmin(3)=min(domain(id)%xmin(3),z(i,j,k))
            !
            domain(id)%xmax(1)=max(domain(id)%xmax(1),x(i,j,k))
            domain(id)%xmax(2)=max(domain(id)%xmax(2),y(i,j,k))
            domain(id)%xmax(3)=max(domain(id)%xmax(3),z(i,j,k))
          endif
          !
        enddo
        enddo
        enddo
        !
        if(id==1) then
          print*,' ** root domain id:',id,' level: ',domain(id)%level
          print*,' **  i extent:',domain(id)%i0,' ~ ',domain(id)%im
          print*,' **  j extent:',domain(id)%j0,' ~ ',domain(id)%jm
          print*,' **  k extent:',domain(id)%k0,' ~ ',domain(id)%km
          print*,' **  x extent:',domain(id)%xmin(1),' ~ ',domain(id)%xmax(1)
          print*,' **  y extent:',domain(id)%xmin(2),' ~ ',domain(id)%xmax(2)
          print*,' **  z extent:',domain(id)%xmin(3),' ~ ',domain(id)%xmax(3)
        endif
        if(id==nalldom) then
          print*,' ** highest domain id:',id,' level: ',domain(id)%level
          print*,' **  i extent:',domain(id)%i0,' ~ ',domain(id)%im
          print*,' **  j extent:',domain(id)%j0,' ~ ',domain(id)%jm
          print*,' **  k extent:',domain(id)%k0,' ~ ',domain(id)%km
          print*,' **  x extent:',domain(id)%xmin(1),' ~ ',domain(id)%xmax(1)
          print*,' **  y extent:',domain(id)%xmin(2),' ~ ',domain(id)%xmax(2)
          print*,' **  z extent:',domain(id)%xmin(3),' ~ ',domain(id)%xmax(3)
        endif
        !
      enddo
      !
    endif
    !
    ! now the search the subdomain that possibly contain p
    !
    allocate(idsaver(ndomai(nlevel)),idstemp(ndomai(nlevel)))
    !
    id=1
    ilevel=1
    ncou=1
    idsaver(1)=1
    !
    if( p(1)>=domain(id)%xmin(1) .and. p(1)<=domain(id)%xmax(1) .and. &
        p(2)>=domain(id)%xmin(2) .and. p(2)<=domain(id)%xmax(2) .and. &
        p(3)>=domain(id)%xmin(3) .and. p(3)<=domain(id)%xmax(3) ) then
      !
      ! write(*,'(2(A,I0))')'  **  found p at domain: ',id,', level: ',domain(id)%level
      ! !
      ! print*,' **  x extent:',domain(id)%xmin(1),' ~ ',domain(id)%xmax(1)
      ! print*,' **  y extent:',domain(id)%xmin(2),' ~ ',domain(id)%xmax(2)
      ! print*,' **  z extent:',domain(id)%xmin(3),' ~ ',domain(id)%xmax(3)
      !
      lfound=.true.
      !
    else
      print*,' ** the point is not in the main domain.'
      !
      lfound=.false.
      !
      ip=0
      jp=0
      kp=0
      return
    endif
    !
    do while(ilevel<nlevel .and. lfound)
      !
      ! check the extent of the subdomain
      lfound=.false.
      !
      nco2=0
      idstemp=0
      !
      do n=1,ncou
        !
        id=idsaver(n)
        !
        do isub=1,8
          !
          ic=domain(id)%next(isub)
          !
          if( p(1)>=domain(ic)%xmin(1) .and. p(1)<=domain(ic)%xmax(1) .and. &
              p(2)>=domain(ic)%xmin(2) .and. p(2)<=domain(ic)%xmax(2) .and. &
              p(3)>=domain(ic)%xmin(3) .and. p(3)<=domain(ic)%xmax(3) ) then
            !
            ! write(*,'(2(A,I0))')'  **  found p at domain: ',ic,', level: ',domain(ic)%level
            ! !
            ! print*,' **  x extent:',domain(ic)%xmin(1),' ~ ',domain(ic)%xmax(1)
            ! print*,' **  y extent:',domain(ic)%xmin(2),' ~ ',domain(ic)%xmax(2)
            ! print*,' **  z extent:',domain(ic)%xmin(3),' ~ ',domain(ic)%xmax(3)
            !
            nco2=nco2+1
            !
            idstemp(nco2)=ic
            !
            lfound=.true.
            ilevel=domain(ic)%level
            !
          endif
          !
        enddo
        !
      enddo
      !
      ncou=nco2
      idsaver=idstemp
      !
    enddo
    !
    ! now find p
    do n=1,ncou
      !
      id=idsaver(n)
      !
      is=domain(id)%i0
      ie=domain(id)%im
      js=domain(id)%j0
      je=domain(id)%jm
      ks=domain(id)%k0
      ke=domain(id)%km
      !
      do k=domain(id)%k0+1,domain(id)%km
      do j=domain(id)%j0+1,domain(id)%jm
      do i=domain(id)%i0+1,domain(id)%im
        !
        if( p(1)>=x(i-1,j,k) .and. p(1)<=x(i,j,k) .and. &
            p(2)>=y(i,j-1,k) .and. p(2)<=y(i,j,k) .and. &
            p(3)>=z(i,j,k-1) .and. p(3)<=z(i,j,k) ) then
          !
          ip=i
          jp=j
          kp=k
          !
          ! print*,' ** found p at domain:',id,'i,j,k',i,j,k
          !
          return
          !
        endif
        !
      enddo
      enddo
      enddo
      !
    enddo
    !
  end subroutine ijksearch
  !+-------------------------------------------------------------------+
  !| The end of the subroutine ijksearch.                              |
  !+-------------------------------------------------------------------+
  !!
  subroutine binterp_2dmesh_1dmesh(x_src,y_src,v1_src,v2_src,v3_src, &
                     x_dst,y_dst,v1_dst,v2_dst,v3_dst)
    !
    real(8),intent(in) :: x_src(:,:),y_src(:,:),x_dst(:),y_dst(:)
    real(8),intent(in) :: v1_src(:,:),v2_src(:,:),v3_src(:,:)
    real(8),intent(out) :: v1_dst(:),v2_dst(:),v3_dst(:)
    !
    integer :: im,jm,in,i,j,i1
    real(8) :: var1,var2,var3,var4,var5
    !
    im=size(x_src,1)
    jm=size(x_src,2)
    in=size(x_dst,1)
    write(*,"(A,3(I0,A))")' ** ',im,'x',jm,' mesh ->',in,' mesh'
    !
    do j=2,jm
    do i=2,im
      !
      do i1=1,in
        !
        if( x_dst(i1) >= x_src(i-1,j) .and. &
            x_dst(i1) <= x_src(i,j)   .and. &
            y_dst(i1) >= y_src(i,j-1) .and. &
            y_dst(i1) <= y_src(i,j) ) then
          !
          var1=1.d0/(x_src(i,j)-x_src(i-1,j))/(y_src(i,j)-y_src(i,j-1))
          var2=(x_src(i,j)-x_dst(i1))  *(y_src(i,j)-y_dst(i1))
          var3=(x_dst(i1)-x_src(i-1,j))*(y_src(i,j)-y_dst(i1))
          var4=(x_src(i,j)-x_dst(i1))  *(y_dst(i1)-y_src(i,j-1))
          var5=(x_dst(i1)-x_src(i-1,j))*(y_dst(i1)-y_src(i,j-1))
          !
          v1_dst(i1)=v1_src(i-1,j-1)*var2+v1_src(i,j-1)*var3+v1_src(i-1,j)*var4+v1_src(i,j)*var5
          v2_dst(i1)=v2_src(i-1,j-1)*var2+v2_src(i,j-1)*var3+v2_src(i-1,j)*var4+v2_src(i,j)*var5
          v3_dst(i1)=v3_src(i-1,j-1)*var2+v3_src(i,j-1)*var3+v3_src(i-1,j)*var4+v3_src(i,j)*var5
          !
          v1_dst(i1)=v1_dst(i1)*var1
          v2_dst(i1)=v2_dst(i1)*var1
          v3_dst(i1)=v3_dst(i1)*var1
          !
          exit
          !
        endif
        !
      enddo
      !
    enddo
    enddo
    !
    return
    !
  end subroutine binterp_2dmesh_1dmesh
  !!
  subroutine binterp_2dmesh_1dmesh_homz(x_src,y_src,v1_src,v2_src,v3_src, &
                     x_dst,y_dst,v1_dst,v2_dst,v3_dst,km)
    !
    real(8),intent(in) :: x_src(:,:),y_src(:,:),x_dst(:),y_dst(:)
    real(8),intent(in) :: v1_src(:,:,:),v2_src(:,:,:),v3_src(:,:,:)
    real(8),intent(out) :: v1_dst(:,:),v2_dst(:,:),v3_dst(:,:)
    integer,intent(in) :: km
    !
    integer :: im,jm,in,i,j,i1,k
    real(8) :: var1,var2,var3,var4,var5
    !
    im=size(x_src,1)
    jm=size(x_src,2)
    in=size(x_dst,1)
    write(*,"(A,3(I0,A))")' ** ',im,'x',jm,' mesh ->',in,' mesh'
    !
    do j=2,jm
    do i=2,im
      !
      do i1=1,in
        !
        if( x_dst(i1) >= x_src(i-1,j) .and. &
            x_dst(i1) <= x_src(i,j)   .and. &
            y_dst(i1) >= y_src(i,j-1) .and. &
            y_dst(i1) <= y_src(i,j) ) then
          !
          var1=1.d0/(x_src(i,j)-x_src(i-1,j))/(y_src(i,j)-y_src(i,j-1))
          var2=(x_src(i,j)-x_dst(i1))  *(y_src(i,j)-y_dst(i1))
          var3=(x_dst(i1)-x_src(i-1,j))*(y_src(i,j)-y_dst(i1))
          var4=(x_src(i,j)-x_dst(i1))  *(y_dst(i1)-y_src(i,j-1))
          var5=(x_dst(i1)-x_src(i-1,j))*(y_dst(i1)-y_src(i,j-1))
          !
          do k=1,km+1
            v1_dst(i1,k)=v1_src(i-1,j-1,k)*var2+v1_src(i,j-1,k)*var3+v1_src(i-1,j,k)*var4+v1_src(i,j,k)*var5
            v2_dst(i1,k)=v2_src(i-1,j-1,k)*var2+v2_src(i,j-1,k)*var3+v2_src(i-1,j,k)*var4+v2_src(i,j,k)*var5
            v3_dst(i1,k)=v3_src(i-1,j-1,k)*var2+v3_src(i,j-1,k)*var3+v3_src(i-1,j,k)*var4+v3_src(i,j,k)*var5
            !
            v1_dst(i1,k)=v1_dst(i1,k)*var1
            v2_dst(i1,k)=v2_dst(i1,k)*var1
            v3_dst(i1,k)=v3_dst(i1,k)*var1
          enddo
          !
          exit
          !
        endif
        !
      enddo
      !
    enddo
    enddo
    !
    return
    !
  end subroutine binterp_2dmesh_1dmesh_homz
  !!
end module interpolation
!+---------------------------------------------------------------------+
!| The end of the module interpolation.                                |
!+---------------------------------------------------------------------+
