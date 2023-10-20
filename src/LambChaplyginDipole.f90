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