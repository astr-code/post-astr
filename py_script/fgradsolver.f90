module fgradsolver
    !
    implicit none
    !
    real(8),parameter :: num2d3  =2.d0/3.d0,   num1d12 =1.d0/12.d0
    !
    contains
    !
    ! subroutine test(im,jm,f,df)
    !     !
    !     integer,intent(in) :: im,jm
    !     real(8),intent(in) :: f(:,:)
    !     real(8),intent(out) :: df(2,im,jm)
        
    !     df=0.d0
    !     !
    ! end subroutine test
    !
    subroutine grad_xy(im,jm,f,x,y,df)
      !
      ! arguments
      integer,intent(in) :: im,jm
      real(8),intent(in) :: f(:,:),x(:,:),y(:,:)
      real(8),intent(out) :: df(2,im,jm)
      !
      ! local data
      integer :: i,j
      real(8),allocatable :: vtemp(:)
      real(8),allocatable :: ddi(:,:,:),ddj(:,:,:)
      !
      allocate(ddi(2,im,jm),ddj(2,im,jm))
      !
      call gridjacobian_xy(x,y,ddi,ddj,im,jm)
      !
      allocate(vtemp(im))
      do j=1,jm
        vtemp=dfdi4c(f(:,j),im)
        df(1,:,j)=vtemp(:)*ddi(1,:,j)
        df(2,:,j)=vtemp(:)*ddi(2,:,j)
      end do
      deallocate(vtemp)
      !
      allocate(vtemp(jm))
      do i=1,im
        vtemp=dfdi4c(f(i,:),jm)
        df(1,i,:)=df(1,i,:)+vtemp(:)*ddj(1,i,:)
        df(2,i,:)=df(2,i,:)+vtemp(:)*ddj(2,i,:)
      end do
      deallocate(vtemp)
      !
      print*, ' ** 2d gradient calculated'
      !
    end subroutine grad_xy
    !
    !+-------------------------------------------------------------------+
    !|this subroutine calculates grid jacobian matrix.                   |
    !+-------------------------------------------------------------------+
    !|Writen by Jian Fang , 2020-04-18.                                  |
    !+-------------------------------------------------------------------+
    subroutine gridjacobian_xy(x,y,ddi,ddj,im,jm)
      !
      ! arguments
      integer,intent(in) :: im,jm
      real(8),intent(in) :: x(:,:),y(:,:)
      real(8),intent(out) :: ddi(1:2,1:im,1:jm), &
                             ddj(1:2,1:im,1:jm)
      !
      ! local data
      integer :: i,j
      real(8) :: var1
      real(8),allocatable :: dx(:,:,:),dy(:,:,:)
      !
      allocate(dx(2,im,jm),dy(2,im,jm))
      !
      do j=1,jm
        dx(1,:,j)=dfdi4c(x(:,j),im)
        dy(1,:,j)=dfdi4c(y(:,j),im)
      end do
      !
      do i=1,im
        dx(2,i,:)=dfdi4c(x(i,:),jm)
        dy(2,i,:)=dfdi4c(y(i,:),jm)
      end do
      !
      ! print*,im,jm
      ! print*,size(x,1),size(x,2)
      ! print*,size(ddi,1),size(ddi,2),size(ddi,3)
      do j=1,jm
      do i=1,im
        !
        var1=dx(1,i,j)*dy(2,i,j)-dx(2,i,j)*dy(1,i,j)
        !
        var1=1.d0/var1
        !
        ddi(1,i,j)= var1*dy(2,i,j)
        ddi(2,i,j)=-var1*dx(2,i,j)
        !
        ddj(1,i,j)=-var1*dy(1,i,j)
        ddj(2,i,j)= var1*dx(1,i,j)
        !
      end do
      end do
      !
      deallocate(dx,dy)
      !
      print*, ' ** Grid Jacobian matrix in a x-y plane is calculated'
      !
    end subroutine gridjacobian_xy
    !
    !+-------------------------------------------------------------------+
    !|these functions are used to get spatial derivatives                |
    !+-------------------------------------------------------------------+
    function dfdi4c(vin,dim) result(dfdi)
      !
      real(8) :: dfdi(dim)
      integer,intent(in) :: dim
      real(8),intent(in) :: vin(dim)
      !
      integer :: i
      !
      dfdi(1)=-0.5d0*vin(3)+2.d0*vin(2)-1.5d0*vin(1)
      dfdi(2)=0.5d0*(vin(3)-vin(1))
      !
      !
      do i=3,dim-2
        dfdi(i)  =num2d3*(vin(i+1)-vin(i-1))-                         &
                  num1d12*(vin(i+2)-vin(i-2))
      enddo
      dfdi(dim-1)=0.5d0*(vin(dim)-vin(dim-2))
      dfdi(dim)  =0.5d0*vin(dim-2)-2.d0*vin(dim-1)+1.5d0*vin(dim)
      !
      return
      !
    end function dfdi4c
    !
end module fgradsolver