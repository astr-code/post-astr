!+---------------------------------------------------------------------+
!| This module contains subroutiens to solve linear  algebra.          |
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2022-05-08                                     |
!+---------------------------------------------------------------------+
module LinearAlgegra
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to solve the tridiagonal matrix using          |
  !| thomas algorithm.                                                 |
  !+-------------------------------------------------------------------+
  !| ref: 数值分析，第四版，北京航空航天出版社                            | 
  !+-------------------------------------------------------------------+
  !| Writen by Jian Fang, 2022-05-09                                   |
  !+-------------------------------------------------------------------+
  function tds(b,scheme,bf,opt) result(x)
    !
    ! arguments
    real(8),intent(in) :: b(:)
    real(8) :: x(size(b))
    real(8),intent(in),optional :: bf
    character(len=*),intent(in) :: scheme
    character(len=*),intent(in),optional :: opt
    !
    ! local data
    integer :: m,i
    real(8) :: var1
    real(8),allocatable,dimension(:),save :: c,d,p,q
    real(8),allocatable :: y(:)
    logical,save :: linit=.true.
    !
    m=size(b)
    !
    if(linit) then
      !
      print*,'m=',m
      !
      call matrix_coef(c,d,m,scheme,bf)
      !
      allocate(p(1:m),q(1:m-1))
      !
      p(1)=1.d0
      do i=1,m-1
        q(i)=c(i)/p(i)
        p(i+1)=1.d0-d(i+1)*q(i)
      enddo
      !
      linit=.false.
      !
    endif
    !
    allocate(y(1:m))
    !
    y(1)=b(1)/p(1)
    do i=2,m
      y(i)=(b(i)-d(i)*y(i-1))/p(i)
    enddo
    !
    x(m)=y(m)
    do i=m-1,1,-1
      x(i)=y(i)-q(i)*x(i+1)
    enddo
    !
    deallocate(y)
    !
    if(present(opt) .and. opt=='once') then
      !
      deallocate(p,q)
      linit=.true.
      !
      return
      !
    endif
    !
    return
    !
  end function tds
  !+-------------------------------------------------------------------+
  !| This subroutine is to solve the quasi-tridiagonal matrix using    |
  !| thomas algorithm.                                                 |
  !+-------------------------------------------------------------------+
  !| ref: 数值分析，第四版，北京航空航天出版社                            | 
  !+-------------------------------------------------------------------+
  !| Writen by Jian Fang, 2022-05-08                                   |
  !+-------------------------------------------------------------------+
  function qtds(b,scheme,bf,opt) result(x)
    !
    ! arguments
    real(8),intent(in) :: b(:)
    real(8) :: x(size(b))
    real(8),intent(in),optional :: bf
    character(len=*),intent(in) :: scheme
    character(len=*),intent(in),optional :: opt
    !
    ! local data
    integer :: m,i
    real(8) :: var1
    real(8),allocatable,dimension(:),save :: c,d,p,q,s,r
    real(8),allocatable :: y(:)
    logical,save :: linit=.true.
    !
    m=size(b)
    !
    if(linit) then
      !
      print*,'m=',m
      !
      call matrix_coef(c,d,m,scheme,bf)
      !
      allocate(p(1:m-1),q(1:m-2),s(1:m-1),r(1:m))
      !
      p(1)=1.d0
      do i=1,m-2
        q(i)=c(i)/p(i)
        p(i+1)=1.d0-d(i+1)*q(i)
      enddo
      !
      s(1)=d(1)/p(1)
      do i=2,m-2
        s(i)=-d(i)*s(i-1)/p(i)
      enddo
      s(m-1)=(c(m-1)-d(m-1)*s(m-2))/p(m-1)
      !
      r(1)=c(m)
      do i=2,m-2
        r(i)=-r(i-1)*q(i-1)
      enddo
      r(m-1)=d(m)-r(m-2)*q(m-2)
      !
      r(m)=1.d0
      do i=1,m-1
        r(m)=r(m)-r(i)*s(i)
      enddo
      !
      linit=.false.
      !
    endif
    !
    allocate(y(1:m))
    !
    y(1)=b(1)/p(1)
    do i=2,m-1
      y(i)=(b(i)-d(i)*y(i-1))/p(i)
    enddo
    !
    y(m)=b(m)
    do i=1,m-1
      y(m)=y(m)-r(i)*y(i)
    enddo
    y(m)=y(m)/r(m)
    !
    x(m)=y(m)
    x(m-1)=y(m-1)-s(m-1)*x(m)
    do i=m-2,1,-1
      x(i)=y(i)-q(i)*x(i+1)-s(i)*x(m)
    enddo
    !
    deallocate(y)
    !
    if(present(opt) .and. opt=='once') then
      !
      deallocate(c,d,p,q,s,r)
      linit=.true.
      !
      return
      !
    endif
    !
    return
    !
  end function qtds
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to set A matrix.                               |
  !+-------------------------------------------------------------------+
  subroutine matrix_coef(c,d,m,scheme,bf)
    !
    integer,intent(in) :: m
    character(len=*),intent(in) :: scheme
    real(8),intent(in),optional :: bf
    real(8),allocatable,intent(out) :: c(:),d(:)
    !
    allocate(c(1:m),d(1:m))
    !
    if(scheme=='cu5') then
      !
      d(:)=0.5d0
      c(:)=1.d0/6.d0
      !
    elseif(scheme=='cu5-ld') then
      !
      d(:)=0.5d0-1.d0/6.d0*bf
      c(:)=1.d0/6.d0+1.d0/6.d0*bf
      !
    else
      stop ' !! scheme not defined @ qtds_a_init'
    endif
    !
    return
    !
  end subroutine matrix_coef
  !+-------------------------------------------------------------------+
  !| The end of the subroutine qtds  system                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is to do the cross product for a 3-D vector.        | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  pure function cross_product(a,b)
    !
    real(8) :: cross_product(3)
    real(8),dimension(3),intent(in) :: a(3), b(3)
  
    cross_product(1) = a(2) * b(3) - a(3) * b(2)
    cross_product(2) = a(3) * b(1) - a(1) * b(3)
    cross_product(3) = a(1) * b(2) - a(2) * b(1)
    !
    return
    !
  end function cross_product
  !+-------------------------------------------------------------------+
  !| The end of the function cross_product.                            |
  !+-------------------------------------------------------------------+
  !
end module LinearAlgegra
!+---------------------------------------------------------------------+
!| The end of the module LinearAlgegra                                 |
!+---------------------------------------------------------------------+