!+---------------------------------------------------------------------+
!| This is the post-process porgram for fluid dynamics.                |
!+---------------------------------------------------------------------+
!| Initilized on 2015-11-11, Writen By Fang Jian, fangjian19@gmail.com |
!+---------------------------------------------------------------------+
program postastr
  !
  use postentrance
  !
  implicit none
  !
  integer :: npost
  !
  print*,'+=========================Statement=========================+'
  print*,'|                                                           |'
  print*,'|      This is the Post-Process Program for code ASTR       |'
  print*,'|                    Writen by Fang Jian                    |'
  print*,'|                   fangjian19@gmail.com                    |'
  print*,'|             <ASTR>  Copyright Resvered <2011>             |'
  print*,'|                                                           |'
  print*,'+=========================Statement=========================+'
  print*
  ! 
  call gotoprocess
  !
end program postastr