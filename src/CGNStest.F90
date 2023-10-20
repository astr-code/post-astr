!+---------------------------------------------------------------------+
!| This module contains subroutines to test CGNS lib.                  |
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2020-04-24                                     |
!+---------------------------------------------------------------------+
module cgnstest
  !
  use cgns
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| https://cgns.github.io/CGNS_docs_current/user/started.html        |
  !|   Single-Zone Structured Grid                                     |
  !+-------------------------------------------------------------------+
  subroutine cgnswrite_grid_str
    !
    !Creates simple 3-D structured grid and writes it to a
    !CGNS file.
    !
    !This program uses the fortran convention that all
    !variables beginning with the letters i-n are integers,
    !by default, and all others are real
    !
    !Example compilation for this program is (change paths!):
    !
    !ifort -I ../CGNS_CVS/cgnslib -c write_grid_str.f
    !ifort -o write_grid_str write_grid_str.o -L ../CGNS_CVS/cgnslib/LINUX -lcgns
    !
    !(../CGNS_CVS/cgnslib/LINUX/ is the location where the compiled
    !library libcgns.a is located)
    ! 
    ! cgnslib_f.h file must be located in directory specified by -I during compile:
    !
    !dimension statements (note that tri-dimensional arrays
    !x,y,z must be dimensioned exactly as (21,17,N) (N>=9)
    !for this particular case or else they will be written to
    !the CGNS file incorrectly!  Other options are to use 1-D
    !arrays, use dynamic memory, or pass index values to a
    !subroutine and dimension exactly there):
    real(8) :: x(21,17,9),y(21,17,9),z(21,17,9)
    integer(4) :: isize(3,3)
    integer :: ni,nj,nk,i,j,k
    integer :: icelldim,iphysdim,ier
    integer :: index_file,index_base,index_zone,index_coord
    character(len=32) :: basename,zonename
    !
    ni=21
    nj=17
    nk=9
    do k=1,nk
      do j=1,nj
        do i=1,ni
          x(i,j,k)=real(i-1,8)
          y(i,j,k)=real(j-1,8)
          z(i,j,k)=real(k-1,8)
        enddo
      enddo
    enddo
    print*,' ** created simple 3-D grid points using CGNS interface'
    !
    !   WRITE X, Y, Z GRID POINTS TO CGNS FILE
    !   open CGNS file for write
    call cg_open_f('grid.cgns',CG_MODE_WRITE,index_file,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
    ! create base (user can give any name)
    !
    basename='Base'
    icelldim=3
    iphysdim=3
    call cg_base_write_f(index_file,basename,icelldim,iphysdim,index_base,ier)
    ! define zone name (user can give any name)
    zonename = 'Zone  1'
    !
    ! vertex size
    isize(1,1)=21
    isize(2,1)=17
    isize(3,1)=9
    !
    ! cell size
    isize(1,2)=isize(1,1)
    isize(2,2)=isize(2,1)
    isize(3,2)=isize(3,1)
    !
    ! boundary vertex size (always zero for structured grids)
    isize(1,3)=0
    isize(2,3)=0
    isize(3,3)=0
    !
    ! create zone
    call cg_zone_write_f(index_file,index_base,zonename,isize,structured,index_zone,ier)
    !
    ! write grid coordinates (user must use SIDS-standard names here)
    call cg_coord_write_f(index_file,index_base,index_zone,realdouble, &
                                        'coordinatex',x,index_coord,ier)
    print*,index_coord
    call cg_coord_write_f(index_file,index_base,index_zone,RealDouble, &
                                        'CoordinateY',y,index_coord,ier)
    print*,index_coord
    call cg_coord_write_f(index_file,index_base,index_zone,RealDouble, &
                                        'CoordinateZ',z,index_coord,ier)
    print*,index_coord
    !
    print*,x(:,1,1)
    ! close CGNS file
    call cg_close_f(index_file,ier)
    !
    print*,' ** Successfully wrote grid to file grid.cgns'
    !
  end subroutine cgnswrite_grid_str
  !+-------------------------------------------------------------------+
  !| The end of the subroutine cgnswrite_grid_str.                     |
  !+-------------------------------------------------------------------+
  subroutine cgnread_grid_str
    !
    !Creates simple 3-D structured grid and writes it to a
    !CGNS file.
    !
    !This program uses the fortran convention that all
    !variables beginning with the letters i-n are integers,
    !by default, and all others are real
    !
    !Example compilation for this program is (change paths!):
    !
    !ifort -I ../CGNS_CVS/cgnslib -c write_grid_str.f
    !ifort -o write_grid_str write_grid_str.o -L ../CGNS_CVS/cgnslib/LINUX -lcgns
    !
    !(../CGNS_CVS/cgnslib/LINUX/ is the location where the compiled
    !library libcgns.a is located)
    ! 
    ! cgnslib_f.h file must be located in directory specified by -I during compile:
    !
    !dimension statements (note that tri-dimensional arrays
    !x,y,z must be dimensioned exactly as (21,17,N) (N>=9)
    !for this particular case or else they will be written to
    !the CGNS file incorrectly!  Other options are to use 1-D
    !arrays, use dynamic memory, or pass index values to a
    !subroutine and dimension exactly there):
    real(8) :: x(21,17,9),y(21,17,9),z(21,17,9)
    integer :: irmin(3),irmax(3)
    integer(4) :: isize(3,3)
    integer :: ni,nj,nk,i,j,k
    integer :: icelldim,iphysdim,ier
    integer :: index_file,index_base,index_zone,index_coord
    character(len=32) :: basename,zonename
    !
    !   WRITE X, Y, Z GRID POINTS TO CGNS FILE
    !   open CGNS file for write
    call cg_open_f('grid.cgns',CG_MODE_READ,index_file,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
    ! create base (user can give any name)
    !
    ! we know there is only one base (real working code would check!)
    index_base=1
    ! we know there is only one zone (real working code would check!)
    index_zone=1
    !
    ! get zone size (and name - although not needed here)
    call cg_zone_read_f(index_file,index_base,index_zone,zonename,isize,ier)
    !
    ! lower range index
    irmin(1)=1
    irmin(2)=1
    irmin(3)=1
    !
    ! upper range index of vertices
    irmax(1)=isize(1,1)
    irmax(2)=isize(2,1)
    irmax(3)=isize(3,1)
    !
    ! read grid coordinates
    call cg_coord_read_f(index_file,index_base,index_zone,'coordinatex',RealSingle,irmin,irmax,x,ier)
    call cg_coord_read_f(index_file,index_base,index_zone,'coordinatey',RealSingle,irmin,irmax,y,ier)
    call cg_coord_read_f(index_file,index_base,index_zone,'coordinatez',RealSingle,irmin,irmax,z,ier)
    !
    print*,irmax
    print*,x(:,1,1)
    !
    ! close CGNS file
    call cg_close_f(index_file,ier)
    print*,' ** Successfully read grid from file grid.cgns'
    !
  end subroutine cgnread_grid_str
  !+-------------------------------------------------------------------+
  !| The end of the subroutine cgnswrite_grid_str.                     |
  !+-------------------------------------------------------------------+
  !
end module cgnstest
!+---------------------------------------------------------------------+
!| The end of the module cgnstest.                                     |
!+---------------------------------------------------------------------+