subroutine h5_read1int(varin,vname,fname,explicit)
  !
  use hdf5
  use h5lt
  !
  integer,intent(out) :: varin
  character(len=*),intent(in) :: vname,fname
  logical,intent(in), optional:: explicit
  logical :: lexplicit
  !
  integer(hid_t) :: file_id
  ! file identifier
  integer(hid_t) :: dset_id1
  ! dataset identifier
  integer :: v(1)
  integer :: h5error ! error flag
  integer(hsize_t) :: dimt(1)
  !
  if (present(explicit)) then
     lexplicit = explicit
  else
     lexplicit = .true.
  end if
  !
  dimt=(/1/)
  !
  call h5open_f(h5error)
  print*,' ** open hdf5 interface'
  !
  call h5fopen_f(fname,h5f_acc_rdwr_f,file_id,h5error)
  !
  call h5ltread_dataset_f(file_id,vname,h5t_native_integer,v,dimt,h5error)
  !
  call h5fclose_f(file_id,h5error)
  !
  if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 1'
  !
  ! close fortran interface.
  call h5close_f(h5error)
  !
  varin=v(1)
  if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 2'
  !
  if(lexplicit)  print*,' >> ',vname,' from ',fname,' ... done'
  !
end subroutine h5_read1int