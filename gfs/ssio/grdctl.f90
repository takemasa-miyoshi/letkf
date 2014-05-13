!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program grdctl

  use common_gfs
  use ssio_tools

  implicit none

  character(len=120) :: dset, options, t_start, t_int
  integer :: t_num

  character(len=120) :: comm_tmp
  integer :: comm_num, comm_len, comm_stat, iost

!-------------------------------------------------------------------------------

  comm_num = COMMAND_ARGUMENT_COUNT()
  if (comm_num < 6) then
    write (*, *)
    write (*, '(A)') 'Usage: ./grdctl DSET OPTIONS T_START T_INT T_NUM GRD_TYPE'
    write (*, *)
    write (*, '(A)') "  GRD_TYPE: 's': Sigma-level output"
    write (*, '(A)') "            'x': Extended sigma-level output"
    write (*, '(A)') "            'p': Pressure-level output"
    write (*, *)
    stop
  end if

  call GET_COMMAND_ARGUMENT (1, comm_tmp, comm_len, comm_stat)
  if (comm_stat /= 0) stop
  dset = comm_tmp(1:comm_len)

  call GET_COMMAND_ARGUMENT (2, comm_tmp, comm_len, comm_stat)
  if (comm_stat /= 0) stop
  options = comm_tmp(1:comm_len)

  call GET_COMMAND_ARGUMENT (3, comm_tmp, comm_len, comm_stat)
  if (comm_stat /= 0) stop
  t_start = comm_tmp(1:comm_len)

  call GET_COMMAND_ARGUMENT (4, comm_tmp, comm_len, comm_stat)
  if (comm_stat /= 0) stop
  t_int = comm_tmp(1:comm_len)

  call GET_COMMAND_ARGUMENT (5, comm_tmp, comm_len, comm_stat)
  if (comm_stat /= 0) stop
  read (comm_tmp(1:comm_len), '(I)', iostat=iost) t_num
  if (iost /= 0) stop

  call GET_COMMAND_ARGUMENT (6, comm_tmp, comm_len, comm_stat)
  if (comm_stat /= 0) stop

  call set_common_gfs
  select case (comm_tmp(1:comm_len))
    case ('s')
      call write_ctl(6, trim(dset), trim(options), trim(t_start), trim(t_int), t_num)
    case ('x')
      call write_ctlx(6, trim(dset), trim(options), trim(t_start), trim(t_int), t_num)
    case ('p')
      call write_ctlp(6, trim(dset), trim(options), trim(t_start), trim(t_int), t_num)
    case default
      stop
  end select

!-------------------------------------------------------------------------------

end program grdctl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
