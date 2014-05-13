program datetime
!=======================================================================
!
! [PURPOSE:] Date and time computation
!
!=======================================================================
  use common, only: com_datetime_reg, com_mdays

  implicit none
  integer :: yyyy, mm, dd, hh, ii, ss, inc
  character (len=1) :: sunit

  read (*, '(I)') yyyy
  read (*, '(I)') mm
  read (*, '(I)') dd
  read (*, '(I)') hh
  read (*, '(I)') ii
  read (*, '(I)') ss
  read (*, '(I)') inc
  read (*, '(A1)') sunit

  select case (sunit)
  case ('y')
    yyyy = yyyy + inc
  case ('m')
    mm = mm + inc
  case ('d')
    dd = dd + inc
  case ('h')
    hh = hh + inc
  case ('i')
    ii = ii + inc
  case default
    ss = ss + inc
  end select

  call com_datetime_reg(yyyy, mm, dd, hh, ii, ss)

  write (*, '(I4.4,5I2.2)') yyyy, mm, dd, hh, ii, ss

end program datetime
