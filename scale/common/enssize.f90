program enssize
!=======================================================================
!
! [PURPOSE:] Print the ensemble size (for Makefile's use)
!
!=======================================================================
  use common_letkf, only: nbv
  implicit none

  write (*, '(I4.4)') nbv

end program enssize
