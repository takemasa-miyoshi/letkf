MODULE h_ope
!=======================================================================
! observation operator
!=======================================================================
  USE common
  USE lorenz96

  IMPLICIT NONE

  INTEGER,PARAMETER :: ny=13
  REAL(r_size),SAVE :: h(ny,nx)

CONTAINS
SUBROUTINE set_h(x)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: x(nx)
  INTEGER :: i
  INTEGER :: idx(nx)

  idx( 1)= 1
  idx( 2)= 4
  idx( 3)= 7
  idx( 4)=10
  idx( 5)=13
  idx( 6)=16
  idx( 7)=19
  idx( 8)=22
  idx( 9)=25
  idx(10)=28
  idx(11)=31
  idx(12)=34
  idx(13)=37
  h = 0.0d0
  DO i=1,ny
    h(i,idx(i)) = 1.0d0
  END DO

  RETURN
END SUBROUTINE set_h

END MODULE h_ope
