PROGRAM nature
  USE common
  USE lorenz96
!  USE lorenz96_oro

  IMPLICIT NONE

  INTEGER,PARAMETER :: ndays=3600 ! 10 years
  REAL(r_size) :: x(nx)
  REAL(r_sngl) :: x4(nx)
  INTEGER :: i,ktoneday
  INTEGER :: ktcyc

  dt=0.005d0
  force=8.0d0
  oneday=0.2d0

  ktoneday = INT(oneday/dt)
  ktcyc = ktoneday/4

  READ(10) x

  DO i=1,ndays*4
    x4 = x
    WRITE(90) x4
    CALL tinteg_rk4(ktcyc,x,x)
  END DO

  STOP
END PROGRAM nature
