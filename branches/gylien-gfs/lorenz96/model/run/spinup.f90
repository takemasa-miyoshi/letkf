PROGRAM spinup
  USE common
  USE lorenz96
!  USE lorenz96_oro

  IMPLICIT NONE

  REAL(r_size) :: x(nx)
  INTEGER :: i,ktoneday

  dt=0.005d0
  force=8.0d0
  oneday=0.2d0

  ktoneday = INT(oneday/dt)

  CALL com_randn(nx,x)
  x = x * 5.0d0

  CALL tinteg_rk4(ktoneday*360*100,x,x) ! 100 years integration

  WRITE(90) x

  STOP
END PROGRAM spinup
