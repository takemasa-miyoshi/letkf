PROGRAM obsmake
!=======================================================================
! simulate observation data
!=======================================================================
  USE common
  USE lorenz96
  USE h_ope

  IMPLICIT NONE

  INTEGER,PARAMETER :: ndays=3600
  INTEGER,PARAMETER :: nt=ndays*4
  REAL(r_size),PARAMETER :: obserr=1.0d0
  REAL(r_sngl) :: x4(nx)
  REAL(r_size) :: x(nx)
  REAL(r_size) :: ober(ny*nt)
  REAL(r_size) :: y(ny)
  REAL(r_sngl) :: y4(ny)
  INTEGER :: it
  INTEGER :: i,j,k

  CALL com_randn(ny*nt,ober)
  ober = ober * obserr
  k=0
  DO it=1,nt
    !
    ! nature run <- fort.10
    !
    READ(10) x4
    x = REAL(x4,r_size)
    CALL set_h(x)
    !
    ! Hx + ober
    !
    DO j=1,ny
      k = k+1
      y(j) = ober(k)
      DO i=1,nx
        y(j) = y(j) + h(j,i) * x(i)
      END DO
    END DO
    !
    ! obs data -> fort.91
    !
    y4 = y
    WRITE(91) y4
  END DO

  STOP
END PROGRAM obsmake
