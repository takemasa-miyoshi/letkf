MODULE common_enkf
!=======================================================================
!
! [PURPOSE:] Ensemble Kalman Filter (EnKF)
!            Model Independent Core Module
!
! [REFERENCES:]
!  Miyoshi, T., 2004: Application of ensemble Kalman filtering
!    including model bias estimation to a primitive-equation global
!    model. PhD Proposal, University of Maryland, College Park, 27pp.
!
! [HISTORY:]
!  07/20/2004 Takemasa Miyoshi  Created at University of Maryland, College Park
!  02/16/2005 Takemasa Miyoshi  Bias analysis (enkf_banl0) is added
!
!=======================================================================
  USE common
  USE common_mtx

  IMPLICIT NONE

  PUBLIC

CONTAINS
!=======================================================================
! Analysis
!=======================================================================
!-----------------------------------------------------------------------
! Analysis for ny=1
!  K = E (HE)^T [ HE (HE)^T + R ]^(-1)
!-----------------------------------------------------------------------
SUBROUTINE enkf_anal0(nx,m,hdxf,r,dxf,gain)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx               ! dimension of analysis space
  INTEGER,INTENT(IN) :: m                ! ensemble size
  REAL(r_size),INTENT(IN) :: hdxf(m)     ! HE
  REAL(r_size),INTENT(IN) :: r           ! observational error variance
  REAL(r_size),INTENT(IN) :: dxf(nx,m)   ! background ensemble ptb
  REAL(r_size),INTENT(OUT) :: gain(nx)   ! Kalman gain
  REAL(r_size) :: wk1d1(m)
  REAL(r_size) :: wk0d1
!
! [ HE (HE)^T + R ] -> wk0d1
!
  wk0d1 = SUM(hdxf**2) / REAL(m-1,r_size) + r
!
! K = E (HE)^T [ HE (HE)^T + R ]^(-1)
!
  wk1d1 = hdxf / wk0d1
  gain = MATMUL(dxf,wk1d1) / REAL(m-1,r_size)

  RETURN
END SUBROUTINE enkf_anal0
!-----------------------------------------------------------------------
! Bias Analysis for ny=1
!  Kb = Eb (HEb)^T [ HEb (HEb)^T + HE (HE)^T + R ]^(-1)
!-----------------------------------------------------------------------
SUBROUTINE enkf_bias_anal0(nx,m,hdbf,hdxf,r,dbf,gain)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx               ! dimension of analysis space
  INTEGER,INTENT(IN) :: m                ! ensemble size
  REAL(r_size),INTENT(IN) :: hdbf(m)     ! HEb
  REAL(r_size),INTENT(IN) :: hdxf(m)     ! HE
  REAL(r_size),INTENT(IN) :: r           ! observational error variance
  REAL(r_size),INTENT(IN) :: dbf(nx,m)   ! background bias ensemble ptb
  REAL(r_size),INTENT(OUT) :: gain(nx)   ! Kalman gain
  REAL(r_size) :: wk1d1(m)
  REAL(r_size) :: wk0d1
!
! [ HEb (HEb)^T + HE (HE)^T + R ] -> wk0d1
!
  wk0d1 = ( SUM(hdbf**2) + SUM(hdxf**2) ) / REAL(m-1,r_size) + r
!
! Kb = Eb (HEb)^T [ HEb (HEb)^T + HE (HE)^T + R ]^(-1)
!
  wk1d1 = hdbf / wk0d1
  gain = MATMUL(dbf,wk1d1) / REAL(m-1,r_size)

  RETURN
END SUBROUTINE enkf_bias_anal0
!-----------------------------------------------------------------------
! Analysis for small ny; usually used serially
!  Matrix inversion of the size [ny x ny] is required
!  K = E (HE)^T [ HE (HE)^T + R ]^(-1)
!-----------------------------------------------------------------------
SUBROUTINE enkf_anal1(nx,ny,m,hdxf,dep,r,dxf,xfbar,xabar,gain)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx               ! dimension of analysis space
  INTEGER,INTENT(IN) :: ny               ! number of observation
  INTEGER,INTENT(IN) :: m                ! ensemble size
  REAL(r_size),INTENT(IN) :: hdxf(ny,m)  ! HE
  REAL(r_size),INTENT(IN) :: dep(ny)     ! departure (y-Hx)
  REAL(r_size),INTENT(IN) :: r(ny,ny)    ! observational error covariance
  REAL(r_size),INTENT(IN) :: dxf(nx,m)   ! background ensemble ptb
  REAL(r_size),INTENT(IN) :: xfbar(nx)   ! background ensemble mean
  REAL(r_size),INTENT(OUT) :: xabar(nx)  ! analysis ensemble mean
  REAL(r_size),INTENT(OUT) :: gain(nx,ny)! Kalman gain
  REAL(r_size) :: wk2d1(ny,ny),wk2d2(ny,ny),wk2d3(m,ny)
!
! [ HE (HE)^T + R ]^(-1) -> wk2d2
!
  wk2d1 = matmul(hdxf,transpose(hdxf)) / REAL(m-1,r_size) + r
  CALL mtx_inv(ny,wk2d1,wk2d2)
!
! K = E (HE)^T [ HE (HE)^T + R ]^(-1)
!
  wk2d3 = matmul(transpose(hdxf),wk2d2)
  gain = matmul(dxf,wk2d3) / REAL(m-1,r_size)
!
! xa = xf + K * dep -> xabar (RETURN)
!
  xabar = xfbar + matmul(gain,dep)

  RETURN
END SUBROUTINE enkf_anal1
!-----------------------------------------------------------------------
! Analysis for not small ny assuming diagonal R
!  Matrix inversion of the size [m x m] is required
!  K = E [ I + (HE)^T R^(-1) (HE) ]^(-1) (HE)^T R^(-1)
!-----------------------------------------------------------------------
SUBROUTINE enkf_anal2(nx,ny,m,hdxf,dep,rdiag,dxf,xfbar,xabar)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx               ! dimension of analysis space
  INTEGER,INTENT(IN) :: ny               ! number of observation
  INTEGER,INTENT(IN) :: m                ! ensemble size
  REAL(r_size),INTENT(IN) :: hdxf(ny,m)  ! HE
  REAL(r_size),INTENT(IN) :: dep(ny)     ! departure (y-Hx)
  REAL(r_size),INTENT(IN) :: rdiag(ny)   ! observational error variance
  REAL(r_size),INTENT(IN) :: dxf(nx,m)   ! background ensemble ptb
  REAL(r_size),INTENT(IN) :: xfbar(nx)   ! background ensemble mean
  REAL(r_size),INTENT(OUT) :: xabar(nx)  ! analysis ensemble mean
  REAL(r_size) :: wk1d1(ny),wk1d2(m),wk1d3(m)
  REAL(r_size) :: wk2d1(ny,m),wk2d2(m,m),wk2d3(m,m)
  INTEGER :: i
!
! [ I + (HE)^T R^(-1) (HE) ]^(-1) -> wk2d3
!
  DO i=1,ny
    wk2d1(i,:) = hdxf(i,:) / rdiag(i)
  END DO
  wk2d2 = matmul(transpose(hdxf),wk2d1) / REAL(m-1,r_size)
  DO i=1,m
    wk2d2(i,i) = wk2d2(i,i) + 1.0d0
  END DO
  CALL mtx_inv(m,wk2d2,wk2d3)
!
! dxa =  E [ I + (HE)^T R^(-1) (HE) ]^(-1) (HE)^T R^(-1) * dep -> xabar
!
  wk1d1 = dep / rdiag
  wk1d2 = matmul(transpose(hdxf),wk1d1)
  wk1d3 = matmul(wk2d3,wk1d2)
  xabar = matmul(dxf,wk1d3) / REAL(m-1,r_size)
!
! xa = xf + dxa -> xabar (RETURN)
!
  xabar = xfbar + xabar

  RETURN
END SUBROUTINE enkf_anal2
!=======================================================================
! Ensemble Update
!=======================================================================
!-----------------------------------------------------------------------
! Serial method: Eq(13) of Whitaker and Hamill (2002) (Potter 1964)
! ny = 1
!-----------------------------------------------------------------------
SUBROUTINE enkf_serial(nx,m,hdxf,r,dxf,gain,dxa)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx               ! dimension of analysis space
  INTEGER,INTENT(IN) :: m                ! ensemble size
  REAL(r_size),INTENT(IN) :: hdxf(m)     ! HE
  REAL(r_size),INTENT(IN) :: r           ! observational error variance
  REAL(r_size),INTENT(IN) :: dxf(nx,m)   ! background ensemble ptb
  REAL(r_size),INTENT(IN) :: gain(nx)    ! Kalman gain
  REAL(r_size),INTENT(OUT) :: dxa(nx,m)  ! analysis ensemble ptb
  REAL(r_size) :: alpha
  INTEGER :: i,j

  alpha = 1 + SQRT( r / ( SUM(hdxf**2)/REAL(m-1,r_size) + r ) )

  DO j=1,m
    DO i=1,nx
      dxa(i,j) = dxf(i,j) - gain(i) * hdxf(j) / alpha
    END DO
  END DO

  RETURN
END SUBROUTINE enkf_serial
!-----------------------------------------------------------------------
! Serial method for bias analysis
! ny = 1
!-----------------------------------------------------------------------
SUBROUTINE enkf_bias_serial(nx,m,hdbf,hdxf,r,dbf,gain,dba)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx               ! dimension of analysis space
  INTEGER,INTENT(IN) :: m                ! ensemble size
  REAL(r_size),INTENT(IN) :: hdbf(m)     ! HEb
  REAL(r_size),INTENT(IN) :: hdxf(m)     ! HE
  REAL(r_size),INTENT(IN) :: r           ! observational error variance
  REAL(r_size),INTENT(IN) :: dbf(nx,m)   ! background bias ensemble ptb
  REAL(r_size),INTENT(IN) :: gain(nx)    ! Kalman gain
  REAL(r_size),INTENT(OUT) :: dba(nx,m)  ! analysis bias ensemble ptb
  REAL(r_size) :: alpha
  INTEGER :: i,j

  alpha = 1 + SQRT( r / ( (SUM(hdbf**2)+SUM(hdxf**2))/REAL(m-1,r_size) + r ) )

  DO j=1,m
    DO i=1,nx
      dba(i,j) = dbf(i,j) - gain(i) * hdbf(j) / alpha
    END DO
  END DO

  RETURN
END SUBROUTINE enkf_bias_serial
!-----------------------------------------------------------------------
! ETKF
!-----------------------------------------------------------------------
SUBROUTINE enkf_etkf(nx,ny,m,hdxf,rdiag,dxf,dxa)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx               ! dimension of analysis space
  INTEGER,INTENT(IN) :: ny               ! number of observation
  INTEGER,INTENT(IN) :: m                ! ensemble size
  REAL(r_size),INTENT(IN) :: hdxf(ny,m)  ! HE
  REAL(r_size),INTENT(IN) :: rdiag(ny)   ! observational error variance
  REAL(r_size),INTENT(IN) :: dxf(nx,m)   ! background ensemble ptb
  REAL(r_size),INTENT(OUT) :: dxa(nx,m)  ! analysis ensemble ptb
  REAL(r_size) :: wk1d1(m),wk2d1(ny,m),wk2d2(m,m),wk2d3(m,m)
  INTEGER :: i,nrank_eff
!
! R^(-1/2) (HE) -> wk2d1
!
  DO i=1,ny
    wk2d1(i,:) = hdxf(i,:) / SQRT( rdiag(i) )
  END DO
!
! [ R^(-1/2) (HE) ]^T [ R^(-1/2) (HE) ] -> wk2d2
!
  wk2d2 = matmul(transpose(wk2d1),wk2d1) / REAL(m-1,r_size)
!
! Eigenvalue decomposition CGC^T
!
  CALL mtx_eigen(1,m,wk2d2,wk1d1,wk2d3,nrank_eff)
!
! T = C [ G + I ]^(-1/2) where CGC^T is eigenvalue decomposition
!
  DO i=1,m
    wk2d2(:,i) = wk2d3(:,i)  / SQRT( wk1d1(i) + 1.0d0 )
  END DO
!
! Ea = Ef T
!
  dxa = matmul(dxf,wk2d2)

  RETURN
END SUBROUTINE enkf_etkf
!=======================================================================
! Localization
!=======================================================================
!-----------------------------------------------------------------------
! Schur Product
!-----------------------------------------------------------------------
SUBROUTINE enkf_schur(scale,dist,factor)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: scale
  REAL(r_size),INTENT(IN) :: dist
  REAL(r_size),INTENT(OUT) :: factor
  REAL(r_size) :: a,b

  a = scale * SQRT(10.0d0/3.0d0)
  b = dist / a

  IF( dist <= a ) THEN
    factor = 1.0d0 -0.25d0*b**5 + 0.5d0*b**4 + 5.0d0/8.0d0*b**3 &
      & - 5.0d0/3.0d0*b**2
  ELSE IF( dist <= 2*a ) THEN
    factor = 1.0d0/12.0d0*b**5 - 0.5d0*b**4 + 5.0d0/8.0d0*b**3 &
      & + 5.0d0/3.0d0*b**2 - 5.0d0*b + 4.0d0 - 2.0d0/3.0d0/b
  ELSE
    factor = 0.0d0
  END IF

  RETURN
END SUBROUTINE enkf_schur

END MODULE common_enkf
