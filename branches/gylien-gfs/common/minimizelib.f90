MODULE minimizelib
!=======================================================================
!
! [PURPOSE:] Minimizer
!
! [ATTENTION:] This module requires 'lbfgs.f'
!
! [HISTORY:]
!   07/31/2002 Yuki HONDA        created
!   03/20/2009 Takemasa MIYOSHI  modified
!
!=======================================================================
  USE common
  IMPLICIT NONE

  INTEGER,PARAMETER :: nsave = 5
  REAL(r_dble),PARAMETER :: epsln = 1.0d-4
  INTEGER,SAVE,PRIVATE :: iter
  INTEGER,SAVE,PRIVATE :: iprint(2)
  INTEGER,SAVE,PRIVATE :: point
  LOGICAL,SAVE,PRIVATE :: diagco
  REAL(r_dble),SAVE,PRIVATE,ALLOCATABLE :: zdiag(:)
  REAL(r_dble),SAVE,PRIVATE,ALLOCATABLE :: zs(:), zy(:), zw(:)

CONTAINS
!-----------------------------------------------------------------------
! INITIALIZE_MINIMIZER: Initialize Minimizer parameters
!-----------------------------------------------------------------------
SUBROUTINE initialize_minimizer(vsiz)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: vsiz

  iter = 0

  iprint(1) = 1; iprint(2) = 0
  diagco = .FALSE.

  ALLOCATE( zdiag(1:vsiz) )
  ALLOCATE( zs(1:vsiz * nsave) )
  ALLOCATE( zy(1:vsiz * nsave) )
  !ALLOCATE( zw(1:vsiz * nsave) )
  ALLOCATE( zw(1:vsiz + 2 * nsave) )

  zdiag = 0.0d0
  zs    = 0.0d0
  zy    = 0.0d0
  zw    = 0.0d0

  RETURN
END SUBROUTINE initialize_minimizer
!-----------------------------------------------------------------------
! MINIMIZE: Minimize by VA15AD (LBFGS)
!      iflag >= 0 : normal termination
!            == - maxiter : reach maximum iteration
!            <  0 : abnormal termination except number of iteration reaches 
!                   maximum number
!-----------------------------------------------------------------------
SUBROUTINE minimize(vsiz, xctl, costf, costg, maxiter, iflag)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: vsiz
  INTEGER,INTENT(IN) :: maxiter
  INTEGER,INTENT(OUT) :: iflag
  REAL(r_size),INTENT(INOUT) :: xctl(vsiz)
  REAL(r_size),INTENT(IN) :: costf
  REAL(r_size),INTENT(IN) :: costg(vsiz)

  CALL va15ad(vsiz, nsave, xctl, costf, costg, diagco, zdiag, iprint,        &
    &         epsln, zs, zy, point, zw, iflag, iter)

  IF(iter > maxiter) THEN
    iflag = - maxiter
  END IF

  RETURN
END SUBROUTINE minimize
!-----------------------------------------------------------------------
! TERMINATE_MINIMIZER: Terminate Minimizer
!-----------------------------------------------------------------------
SUBROUTINE terminate_minimizer
  IMPLICIT NONE

  DEALLOCATE( zdiag )
  DEALLOCATE( zs )
  DEALLOCATE( zy )
  DEALLOCATE( zw )

  RETURN
END SUBROUTINE terminate_minimizer

END MODULE minimizelib
