MODULE common_letkf
!=======================================================================
!
! [PURPOSE:] Local Ensemble Transform Kalman Filtering (LETKF)
!            Model Independent Core Module
!
! [REFERENCES:]
!  [1] Ott et al., 2004: A local ensemble Kalman filter for atmospheric
!    data assimilation. Tellus, 56A, 415-428.
!  [2] Hunt et al., 2007: Efficient Data Assimilation for Spatiotemporal
!    Chaos: A Local Ensemble Transform Kalman Filter. Physica D, 230,
!    112-126.
!
! [HISTORY:]
!  01/21/2009 Takemasa Miyoshi  Created at U. of Maryland, College Park
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mtx

  IMPLICIT NONE

  PUBLIC
!=======================================================================
!  LEKF Model Independent Parameters
!=======================================================================
  INTEGER,PARAMETER :: nbv=20    ! ensemble size

CONTAINS
!=======================================================================
!  Main Subroutine of LETKF Core
!   INPUT
!     ne               : ensemble size                                           !GYL
!     nobs             : array size, but only first nobsl elements are used
!     nobsl            : total number of observation assimilated at the point
!     hdxb(nobs,ne)    : obs operator times fcst ens perturbations
!     rdiag(nobs)      : observation error variance
!     rloc(nobs)       : localization weighting function
!     dep(nobs)        : observation departure (yo-Hxb)
!     parm_infl        : covariance inflation parameter
!     rdiag_wloc       : (optional) flag indicating that rdiag = rdiag / rloc    !GYL
!     minfl            : (optional) minimum covariance inflation parameter       !GYL
!   OUTPUT
!     parm_infl        : updated covariance inflation parameter
!     trans(ne,ne)     : transformation matrix
!     transm(ne)       : (optional) transformation matrix mean                   !GYL
!     pao(ne,ne)       : (optional) analysis covariance matrix in ensemble space !GYL
!=======================================================================
SUBROUTINE letkf_core(ne,nobs,nobsl,hdxb,rdiag,rloc,dep,parm_infl,trans,transm,pao,rdiag_wloc,minfl)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ne                      !GYL
  INTEGER,INTENT(IN) :: nobs
  INTEGER,INTENT(IN) :: nobsl
  REAL(r_size),INTENT(IN) :: hdxb(1:nobs,1:ne)
  REAL(r_size),INTENT(IN) :: rdiag(1:nobs)
  REAL(r_size),INTENT(IN) :: rloc(1:nobs)
  REAL(r_size),INTENT(IN) :: dep(1:nobs)
  REAL(r_size),INTENT(INOUT) :: parm_infl
  REAL(r_size),INTENT(OUT) :: trans(ne,ne)
  REAL(r_size),INTENT(OUT),OPTIONAL :: transm(ne)
  REAL(r_size),INTENT(OUT),OPTIONAL :: pao(ne,ne)
  LOGICAL,INTENT(IN),OPTIONAL :: rdiag_wloc     !GYL
  REAL(r_size),INTENT(IN),OPTIONAL :: minfl     !GYL

  REAL(r_size) :: hdxb_rinv(nobsl,ne)
  REAL(r_size) :: eivec(ne,ne)
  REAL(r_size) :: eival(ne)
  REAL(r_size) :: pa(ne,ne)
  REAL(r_size) :: work1(ne,ne)
  REAL(r_size) :: work2(ne,nobsl)
  REAL(r_size) :: work3(ne)
  REAL(r_size) :: rho
  REAL(r_size) :: parm(4),sigma_o,gain
  REAL(r_size),PARAMETER :: sigma_b = 0.04d0 !error stdev of parm_infl
  LOGICAL :: rdiag_wloc_
  INTEGER :: i,j,k

  rdiag_wloc_ = .FALSE.                            !GYL
  IF(present(rdiag_wloc)) rdiag_wloc_ = rdiag_wloc !GYL

  IF(nobsl == 0) THEN
    trans = 0.0d0
    DO i=1,ne
      trans(i,i) = SQRT(parm_infl)
    END DO
    IF (PRESENT(transm)) THEN   !GYL
      transm = 0.0d0            !GYL
    END IF                      !GYL
    IF (PRESENT(pao)) THEN                        !GYL
      pao = 0.0d0                                 !GYL
      DO i=1,ne                                   !GYL
        pao(i,i) = parm_infl / REAL(ne-1,r_size)  !GYL
      END DO                                      !GYL
    END IF                                        !GYL
    RETURN
  ELSE
!-----------------------------------------------------------------------
!  hdxb Rinv
!-----------------------------------------------------------------------
  IF(rdiag_wloc_) THEN                            !GYL
    DO j=1,ne                                     !GYL
      DO i=1,nobsl                                !GYL
        hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i)     !GYL
      END DO                                      !GYL
    END DO                                        !GYL
  ELSE                                            !GYL
    DO j=1,ne
      DO i=1,nobsl
        hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i) * rloc(i)
      END DO
    END DO
  END IF                                          !GYL
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb
!-----------------------------------------------------------------------
  CALL dgemm('t','n',ne,ne,nobsl,1.0d0,hdxb_rinv,nobsl,hdxb(1:nobsl,:),&
    & nobsl,0.0d0,work1,ne)
!  DO j=1,ne
!    DO i=1,ne
!      work1(i,j) = hdxb_rinv(1,i) * hdxb(1,j)
!      DO k=2,nobsl
!        work1(i,j) = work1(i,j) + hdxb_rinv(k,i) * hdxb(k,j)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb + (m-1) I / rho (covariance inflation)
!-----------------------------------------------------------------------
  IF (PRESENT(minfl)) THEN                           !GYL
    IF (minfl > 0.0d0 .AND. parm_infl < minfl) THEN  !GYL
      parm_infl = minfl                              !GYL
    END IF                                           !GYL
  END IF                                             !GYL
  rho = 1.0d0 / parm_infl
  DO i=1,ne
    work1(i,i) = work1(i,i) + REAL(ne-1,r_size) * rho
  END DO
!-----------------------------------------------------------------------
!  eigenvalues and eigenvectors of [ hdxb^T Rinv hdxb + (m-1) I ]
!-----------------------------------------------------------------------
  CALL mtx_eigen(1,ne,work1,eival,eivec,i)
!-----------------------------------------------------------------------
!  Pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv
!-----------------------------------------------------------------------
  DO j=1,ne
    DO i=1,ne
      work1(i,j) = eivec(i,j) / eival(j)
    END DO
  END DO
  CALL dgemm('n','t',ne,ne,ne,1.0d0,work1,ne,eivec,&
    & ne,0.0d0,pa,ne)
!  DO j=1,ne
!    DO i=1,ne
!      pa(i,j) = work1(i,1) * eivec(j,1)
!      DO k=2,ne
!        pa(i,j) = pa(i,j) + work1(i,k) * eivec(j,k)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T
!-----------------------------------------------------------------------
  CALL dgemm('n','t',ne,nobsl,ne,1.0d0,pa,ne,hdxb_rinv,&
    & nobsl,0.0d0,work2,ne)
!  DO j=1,nobsl
!    DO i=1,ne
!      work2(i,j) = pa(i,1) * hdxb_rinv(j,1)
!      DO k=2,ne
!        work2(i,j) = work2(i,j) + pa(i,k) * hdxb_rinv(j,k)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
  DO i=1,ne
    work3(i) = work2(i,1) * dep(1)
    DO j=2,nobsl
      work3(i) = work3(i) + work2(i,j) * dep(j)
    END DO
  END DO
!-----------------------------------------------------------------------
!  T = sqrt[(m-1)Pa]
!-----------------------------------------------------------------------
  DO j=1,ne
    rho = SQRT( REAL(ne-1,r_size) / eival(j) )
    DO i=1,ne
      work1(i,j) = eivec(i,j) * rho
    END DO
  END DO
  CALL dgemm('n','t',ne,ne,ne,1.0d0,work1,ne,eivec,&
    & ne,0.0d0,trans,ne)
!  DO j=1,ne
!    DO i=1,ne
!      trans(i,j) = work1(i,1) * eivec(j,1)
!      DO k=2,ne
!        trans(i,j) = trans(i,j) + work1(i,k) * eivec(j,k)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  T + Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
  IF (PRESENT(transm)) THEN                !GYL - if transm is present,
    transm = work3                         !GYL - return both trans and transm without adding them
  ELSE                                     !GYL
    DO j=1,ne
      DO i=1,ne
        trans(i,j) = trans(i,j) + work3(i)
      END DO
    END DO
  END IF                                   !GYL
  IF (PRESENT(pao)) pao = pa               !GYL
!-----------------------------------------------------------------------
!  Inflation estimation
!-----------------------------------------------------------------------
  parm = 0.0d0
  IF(rdiag_wloc_) THEN                            !GYL
    DO i=1,nobsl                                  !GYL
      parm(1) = parm(1) + dep(i)*dep(i)/rdiag(i)  !GYL
    END DO                                        !GYL
  ELSE                                            !GYL
    DO i=1,nobsl
      parm(1) = parm(1) + dep(i)*dep(i)/rdiag(i) * rloc(i)
    END DO
  END IF                                          !GYL
  DO j=1,ne
    DO i=1,nobsl
      parm(2) = parm(2) + hdxb_rinv(i,j) * hdxb(i,j)
    END DO
  END DO
  parm(2) = parm(2) / REAL(ne-1,r_size)
  parm(3) = SUM(rloc(1:nobsl))
  parm(4) = (parm(1)-parm(3))/parm(2) - parm_infl
!  sigma_o = 1.0d0/REAL(nobsl,r_size)/MAXVAL(rloc(1:nobsl))
  sigma_o = 2.0d0/parm(3)*((parm_infl*parm(2)+parm(3))/parm(2))**2
  gain = sigma_b**2 / (sigma_o + sigma_b**2)
  parm_infl = parm_infl + gain * parm(4)

  RETURN
  END IF
END SUBROUTINE letkf_core

END MODULE common_letkf
