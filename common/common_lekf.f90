MODULE common_lekf
!=======================================================================
!
! [PURPOSE:] Local Ensemble Kalman Filtering (LEKF)
!    and     Local Ensemble Transform Kalman Filtering (LETKF)
!            Model Independent Core Module
!
! [ATTENTION:]
!  + Comments
!    - this module considers only local 1D-array, allowing model independence
!    - local patch projection should be done beforehand
!    - LEKF and LETKF has the same interface, easy to switch
!  + Input to LEKF/LETKF
!    - background ensemble perturbations
!    - background ensemble mean
!    - observation increment (a.k.a. departure)
!    - observation error variance
!    - covariance inflation parameter
!  + Output by LEKF/LETKF
!    - analyzed ensemble
!
! [REFERENCES:]
!  [1] Ott et al., 2004: A local ensemble Kalman filter for atmospheric
!    data assimilation. Tellus, 56A, 415-428.
!  [2] Hunt, 2005: Efficient Data Assimilation for Spatiotemporal Chaos: 
!    a Local Ensemble Transform Kalman Filter. arXiv:physics/0511236, 
!    http:arxiv.org/abs/physics/0511236
!
! [TESTED WITH:]
!  LORENZ-96 model, SPEEDY model, AFES model
!
! [HISTORY:]
!  10/31/2003 Takemasa Miyoshi  Created at University of Maryland, College Park
!  10/31/2004 Takemasa Miyoshi  Version up by reference [1]
!  10/07/2005 Takemasa Miyoshi  Avoid ALLOCATE and MATMUL
!  10/19/2005 Takemasa Miyoshi  Use H^T instead of H for vector optimization
!  11/16/2005 Takemasa Miyoshi  Bug fix, this bug was very important
!  11/16/2005 Takemasa Miyoshi  Add LETKF
!  12/02/2005 Takemasa Miyoshi  More efficient LETKF by reference [2]
!  12/02/2005 Takemasa Miyoshi  Change interfaces (xb -> dxb, xb_bar)
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mtx

  IMPLICIT NONE

  PRIVATE
!=======================================================================
!  LEKF Model Independent Parameters
!=======================================================================
  INTEGER,SAVE :: nbv    ! number of BV members
  INTEGER,SAVE :: nrank  ! rank of internal coordinate system
  INTEGER,SAVE :: ndiml  ! total dimension of state variable in local patch
  INTEGER,SAVE :: nobsl  ! total number of observation in local patch
  LOGICAL,SAVE :: msw_inflate_enhanced=.TRUE. ! enhanced inflation
  LOGICAL,SAVE :: msw_orthogonal=.TRUE. ! Consideration of Orthogonal Component

!$OMP THREADPRIVATE(nrank,ndiml,nobsl)
  PUBLIC :: nbv, ndiml, nobsl
  PUBLIC :: letkf_core, lekf_core, set_common_lekf, calc_edim_local, calc_expvar_local

CONTAINS
!=======================================================================
!  Set LEKF Parameters
!=======================================================================
SUBROUTINE set_common_lekf
  INTEGER :: ios

  NAMELIST/common_lekf/nbv,msw_inflate_enhanced,msw_orthogonal
  OPEN(10,FILE='parm.txt',IOSTAT=ios)
  READ(10,NML=common_lekf,IOSTAT=ios)
  CLOSE(10)
  IF( ios/=0 ) THEN
    PRINT *,'WARNING (set_common_lekf): NAMELIST READ ERROR ',ios
    STOP 1
  END IF

  RETURN
END SUBROUTINE set_common_lekf
!=======================================================================
!  Main Subroutine of LETKF Core
!   INPUT
!     ht(ndiml,nobsl)   : transpose of linear observation operator
!     rdiag(nobsl)      : observation error variance
!     dep(nobsl)        : observation departure (yo-Hxb)
!     dxb(ndiml,nbv)    : background ensemble perturbations
!     xb_bar(ndiml,nbv) : background ensemble mean
!     parm_infl         : covariance inflation parameter
!   OUTPUT
!     xa(ndiml,nbv)     : analysis ensembles
!=======================================================================
SUBROUTINE letkf_core(ht,rdiag,dep,dxb,xb_bar,xa,parm_infl)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: ht(1:ndiml,1:nobsl)
  REAL(r_size),INTENT(IN) :: rdiag(1:nobsl)
  REAL(r_size),INTENT(IN) :: dep(1:nobsl)
  REAL(r_size),INTENT(IN) :: dxb(1:ndiml,1:nbv)
  REAL(r_size),INTENT(IN) :: xb_bar(1:ndiml)
  REAL(r_size),INTENT(OUT) :: xa(1:ndiml,1:nbv)
  REAL(r_size),INTENT(IN) :: parm_infl
  REAL(r_size) :: hdxb(nobsl,nbv)
  REAL(r_size) :: hdxb_rinv(nobsl,nbv)
  REAL(r_size) :: eivec(nbv,nbv)
  REAL(r_size) :: eival(nbv)
  REAL(r_size) :: pa(nbv+1,nbv)
  REAL(r_size) :: trans(nbv,nbv)
  REAL(r_size) :: work1(nbv,nbv)
  REAL(r_size) :: work2(nbv+1,nobsl)
  REAL(r_size) :: work3(nbv)
  REAL(r_size) :: rho
  INTEGER :: i,j,k
!-----------------------------------------------------------------------
!  H dxb
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nobsl
      hdxb(i,j) = ht(1,i) * dxb(1,j)
      DO k=2,ndiml
        hdxb(i,j) = hdxb(i,j) + ht(k,i) * dxb(k,j)
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
!  hdxb Rinv
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nobsl
      hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i)
    END DO
  END DO
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nbv
      work1(i,j) = hdxb_rinv(1,i) * hdxb(1,j)
      DO k=2,nobsl
        work1(i,j) = work1(i,j) + hdxb_rinv(k,i) * hdxb(k,j)
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb + (m-1) I / rho (covariance inflation)
!-----------------------------------------------------------------------
  rho = 1.0d0 / (1.0d0 + parm_infl) / (1.0d0 + parm_infl)
  DO i=1,nbv
    work1(i,i) = work1(i,i) + REAL(nbv-1,r_size) * rho
  END DO
!-----------------------------------------------------------------------
!  eigenvalues and eigenvectors of [ hdxb^T Rinv hdxb + (m-1) I ]
!-----------------------------------------------------------------------
  CALL mtx_eigen(1,nbv,work1,eival,eivec,i)
!-----------------------------------------------------------------------
!  Pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nbv
      work1(i,j) = eivec(i,j) / eival(j)
    END DO
  END DO
  DO j=1,nbv
    DO i=1,nbv
      pa(i,j) = work1(i,1) * eivec(j,1)
      DO k=2,nbv
        pa(i,j) = pa(i,j) + work1(i,k) * eivec(j,k)
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T
!-----------------------------------------------------------------------
  DO j=1,nobsl
    DO i=1,nbv
      work2(i,j) = pa(i,1) * hdxb_rinv(j,1)
      DO k=2,nbv
        work2(i,j) = work2(i,j) + pa(i,k) * hdxb_rinv(j,k)
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
  DO i=1,nbv
    work3(i) = work2(i,1) * dep(1)
    DO j=2,nobsl
      work3(i) = work3(i) + work2(i,j) * dep(j)
    END DO
  END DO
!-----------------------------------------------------------------------
!  T = sqrt[(m-1)Pa]
!-----------------------------------------------------------------------
  DO j=1,nbv
    rho = SQRT( REAL(nbv-1,r_size) / eival(j) )
    DO i=1,nbv
      work1(i,j) = eivec(i,j) * rho
    END DO
  END DO
  DO j=1,nbv
    DO i=1,nbv
      trans(i,j) = work1(i,1) * eivec(j,1)
      DO k=2,nbv
        trans(i,j) = trans(i,j) + work1(i,k) * eivec(j,k)
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
!  T + Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nbv
      trans(i,j) = trans(i,j) + work3(i)
    END DO
  END DO
!-----------------------------------------------------------------------
!  xa = xb_bar + dxb T
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,ndiml
      xa(i,j) = xb_bar(i)
      DO k=1,nbv
        xa(i,j) = xa(i,j) + dxb(i,k) * trans(k,j)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE letkf_core
!=======================================================================
!  Main Subroutine of LEKF Core
!   INPUT
!     ht(ndiml,nobsl)   : transpose of linear observation operator
!     rdiag(nobsl)      : observation error variance
!     dep(nobsl)        : observation departure (yo-Hxb)
!     dxb(ndiml,nbv)    : background ensemble perturbations
!     xb_bar(ndiml,nbv) : background ensemble mean
!     parm_infl         : covariance inflation parameter
!   OUTPUT
!     xa(ndiml,nbv)     : analysis ensembles
!=======================================================================
SUBROUTINE lekf_core(ht,rdiag,dep,dxb,xb_bar,xa,parm_infl)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: ht(1:ndiml,1:nobsl)
  REAL(r_size),INTENT(IN) :: rdiag(1:nobsl)
  REAL(r_size),INTENT(IN) :: dep(1:nobsl)
  REAL(r_size),INTENT(INOUT) :: dxb(1:ndiml,1:nbv)
  REAL(r_size),INTENT(IN) :: xb_bar(1:ndiml)
  REAL(r_size),INTENT(OUT) :: xa(1:ndiml,1:nbv)
  REAL(r_size),INTENT(IN) :: parm_infl
  REAL(r_size) :: xa_bar(1:ndiml)
  REAL(r_size) :: dxa(1:ndiml,1:nbv)
  REAL(r_size) :: lambda(1:nbv)
  REAL(r_size) :: q(1:ndiml,1:nbv)
  REAL(r_size) :: y(1:nbv,1:nbv)
  REAL(r_size) :: dxb_hat(nbv-1,nbv)
  REAL(r_size) :: dxa_hat(nbv-1,nbv)
  REAL(r_size) :: pa_hat(nbv-1,nbv-1)
  REAL(r_size) :: h_hat(nobsl,nbv-1)
  REAL(r_size) :: dxa_bar(nbv-1)
  INTEGER :: i,j,k
!-----------------------------------------------------------------------
!  Covariance Inflation (Regular Covariance Inflation)
!-----------------------------------------------------------------------
  IF(parm_infl /= 0) THEN
    IF(.NOT.msw_inflate_enhanced) THEN
      dxb = dxb * (1.0d0 + parm_infl)
    ENDIF
  ENDIF
!-----------------------------------------------------------------------
!  eigenvectors and eigenvalues
!  i.e. Q and lambda = diag(Pb^)
!-----------------------------------------------------------------------
  CALL calc_pb(dxb,lambda,q,nrank)
!  CALL monit_pb_hat(lambda)
!-----------------------------------------------------------------------
!  H^ = H Q
!-----------------------------------------------------------------------
!  h_hat(1:nobsl,1:nrank) = MATMUL(h(1:nobsl,1:ndiml),q(1:ndiml,1:nrank))
  DO j=1,nrank
!CDIR NOVECTOR
    DO i=1,nobsl
      h_hat(i,j) = ht(1,i)*q(1,j)
      DO k=2,ndiml
        h_hat(i,j) = h_hat(i,j) + ht(k,i)*q(k,j)
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
!  Pa^ = Pb^ [I + H^T Rinv H^ Pb^]inv -- Eq.(23)
!-----------------------------------------------------------------------
  CALL calc_pa(lambda,h_hat,rdiag,pa_hat)
!  CALL monit_p_hat(pa_hat)
!-----------------------------------------------------------------------
! Ott's Enhanced Variance Inflation Method -- Eq.(42)
!-----------------------------------------------------------------------
  IF(parm_infl /= 0) THEN
    IF(msw_inflate_enhanced) THEN  ! Ott's Enhanced Variance Inflation
      CALL inflate(pa_hat(1:nrank,1:nrank),parm_infl)
    ENDIF
  ENDIF
!-----------------------------------------------------------------------
!  dxa^_bar = Pa^ H^T Rinv * departure -- Eq.(22)
!-----------------------------------------------------------------------
  CALL calc_dxa_bar(pa_hat,h_hat,rdiag,dep,dxa_bar)
!-----------------------------------------------------------------------
!! Analysis state !!
!  xa_bar = xb_bar + Q dxa^_bar -- Eq.(24)
!-----------------------------------------------------------------------
!  xa_bar = xb_bar + MATMUL(q(1:ndiml,1:nrank),dxa_bar(1:nrank))
  DO i=1,ndiml
    xa_bar(i) = xb_bar(i) + q(i,1) * dxa_bar(1)
    DO k=2,nrank
      xa_bar(i) = xa_bar(i) + q(i,k) * dxa_bar(k)
    END DO
  END DO
!-----------------------------------------------------------------------
!  dxb^ = QT dxb
!-----------------------------------------------------------------------
!  dxb_hat = MATMUL(TRANSPOSE(q(1:ndiml,1:nrank)),dxb(1:ndiml,1:nbv))
  DO j=1,nbv
!CDIR NOVECTOR
    DO i=1,nrank
      dxb_hat(i,j) = q(1,i)*dxb(1,j)
      DO k=2,ndiml
        dxb_hat(i,j) = dxb_hat(i,j) + q(k,i)*dxb(k,j)
      END DO
    END DO
  END DO
!  CALL monit_p_hat( MATMUL(dxb_hat,transpose(dxb_hat))/REAL(nbv-1,r_size) )
!-----------------------------------------------------------------------
!  Y (Definition of Y: dxa^ = dxb^ Y) -- Eq.(41)
!-----------------------------------------------------------------------
  CALL calc_y(dxb_hat,lambda,pa_hat,y)
!-----------------------------------------------------------------------
!  dxa^ = dxb^ Y -- Eq.(34)
!-----------------------------------------------------------------------
!  dxa_hat = MATMUL(dxb_hat(1:nrank,1:nbv),y(1:nbv,1:nbv))
  DO j=1,nbv
    DO i=1,nrank
      dxa_hat(i,j) = dxb_hat(i,1)*y(1,j)
      DO k=2,nbv
        dxa_hat(i,j) = dxa_hat(i,j) + dxb_hat(i,k)*y(k,j)
      END DO
    END DO
  END DO
!  CALL monit_p_hat( MATMUL(dxa_hat,transpose(dxa_hat))/REAL(nbv-1,r_size) )
!-----------------------------------------------------------------------
!  dxa = Q dxa^
!-----------------------------------------------------------------------
!  dxa = MATMUL(q(1:ndiml,1:nrank),dxa_hat(1:nrank,1:nbv))
  DO j=1,nbv
    DO i=1,ndiml
      dxa(i,j) = q(i,1)*dxa_hat(1,j)
      DO k=2,nrank
        dxa(i,j) = dxa(i,j) + q(i,k)*dxa_hat(k,j)
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
!  dxa = Q dxa^ + [I - Q QT] dxb ! Orthogonal Component Eq.(14),(15),(16)
!   [ATTN:] This part requires [ndiml x ndiml] matrix, which may be huge!
!-----------------------------------------------------------------------
  IF(msw_orthogonal) THEN
    dxa = dxa + dxb &
      & - MATMUL(MATMUL(q(1:ndiml,1:nrank),TRANSPOSE(q(1:ndiml,1:nrank))),dxb)
  END IF
!-----------------------------------------------------------------------
!  xa = xa_bar + dxa
!-----------------------------------------------------------------------
  DO i=1,nbv
    xa(:,i) = xa_bar(:) + dxa(:,i)
  END DO

  RETURN
END SUBROUTINE lekf_core
!=======================================================================
!  Compute Eigenvectors and Eigenvalues
!  which gives transformation matrix Q and diagonal components of Pb
!    diag(Pb^) = eigenvalues
!    Q = (eivenvectors)
!=======================================================================
SUBROUTINE calc_pb(dxb,lambda,q,nrank_eff)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: dxb(1:ndiml,1:nbv)
  REAL(r_size),INTENT(OUT) :: lambda(1:nbv)
  REAL(r_size),INTENT(OUT) :: q(1:ndiml,1:nbv)
  INTEGER,INTENT(OUT) :: nrank_eff
  REAL(r_size) :: pb(1:nbv,1:nbv)
  REAL(r_size) :: eivec(1:nbv,1:nbv)
  INTEGER :: i,j,k

!  pb = MATMUL(TRANSPOSE(dxb),dxb) / REAL(nbv-1,r_size)
  DO j=1,nbv
!CDIR NOVECTOR
    DO i=1,nbv
      pb(i,j) = dxb(1,i) * dxb(1,j)
      DO k=2,ndiml
        pb(i,j) = pb(i,j) + dxb(k,i)*dxb(k,j)
      END DO
    END DO
  END DO
  pb = pb / REAL(nbv-1,r_size)

  CALL mtx_eigen(1,nbv,pb,lambda,eivec,nrank_eff)

  IF( nrank_eff < nbv) THEN
    lambda(nrank_eff+1:nbv) = 0.0d0
    eivec(:,nrank_eff+1:nbv) = 0.0d0
  END IF

!  q = MATMUL(dxb(1:ndiml,1:nbv),eivec(1:nbv,1:nbv))
  DO j=1,nrank_eff
    DO i=1,ndiml
      q(i,j) = dxb(i,1)*eivec(1,j)
      DO k=2,nbv
        q(i,j) = q(i,j) + dxb(i,k)*eivec(k,j)
      END DO
    END DO
  END DO

  DO i=1,nrank_eff
    q(:,i) = q(:,i) / SQRT( SUM( q(:,i)**2 ) )
  END DO

  RETURN
END SUBROUTINE calc_pb
!=======================================================================
!  Compute Eq.(23) ! analysis covariabce matrix on hat space
!    Pa^ = Pb^ [I + H^T Rinv H^ Pb^]inv
!=======================================================================
SUBROUTINE calc_pa(lambda,h_hat,rdiag,pa_hat)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: lambda(1:nbv)
  REAL(r_size),INTENT(IN) :: h_hat(1:nobsl,1:nbv-1)
  REAL(r_size),INTENT(IN) :: rdiag(1:nobsl)
  REAL(r_size),INTENT(OUT) :: pa_hat(1:nbv-1,1:nbv-1)
!  REAL(r_size) :: pb_hat(1:nrank,1:nrank)
!  REAL(r_size) :: rinv(1:nobsl,1:nobsl)
  REAL(r_size) :: work_h(1:nobsl,1:nrank)
  REAL(r_size) :: work_obs(1:nobsl,1:nrank)
  REAL(r_size) :: work_hat(1:nrank,1:nrank)
  REAL(r_size) :: work_hat_inv(1:nrank,1:nrank)
  INTEGER :: i,j,k
!-----------------------------------------------------------------------
!  Pb^
!-----------------------------------------------------------------------
!  pb_hat(:,:) = 0.0d0
!  DO i=1,nrank
!    pb_hat(i,i) = lambda(i)
!  END DO
!-----------------------------------------------------------------------
!  Rinv
!-----------------------------------------------------------------------
!-----
! For Future Use (in case that observation error correlation is not zero)
!
!  CALL mtx_inv(nobsl,r,rinv)
!-----
!  rinv = 0.0d0
!  DO i=1,nobsl
!    rinv(i,i) = 1.0d0 / rdiag(i)
!  END DO
!-----------------------------------------------------------------------
!  H^ Pb^
!-----------------------------------------------------------------------
!  work_h(1:nobsl,1:nrank) = MATMUL(h_hat(1:nobsl,1:nrank),pb_hat(1:nrank,1:nrank))
  DO j=1,nrank
    DO i=1,nobsl
!      work_h(i,j) = h_hat(i,1)*pb_hat(1,j)
!      DO k=2,nrank
!        work_h(i,j) = work_h(i,j) + h_hat(i,k)*pb_hat(k,j)
!      END DO
      work_h(i,j) = h_hat(i,j)*lambda(j)
    END DO
  END DO
!-----------------------------------------------------------------------
!  Rinv H^ Pb^
!-----------------------------------------------------------------------
!  work_obs(1:nobsl,1:nrank) = MATMUL(rinv(1:nobsl,1:nobsl),work_h(1:nobsl,1:nrank))
  DO j=1,nrank
    DO i=1,nobsl
!      work_obs(i,j) = rinv(i,1)*work_h(1,j)
!      DO k=2,nobsl
!        work_obs(i,j) = work_obs(i,j) + rinv(i,k)*work_h(k,j)
!      END DO
      work_obs(i,j) = work_h(i,j) / rdiag(i)
    END DO
  END DO
!-----------------------------------------------------------------------
!  H^T Rinv H^ Pb^
!-----------------------------------------------------------------------
!  work_hat(1:nrank,1:nrank) = MATMUL(TRANSPOSE(h_hat(1:nobsl,1:nrank)),work_obs(1:nobsl,1:nrank))
  DO j=1,nrank
    DO i=1,nrank
      work_hat(i,j) = h_hat(1,i)*work_obs(1,j)
      DO k=2,nobsl
        work_hat(i,j) = work_hat(i,j) + h_hat(k,i)*work_obs(k,j)
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
!  I + H^T Rinv H^ Pb^
!-----------------------------------------------------------------------
  DO i=1,nrank
    work_hat(i,i) = work_hat(i,i) + 1.0d0
  END DO
!-----------------------------------------------------------------------
!  [I + H^T Rinv H^ Pb^]inv
!-----------------------------------------------------------------------
  CALL mtx_inv(nrank,work_hat,work_hat_inv)
!-----------------------------------------------------------------------
!  Pb^ [I + H^T Rinv H^ Pb^]inv
!-----------------------------------------------------------------------
!  pa_hat = MATMUL(pb_hat,pa_hat)
  DO j=1,nrank
    DO i=1,nrank
!      pa_hat(i,j) = pb_hat(i,1)*work_hat_inv(1,j)
!      DO k=2,nrank
!        pa_hat(i,j) = pa_hat(i,j) + pb_hat(i,k)*work_hat_inv(k,j)
!      END DO
      pa_hat(i,j) = lambda(i)*work_hat_inv(i,j)
    END DO
  END DO

  RETURN
END SUBROUTINE calc_pa
!=======================================================================
!  Compute Eq.(22) : Analysis increment on hat space
!    dxa^ = Pa^ H^T Rinv * departure
!=======================================================================
SUBROUTINE calc_dxa_bar(pa_hat,h_hat,rdiag,dep,dxa_bar)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: pa_hat(1:nbv-1,1:nbv-1)
  REAL(r_size),INTENT(IN) :: h_hat(1:nobsl,1:nbv-1)
  REAL(r_size),INTENT(IN) :: rdiag(1:nobsl)
  REAL(r_size),INTENT(IN) :: dep(1:nobsl)
  REAL(r_size),INTENT(OUT) :: dxa_bar(1:nbv-1)
!  REAL(r_size) :: rinv(1:nobsl,1:nobsl)
  REAL(r_size) :: work1(1:nbv-1,1:nobsl)
  REAL(r_size) :: work2(1:nbv-1,1:nobsl)
  INTEGER :: i,j,k
!-----------------------------------------------------------------------
!  Rinv
!-----------------------------------------------------------------------
!-----
! For Future Use (in case that observation error correlation is not zero)
!
!  CALL mtx_inv(nobsl,r,rinv)
!-----
!  rinv = 0.0d0
!  DO i=1,nobsl
!    rinv(i,i) = 1.0d0 / rdiag(i)
!  END DO
!-----------------------------------------------------------------------
!  H^T Rinv
!-----------------------------------------------------------------------
!  work1(1:nrank,1:nobsl) = MATMUL(TRANSPOSE(h_hat(1:nobsl,1:nrank)),rinv(1:nobsl,1:nobsl))
  DO j=1,nobsl
    DO i=1,nrank
!      work1(i,j) = h_hat(1,i)*rinv(1,j)
!      DO k=2,nobsl
!        work1(i,j) = work1(i,j) + h_hat(k,i)*rinv(k,j)
!      END DO
      work1(i,j) = h_hat(j,i) / rdiag(j)
    END DO
  END DO
!-----------------------------------------------------------------------
!  Pa^ H^T Rinv
!-----------------------------------------------------------------------
!  work2 = MATMUL(pa_hat(1:nrank,1:nrank),work1(1:nrank,1:nobsl))
  DO j=1,nobsl
    DO i=1,nrank
      work2(i,j) = pa_hat(i,1)*work1(1,j)
      DO k=2,nrank
        work2(i,j) = work2(i,j) + pa_hat(i,k)*work1(k,j)
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
!  Pa^ H^T Rinv * departure
!-----------------------------------------------------------------------
!  dxa_bar(1:nrank) = MATMUL(work2(1:nrank,1:nobsl),dep(1:nobsl))
  DO i=1,nrank
    dxa_bar(i) = work2(i,1)*dep(1)
    DO k=2,nobsl
      dxa_bar(i) = dxa_bar(i) + work2(i,k)*dep(k)
    END DO
  END DO

  RETURN
END SUBROUTINE calc_dxa_bar
!=======================================================================
!  Compute Y -- Eq.(41)
!    Definition of Y:  dxa^ = dxb^ Y
!    Y = sqrt[I + Xb^T Pb^inv (Pa^ - Pb^) Pb^inv Xb^]
!=======================================================================
SUBROUTINE calc_y(dxb_hat,lambda,pa_hat,y)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: dxb_hat(1:nbv-1,1:nbv)
  REAL(r_size),INTENT(IN) :: lambda(1:nbv)
  REAL(r_size),INTENT(IN) :: pa_hat(1:nbv-1,1:nbv-1)
  REAL(r_size),INTENT(OUT) :: y(1:nbv,1:nbv)
  REAL(r_size) :: work1(1:nbv,1:nbv)
  REAL(r_size) :: work2(1:nbv-1,1:nrank)
  REAL(r_size) :: work3(1:nbv-1,1:nrank)
  REAL(r_size) :: work4(1:nbv-1,1:nrank)
  REAL(r_size) :: work5(1:nbv-1,1:nbv)
  INTEGER :: i,j,k
!-----------------------------------------------------------------------
!  Pa^ - Pb^ -> work2, Pb^inv -> work3
!-----------------------------------------------------------------------
  work2(1:nrank,1:nrank) = pa_hat(1:nrank,1:nrank)
  work3 = 0.0d0
  DO i=1,nrank
    work2(i,i) = pa_hat(i,i) - lambda(i)
    work3(i,i) = 1.0d0 / lambda(i)
  END DO
!-----------------------------------------------------------------------
!  Pb^inv [Pa^ - Pb^] Pb^inv -> work2
!-----------------------------------------------------------------------
!  work4 = MATMUL(work3,work2)
  DO j=1,nrank
    DO i=1,nrank
      work4(i,j) = work3(i,1)*work2(1,j)
      DO k=2,nrank
        work4(i,j) = work4(i,j) + work3(i,k)*work2(k,j)
      END DO
    END DO
  END DO
!  work2 = MATMUL(work4,work3)
  DO j=1,nrank
    DO i=1,nrank
      work2(i,j) = work4(i,1)*work3(1,j)
      DO k=2,nrank
        work2(i,j) = work2(i,j) + work4(i,k)*work3(k,j)
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
!  Xb^T [Pb^inv [Pa^ - Pb^] Pb^inv] Xb^ -> work1
!-----------------------------------------------------------------------
!  work5 = MATMUL(work2(1:nrank,1:nrank),dxb_hat(1:nrank,1:nbv))
  DO j=1,nbv
    DO i=1,nrank
      work5(i,j) = work2(i,1)*dxb_hat(1,j)
      DO k=2,nrank
        work5(i,j) = work5(i,j) + work2(i,k)*dxb_hat(k,j)
      END DO
    END DO
  END DO
!  work1 = MATMUL(TRANSPOSE(dxb_hat),work5)
  DO j=1,nbv
!CDIR NOVECTOR
    DO i=1,nbv
      work1(i,j) = dxb_hat(1,i)*work5(1,j)
      DO k=2,nrank
        work1(i,j) = work1(i,j) + dxb_hat(k,i)*work5(k,j)
      END DO
    END DO
  END DO
  work1 = work1 / REAL(nbv-1,r_size)
!-----------------------------------------------------------------------
!  I + Xb^T [Pb^inv [Pa^ - Pb^] Pb^inv] Xb^
!-----------------------------------------------------------------------
  DO i=1,nbv
    work1(i,i) = work1(i,i) + 1.0d0
  END DO
!-----------------------------------------------------------------------
!  Y = sqrt[I + Xb^T [Pb^inv [Pa^ - Pb^] Pb^inv] Xb^]
!-----------------------------------------------------------------------
  CALL mtx_sqrt(nbv,work1,y)

  RETURN
END SUBROUTINE calc_y
!=======================================================================
!  Ott's Enhanced Variance Inflation Method
!=======================================================================
SUBROUTINE inflate(pa_hat,parm_infl)
  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: pa_hat(1:nrank,1:nrank)
  REAL(r_size),INTENT(IN) :: parm_infl
  REAL(r_size) :: trace
  INTEGER :: i
!-----------------------------------------------------------------------
!  Trace( Pa^ )
!-----------------------------------------------------------------------
  trace = 0.0d0
  DO i=1,nrank
    trace = trace + pa_hat(i,i)
  END DO
!-----------------------------------------------------------------------
!   Pa^ = Pa^ + e Tr(Pa^) / nrank I
!-----------------------------------------------------------------------
  DO i=1,nrank
    pa_hat(i,i) = pa_hat(i,i) + parm_infl * trace / REAL(nrank,r_size)
  END DO

  RETURN
END SUBROUTINE inflate
!=======================================================================
!  Debug tools
!=======================================================================
SUBROUTINE monit_p_hat(p_hat)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: p_hat(1:nrank,1:nrank)
  INTEGER :: i,n

  n = MIN(nrank,8)
  DO i=1,n
    PRINT '(8F9.4)',p_hat(1:n,i)
  END DO

  RETURN
END SUBROUTINE monit_p_hat

SUBROUTINE monit_pb_hat(lambda)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: lambda(1:nbv)
  REAL(r_size) :: p_hat(1:nrank,1:nrank)
  INTEGER :: i,n

  p_hat = 0.0d0
  DO i=1,nrank
    p_hat(i,i) = lambda(i)
  END DO

  n = MIN(nrank,8)
  DO i=1,n
    PRINT '(8F9.4)',p_hat(1:n,i)
  END DO

  RETURN
END SUBROUTINE monit_pb_hat
!=======================================================================
!  E-dimension
!=======================================================================
SUBROUTINE calc_edim_local(dx,edim)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: dx(ndiml,nbv)
  REAL(r_size),INTENT(OUT) :: edim
  REAL(r_size) :: cov(nbv,nbv)
  REAL(r_size) :: sigma(nbv)
  REAL(r_size) :: wk2d(nbv,nbv)
  INTEGER :: i,n1,n2

  DO n1=1,nbv
    DO n2=1,nbv
      cov(n1,n2) = SUM( dx(:,n1) * dx(:,n2) )
    END DO
  END DO

  CALL mtx_eigen(0,nbv,cov,sigma,wk2d,i)

  edim = SUM( SQRT(sigma(:)) )**2 / SUM(sigma(:))

  RETURN
END SUBROUTINE calc_edim_local
!=======================================================================
!  Explained variance and E-dimension
!=======================================================================
SUBROUTINE calc_expvar_local(dxt,dx,expvar,edim)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: dxt(ndiml) ! xf-xt
  REAL(r_size),INTENT(IN) :: dx(ndiml,nbv) ! ensemble perturbation
  REAL(r_size),INTENT(OUT) :: expvar ! explained variance
  REAL(r_size),INTENT(OUT) :: edim ! E-dimension
  REAL(r_size) :: lambda(nbv)
  REAL(r_size) :: q(ndiml,nbv)
  REAL(r_size) :: wk1d(ndiml)
  INTEGER :: n

  CALL calc_pb(dx,lambda,q,n)
  wk1d = MATMUL(q,MATMUL(TRANSPOSE(q),dxt))
  edim = SUM( SQRT(lambda) )**2 / SUM(lambda)
  expvar = SUM( wk1d**2 ) / SUM( dxt**2 )

  RETURN
END SUBROUTINE calc_expvar_local

END MODULE common_lekf
