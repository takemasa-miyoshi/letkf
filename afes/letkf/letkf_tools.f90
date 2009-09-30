MODULE letkf_tools
!=======================================================================
!
! [PURPOSE:] Module for LETKF with AFES
!
! [HISTORY:]
!   01/26/2009 Takemasa Miyoshi  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_obs
  USE common_afes
  USE common_mpi_afes
  USE common_obs_afes
!TVS  USE common_tvs_afes
  USE common_letkf

  IMPLICIT NONE

  PRIVATE
  PUBLIC ::  das_letkf

  INTEGER,SAVE :: nobstotal

!  INTEGER,PARAMETER :: nlev_dampinfl = 15 ! number of levels for inflation damp
  REAL(r_size),PARAMETER :: sp_inflation_n = 0.10d0 !SQRT of cov infl in NH
  REAL(r_size),PARAMETER :: sp_inflation_s = 0.10d0 !SQRT of cov infl in SH
  REAL(r_size),PARAMETER :: sp_infl_additive = 0.d0 !additive inflation
!TVS  LOGICAL,PARAMETER :: msw_vbc = .FALSE.

CONTAINS
!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
SUBROUTINE das_letkf(gues3d,gues2d,anal3d,anal2d)
  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: gues3d(nij1,nlev,nbv,nv3d) ! background ensemble
  REAL(r_size),INTENT(INOUT) :: gues2d(nij1,nbv,nv2d)      !  output: destroyed
  REAL(r_size),INTENT(OUT) :: anal3d(nij1,nlev,nbv,nv3d) ! analysis ensemble
  REAL(r_size),INTENT(OUT) :: anal2d(nij1,nbv,nv2d)
  REAL(r_size),ALLOCATABLE :: mean3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: mean2d(:,:)
  REAL(r_size),ALLOCATABLE :: hdxf(:,:)
  REAL(r_size),ALLOCATABLE :: rdiag(:)
  REAL(r_size),ALLOCATABLE :: dep(:)
  REAL(r_size),ALLOCATABLE :: work3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: work2d(:,:)
  REAL(r_size),ALLOCATABLE :: logpfm(:,:)
!  REAL(r_size),ALLOCATABLE :: parm_infl(:,:) ! spread inflation parameter
  REAL(r_size) :: parm_infl_n,parm_infl_s
  REAL(r_size) :: parm
  REAL(r_size) :: parm_infl_damp(nlev)
  REAL(r_size) :: trans(nbv,nbv)
  INTEGER :: ij,ilev,n,m,i,j,k,nobsl,ierr

  WRITE(6,'(A)') 'Hello from das_letkf'
  nobstotal = nobs !+ ntvs
  WRITE(6,'(A,I8)') 'Target observation numbers : NOBS=',nobs!,', NTVS=',ntvs
  !
  ! In case of no obs
  !
  IF(nobstotal == 0) THEN
    WRITE(6,'(A)') 'No observation assimilated'
!$OMP PARALLEL WORKSHARE
    anal3d = gues3d
    anal2d = gues2d
!$OMP END PARALLEL WORKSHARE
    RETURN
  END IF
  !
  ! INFLATION PARAMETER ESTIMATION
  !
  parm_infl_n = (1.0d0 + sp_inflation_n)**2 - 1.0d0
  parm_infl_s = (1.0d0 + sp_inflation_s)**2 - 1.0d0
  parm_infl_damp(:) = 1.0d0
!  DO ilev=1,nlev
!    parm_infl_damp(ilev) = REAL(nlev-ilev,r_size)/REAL(nlev_dampinfl,r_size)
!    parm_infl_damp(ilev) = MIN(parm_infl_damp(ilev),1.0d0)
!  END DO
  !
  ! FCST PERTURBATIONS
  !
  ALLOCATE(mean3d(nij1,nlev,nv3d))
  ALLOCATE(mean2d(nij1,nv2d))
  CALL ensmean_grd(nbv,nij1,gues3d,gues2d,mean3d,mean2d)
  DO n=1,nv3d
!$OMP PARALLEL DO PRIVATE(i,j,k)
    DO m=1,nbv
      DO k=1,nlev
        DO i=1,nij1
          gues3d(i,k,m,n) = gues3d(i,k,m,n) - mean3d(i,k,n)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END DO
  DO n=1,nv2d
!$OMP PARALLEL DO PRIVATE(i,j,k)
    DO m=1,nbv
      DO i=1,nij1
        gues2d(i,m,n) = gues2d(i,m,n) - mean2d(i,n)
      END DO
    END DO
!$OMP END PARALLEL DO
  END DO
  !
  ! p_full for background ensemble mean
  !
  ALLOCATE(logpfm(nij1,nlev))
  CALL calc_pfull(nij1,1,mean2d(:,iv2d_ps),logpfm)
!$OMP PARALLEL WORKSHARE
  logpfm = DLOG(logpfm)
!$OMP END PARALLEL WORKSHARE
  !
  ! MAIN ASSIMILATION LOOP
  !
!$OMP PARALLEL PRIVATE(ij,ilev,n,i,k,hdxf,rdiag,dep,trans,parm,nobsl)
  ALLOCATE( hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),dep(1:nobstotal) )
!--- For ILEV = 1
! Remark by YS:
!   Incomprehensible error was occured when
!   the ilev loop was iterated from 1 to NLEV.
!---
  ilev=1
!$OMP DO SCHEDULE(DYNAMIC)
    DO ij=1,nij1

      CALL obs_local(ij,ilev,hdxf,rdiag,dep,nobsl,logpfm)

      IF( nobsl /= 0 ) THEN
!!!!! INFLATION SETTING !!!!!
        IF(MAX(sp_inflation_n,sp_inflation_s) == 0.0) THEN
          parm = 0.0d0
        ELSE
          IF(sp_inflation_n == sp_inflation_s) THEN
            parm = parm_infl_n
          ELSE
            IF(lat1(ij) > 20.d0) THEN
              parm = parm_infl_n
            ELSE IF(lat1(ij) < -20.d0) THEN
              parm = parm_infl_s
            ELSE
              parm = (lat1(ij) + 20.d0) / 40.d0
              parm = parm * parm_infl_n + (1.d0 - parm) * parm_infl_s
            END IF
          END IF
        END IF
        CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,dep,parm,trans)

        DO n=1,nv3d
          DO m=1,nbv
            anal3d(ij,ilev,m,n) = mean3d(ij,ilev,n)
            DO k=1,nbv
              anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
                & + gues3d(ij,ilev,k,n) * trans(k,m)
            END DO
          END DO
        END DO
        DO n=1,nv2d
          DO m=1,nbv
            anal2d(ij,m,n)  = mean2d(ij,n)
            DO k=1,nbv
              anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,k,n) * trans(k,m)
            END DO
          END DO
        END DO

      ELSE ! no observation for this analysis
        DO n=1,nv3d
          DO m=1,nbv
            anal3d(ij,ilev,m,n)  = mean3d(ij,ilev,n) + gues3d(ij,ilev,m,n)
          END DO
        END DO
        DO n=1,nv2d
          DO m=1,nbv
            anal2d(ij,m,n)  = mean2d(ij,n) + gues2d(ij,m,n)
          END DO
        END DO
      END IF

    END DO
!$OMP END DO
!--- For ILEV = 2 - NLEV
!$OMP DO SCHEDULE(DYNAMIC)
  DO ilev=2,nlev
    DO ij=1,nij1

      CALL obs_local(ij,ilev,hdxf,rdiag,dep,nobsl,logpfm)

      IF( nobsl /= 0 ) THEN
!!!!! INFLATION SETTING !!!!!
        IF(MAX(sp_inflation_n,sp_inflation_s) == 0.0) THEN
          parm = 0.0d0
        ELSE
          IF(sp_inflation_n == sp_inflation_s) THEN
            parm = parm_infl_n
          ELSE
            IF(lat1(ij) > 20.d0) THEN
              parm = parm_infl_n
            ELSE IF(lat1(ij) < -20.d0) THEN
              parm = parm_infl_s
            ELSE
              parm = (lat1(ij) + 20.d0) / 40.d0
              parm = parm * parm_infl_n + (1.d0 - parm) * parm_infl_s
            END IF
          END IF
          parm = parm * parm_infl_damp(ilev)
        END IF
        CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,dep,parm,trans)

        DO n=1,nv3d
          DO m=1,nbv
            anal3d(ij,ilev,m,n) = mean3d(ij,ilev,n)
            DO k=1,nbv
              anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
                & + gues3d(ij,ilev,k,n) * trans(k,m)
            END DO
          END DO
        END DO

      ELSE ! no observation for this analysis
        DO n=1,nv3d
          DO m=1,nbv
            anal3d(ij,ilev,m,n)  = mean3d(ij,ilev,n) + gues3d(ij,ilev,m,n)
          END DO
        END DO
      END IF

    END DO
  END DO
!$OMP END DO
  DEALLOCATE(hdxf,rdiag,dep)
!$OMP END PARALLEL
  !
  ! Additive inflation
  !
  IF(sp_infl_additive > 0.0d0) THEN
    CALL read_ens_mpi('addi',nbv,gues3d,gues2d)
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    CALL ensmean_grd(nbv,nij1,gues3d,gues2d,work3d,work2d)
    DO n=1,nv3d
!$OMP PARALLEL DO PRIVATE(i,j,k)
      DO m=1,nbv
        DO k=1,nlev
          DO i=1,nij1
            gues3d(i,k,m,n) = gues3d(i,k,m,n) - work3d(i,k,n)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    END DO
    DO n=1,nv2d
!$OMP PARALLEL DO PRIVATE(i,j,k)
      DO m=1,nbv
        DO i=1,nij1
          gues2d(i,m,n) = gues2d(i,m,n) - work2d(i,n)
        END DO
      END DO
!$OMP END PARALLEL DO
    END DO

    DEALLOCATE(work3d,work2d)
    WRITE(6,'(A)') '===== Additive covariance inflation ====='
    WRITE(6,'(A,F10.4)') '  parameter:',sp_infl_additive
    WRITE(6,'(A)') '========================================='
!    parm = 0.7d0
!    DO ilev=1,nlev
!      parm_infl_damp(ilev) = 1.0d0 + parm &
!        & + parm * REAL(1-ilev,r_size)/REAL(nlev_dampinfl,r_size)
!      parm_infl_damp(ilev) = MAX(parm_infl_damp(ilev),1.0d0)
!    END DO
    DO n=1,nv3d
      DO m=1,nbv
        DO ilev=1,nlev
          DO ij=1,nij1
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
              & + gues3d(ij,ilev,m,n) * sp_infl_additive
          END DO
        END DO
      END DO
    END DO
    DO n=1,nv2d
      DO m=1,nbv
        DO ij=1,nij1
          anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,m,n) * sp_infl_additive
        END DO
      END DO
    END DO
  END IF
  !
  ! Force non-negative q
  !
!$OMP PARALLEL DO PRIVATE(ilev,ij)
  DO m=1,nbv
    DO ilev=1,nlev
      DO ij=1,nij1
        anal3d(ij,ilev,m,iv3d_q) = MAX(anal3d(ij,ilev,m,iv3d_q),0.0d0)
        anal3d(ij,ilev,m,iv3d_cw) = MAX(anal3d(ij,ilev,m,iv3d_cw),0.0d0)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  DEALLOCATE(logpfm,mean3d,mean2d)
  RETURN
END SUBROUTINE das_letkf
!-----------------------------------------------------------------------
! Project global observations to local
!     (hdxf_g,dep_g,rdiag_g) -> (hdxf,dep,rdiag)
!-----------------------------------------------------------------------
SUBROUTINE obs_local(ij,ilev,hdxf,rdiag,dep,nobsl,logpfm)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ij,ilev
  REAL(r_size),INTENT(IN) :: logpfm(nij1,nlev)
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  INTEGER,INTENT(OUT) :: nobsl
  REAL(r_size) :: minlon,maxlon,minlat,maxlat,dist,dlev
  REAL(r_size) :: tmplon,tmplat,tmperr,tmpwgt(nlev)
  INTEGER :: tmpqc
  INTEGER,ALLOCATABLE:: nobs_use(:)
!TVS  INTEGER,ALLOCATABLE:: ntvs_use_prof(:),ntvs_use_inst(:),ntvs_use_slot(:)
  INTEGER :: imin,imax,jmin,jmax,im,ichan
  INTEGER :: n,nn,tvnn
!
! INITIALIZE
!
  IF( nobs > 0 ) THEN
    ALLOCATE(nobs_use(nobs))
  END IF
!TVS  IF( ntvs > 0 ) THEN
!TVS    ALLOCATE(ntvs_use_prof(ntvs))
!TVS    ALLOCATE(ntvs_use_inst(ntvs))
!TVS    ALLOCATE(ntvs_use_slot(ntvs))
!TVS  END IF
!
! data search
!
  minlon = lon1(ij) - dlon_zero(ij)
  maxlon = lon1(ij) + dlon_zero(ij)
  minlat = lat1(ij) - dlat_zero
  maxlat = lat1(ij) + dlat_zero

  DO jmin=1,nlat-2
    IF(minlat < lat(jmin+1)) EXIT
  END DO
  DO jmax=1,nlat-2
    IF(maxlat < lat(jmax+1)) EXIT
  END DO
  nn = 1
!TVS  tvnn = 1
  IF(minlon >= 0 .AND. maxlon <= 360.0) THEN
    DO imin=1,nlon-1
      IF(minlon < lon(imin+1)) EXIT
    END DO
    DO imax=1,nlon-1
      IF(maxlon < lon(imax+1)) EXIT
    END DO
    IF( nobs > 0 ) &
    & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS    IF( ntvs > 0 ) &
!TVS    & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS    &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
  ELSE IF(minlon >= 0 .AND. maxlon > 360.0) THEN
    DO imin=1,nlon-1
      IF(minlon < lon(imin+1)) EXIT
    END DO
    maxlon = maxlon - 360.0d0
    IF(maxlon > 360.0d0) THEN
      imin = 1
      imax = nlon
      IF( nobs > 0 ) &
      & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS      IF( ntvs > 0 ) &
!TVS      & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS      &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
    ELSE
      DO imax=1,nlon-1
        IF(maxlon < lon(imax+1)) EXIT
      END DO
      IF(imax > imin) THEN
        imin = 1
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      ELSE
        imin = 1
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
        DO imin=1,nlon-1
          IF(minlon < lon(imin+1)) EXIT
        END DO
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      END IF
    END IF
  ELSE IF(minlon < 0 .AND. maxlon <= 360.0d0) THEN
    DO imax=1,nlon-1
      IF(maxlon < lon(imax+1)) EXIT
    END DO
    minlon = minlon + 360.0d0
    IF(minlon < 0) THEN
      imin = 1
      imax = nlon
      IF( nobs > 0 ) &
      & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS      IF( ntvs > 0 ) &
!TVS      & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS      &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
    ELSE
      DO imin=1,nlon-1
        IF(minlon < lon(imin+1)) EXIT
      END DO
      IF(imin < imax) THEN
        imin = 1
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      ELSE
        imin = 1
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
        DO imin=1,nlon-1
          IF(minlon < lon(imin+1)) EXIT
        END DO
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      END IF
    END IF
  ELSE
    maxlon = maxlon - 360.0d0
    minlon = minlon + 360.0d0
    IF(maxlon > 360.0 .OR. minlon < 0) THEN
      imin = 1
      imax = nlon
      IF( nobs > 0 ) &
      & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS      IF( ntvs > 0 ) &
!TVS      & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS      &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
    ELSE
      DO imin=1,nlon-1
        IF(minlon < lon(imin+1)) EXIT
      END DO
      DO imax=1,nlon-1
        IF(maxlon < lon(imax+1)) EXIT
      END DO
      IF(imin > imax) THEN
        imin = 1
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      ELSE
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      END IF
    END IF
  END IF
  nn = nn-1
!TVS  tvnn = tvnn -1
!TVS  IF( nn < 1 .AND. tvnn < 1 ) THEN
  IF(nn < 1) THEN
    nobsl = 0
    RETURN
  END IF
!
! CONVENTIONAL
!
  nobsl = 0
  IF(nn > 0) THEN
    DO n=1,nn
      IF(NINT(obselm(nobs_use(n))) == id_ps_obs .AND. ilev > 1) THEN
        dlev = ABS(LOG(obsdat(nobs_use(n))) - logpfm(ij,ilev))
      ELSE IF(NINT(obselm(nobs_use(n))) /= id_ps_obs) THEN
        dlev = ABS(LOG(obslev(nobs_use(n))) - logpfm(ij,ilev))
      ELSE
        dlev = 0.0d0
      END IF
      IF(dlev > dist_zerov) CYCLE

      tmplon=obslon(nobs_use(n))
      tmplat=obslat(nobs_use(n))
      CALL com_distll_1( tmplon, tmplat,lon1(ij), lat1(ij), dist)
      IF(dist > dist_zero ) CYCLE

      nobsl = nobsl + 1
      hdxf(nobsl,:) = obshdxf(nobs_use(n),:)
      dep(nobsl)    = obsdep(nobs_use(n))
      !
      ! Observational localization
      !
      tmperr=obserr(nobs_use(n))
      rdiag(nobsl) = tmperr * tmperr &
        & * exp(0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2))
    END DO
  END IF
!TVS!
!TVS! ATOVS
!TVS!
!TVS  IF(tvnn > 0) THEN
!TVS    DO n=1,tvnn
!TVS      tmplon=tvslon(ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n))
!TVS      tmplat=tvslat(ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n))
!TVS      CALL com_distll_1( tmplon, tmplat, lon1(ij), lat1(ij), dist)
!TVS      IF( dist > dist_zero) CYCLE
!TVS
!TVS      DO ichan=1,ntvsch(ntvs_use_inst(n))
!TVS        tmperr=tvserr(ichan,ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n))
!TVS        tmpqc=tvsqc(ichan,ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n))
!TVS        tmpwgt(:)=tvswgt(:,ichan, &
!TVS                         & ntvs_use_prof(n), &
!TVS                         & ntvs_use_inst(n), &
!TVS                         & ntvs_use_slot(n))
!TVS        IF( tmpqc == 1 .AND. tmpwgt(ilev) > 0.05D0 ) THEN
!TVS          nobsl = nobsl + 1
!TVS          DO im = 1, nbv
!TVS            hdxf(nobsl,im) = tvshdxf(im,ichan, &
!TVS                              & ntvs_use_prof(n), &
!TVS                              & ntvs_use_inst(n), &
!TVS                              & ntvs_use_slot(n))
!TVS          END DO
!TVS          dep(nobsl)    = tvsdep(ichan, &
!TVS                              & ntvs_use_prof(n), &
!TVS                              & ntvs_use_inst(n), &
!TVS                              & ntvs_use_slot(n))
!TVS          rdiag(nobsl)  = tmperr * tmperr &
!TVS                        & * exp(0.5d0 * (dist/sigma_obs)**2) &
!TVS                        & / (tmpwgt(ilev) * tmpwgt(ilev))
!TVS        END IF
!TVS      END DO
!TVS    END DO
!TVS  END IF
!
! DEBUG
! IF( ILEV == 1 .AND. ILON == 1 ) &
! & WRITE(6,*) 'ILEV,ILON,ILAT,NN,TVNN,NOBSL=',ilev,ij,nn,tvnn,nobsl
!
  IF( nobsl > nobstotal ) THEN
    WRITE(6,'(A,I5,A,I5)') 'FATAL ERROR, NOBSL=',nobsl,' > NOBSTOTAL=',nobstotal
    WRITE(6,*) 'IJ,NN,TVNN=', ij, nn, tvnn
    STOP 99
  END IF
!
  IF( nobs > 0 ) THEN
    DEALLOCATE(nobs_use)
  END IF
!TVS  IF( ntvs > 0 ) THEN
!TVS    DEALLOCATE(ntvs_use_prof)
!TVS    DEALLOCATE(ntvs_use_inst)
!TVS    DEALLOCATE(ntvs_use_slot)
!TVS  END IF
!
  RETURN
END SUBROUTINE obs_local

SUBROUTINE obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax
  INTEGER,INTENT(INOUT) :: nn, nobs_use(nobs)
  INTEGER :: j,n,ib,ie,ip

  DO j=jmin,jmax
    IF(imin > 1) THEN
      ib = nobsgrd(imin-1,j)+1
    ELSE
      IF(j > 1) THEN
        ib = nobsgrd(nlon,j-1)+1
      ELSE
        ib = 1
      END IF
    END IF
    ie = nobsgrd(imax,j)
    n = ie - ib + 1
    IF(n == 0) CYCLE
    DO ip=ib,ie
      IF(nn > nobs) THEN
        WRITE(6,*) 'FATALERROR, NN > NOBS', NN, NOBS
      END IF
      nobs_use(nn) = ip
      nn = nn + 1
    END DO
  END DO

  RETURN
END SUBROUTINE obs_local_sub

!TVSSUBROUTINE tvs_local_sub(imin,imax,jmin,jmax,nn,ntvs_prof,ntvs_inst,ntvs_slot)
!TVS  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax
!TVS  INTEGER,INTENT(INOUT) :: nn, ntvs_prof(ntvs), ntvs_inst(ntvs), ntvs_slot(ntvs)
!TVS  INTEGER :: j,n,ib,ie,ip
!TVS  INTEGER :: islot, iinst
!TVS
!TVS  DO j=jmin,jmax
!TVS    DO islot=1,nslots
!TVS      DO iinst=1,ninstrument
!TVS        IF(imin > 1) THEN
!TVS          ib = ntvsgrd(imin-1,j,iinst,islot)+1
!TVS        ELSE
!TVS          IF(j > 1) THEN
!TVS            ib = ntvsgrd(nlon,j-1,iinst,islot)+1
!TVS          ELSE
!TVS            ib = 1
!TVS          END IF
!TVS        END IF
!TVS        ie = ntvsgrd(imax,j,iinst,islot)
!TVS        n = ie - ib + 1
!TVS        IF(n == 0) CYCLE
!TVS        DO ip=ib,ie
!TVS          IF(nn > nobs) THEN
!TVS            WRITE(6,*) 'FATALERROR, NN > NTVS', NN, NTVS
!TVS          END IF
!TVS          ntvs_prof(nn)=ip
!TVS          ntvs_inst(nn)=iinst
!TVS          ntvs_slot(nn)=islot
!TVS          nn = nn + 1
!TVS        END DO
!TVS      END DO
!TVS    END DO
!TVS  END DO
!TVS  RETURN
!TVSEND SUBROUTINE tvs_local_sub
!TVS!-----------------------------------------------------------------------
!TVS! Data Assimilation for VARBC
!TVS!-----------------------------------------------------------------------
!TVSSUBROUTINE das_vbc(um,vm,tm,qm,qlm,psm,vbcf,vbca)
!TVS  USE common_mtx
!TVS  IMPLICIT NONE
!TVS  REAL(r_size),INTENT(IN) :: um(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: vm(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: tm(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: qm(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: qlm(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: psm(nij1)
!TVS  REAL(r_size),INTENT(INOUT) :: vbcf(maxvbc,maxtvsch,ninstrument)
!TVS  REAL(r_size),INTENT(OUT)   :: vbca(maxvbc,maxtvsch,ninstrument)
!TVS  REAL(r_sngl) :: u4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: v4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: t4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: q4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: ql4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: ps4(nlon,nlat)
!TVS  REAL(r_size) :: u(nlon,nlat,nlev)
!TVS  REAL(r_size) :: v(nlon,nlat,nlev)
!TVS  REAL(r_size) :: t(nlon,nlat,nlev)
!TVS  REAL(r_size) :: q(nlon,nlat,nlev)
!TVS  REAL(r_size) :: ql(nlon,nlat,nlev)
!TVS  REAL(r_size) :: ps(nlon,nlat)
!TVS  REAL(r_size) :: p_full(nlon,nlat,nlev)
!TVS  REAL(r_size),ALLOCATABLE :: hx(:,:,:,:)
!TVS  REAL(r_size),ALLOCATABLE :: pred(:,:,:,:,:)
!TVS  INTEGER,ALLOCATABLE :: tmpqc(:,:,:)
!TVS  REAL(r_size),ALLOCATABLE :: tmpwgt(:,:,:,:)
!TVS  REAL(r_size) :: a(maxvbc,maxvbc)
!TVS  REAL(r_size) :: b(maxvbc)
!TVS  REAL(r_size) :: ainv(maxvbc,maxvbc)
!TVS  INTEGER:: ntvschan1(maxtvsch,ninstrument)
!TVS  INTEGER:: i,j,k,n,islot,nn
!TVS
!TVS  PRINT *,'Hello from das_vbc'
!TVS
!TVS  IF(ntvs == 0) THEN
!TVS    PRINT *,'No radiance data: das_vbc skipped..'
!TVS!$OMP PARALLEL WORKSHARE
!TVS    vbca = vbcf
!TVS!$OMP END PARALLEL WORKSHARE
!TVS    RETURN
!TVS  END IF
!TVS
!TVS  CALL gather_grd_mpi(0,um,vm,tm,qm,qlm,psm,u4,v4,t4,q4,ql4,ps4)
!TVS  n = nlon*nlat*nlev
!TVS  CALL MPI_BARRIER(MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(u4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(v4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(t4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(q4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(ql4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  n = nlon*nlat
!TVS  CALL MPI_BCAST(ps4(1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS!$OMP PARALLEL WORKSHARE
!TVS  u = REAL(u4,r_size)
!TVS  v = REAL(v4,r_size)
!TVS  t = REAL(t4,r_size)
!TVS  q = REAL(q4,r_size)
!TVS  ql = REAL(ql4,r_size)
!TVS  ps = REAL(ps4,r_size)
!TVS!$OMP END PARALLEL WORKSHARE
!TVS  CALL calc_pfull(ps,p_full)
!TVS
!TVS  ALLOCATE( hx(maxtvsch,maxtvsprof,ninstrument,nslots) )
!TVS  ALLOCATE( pred(maxvbc,maxtvsch,maxtvsprof,ninstrument,nslots) )
!TVS  ALLOCATE( tmpqc(maxtvsch,maxtvsprof,ninstrument) )
!TVS  ALLOCATE( tmpwgt(nlev,maxtvsch,maxtvsprof,ninstrument) )
!TVS  DO islot=1,nslots
!TVS!    IF(SUM(ntvsprofslots(:,islot)) == 0) CYCLE
!TVS    ntvsprof(:) = ntvsprofslots(:,islot)
!TVS    CALL Trans_XtoY_tvs(u,v,t,q,ql,ps,p_full, &
!TVS      & tvslon(:,:,islot),tvslat(:,:,islot),tvszenith(:,:,islot),&
!TVS      & tvsskin(:,:,islot),tvsstmp(:,:,islot),tvsclw(:,:,islot),&
!TVS      & tvsemis(:,:,:,islot),tmpqc,hx(:,:,:,islot),tmpwgt,pred(:,:,:,:,islot))
!TVS  END DO
!TVS  DEALLOCATE(tmpqc,tmpwgt)
!TVS
!TVS!$OMP PARALLEL PRIVATE(j,k,n,a,b,ainv)
!TVS!$OMP WORKSHARE
!TVS  vbca = 0.0d0
!TVS!$OMP END WORKSHARE
!TVS!$OMP DO SCHEDULE(DYNAMIC)
!TVS  DO k=1,ninstrument
!TVS    DO j=1,maxtvsch
!TVS      !
!TVS      ! Parallel processing
!TVS      !
!TVS      IF(MOD(j+maxtvsch*(k-1)-1,nprocs) /= myrank) CYCLE
!TVS      !
!TVS      ! DATA NUMBER
!TVS      !
!TVS      ntvschan(j,k) = SUM(tvsqc(j,:,k,:))
!TVS      IF(msw_vbc .AND. ntvschan(j,k) /= 0 ) THEN
!TVS        PRINT '(3A,I3,A,I6)',' >> VBC executed for instrument,channel,ntvsl: ',&
!TVS                            & tvsname(k),',',tvsch(j,k),',',ntvschan(j,k)
!TVS        CALL vbc_local(j,k,ntvschan(j,k),hx,pred,a,b)
!TVS        CALL mtx_inv(maxvbc,a,ainv)
!TVS        vbca(:,j,k) = vbcf(:,j,k)
!TVS        DO n=1,maxvbc
!TVS          vbca(:,j,k) = vbca(:,j,k) - ainv(:,n)*b(n) !ATTN: sign for beta
!TVS        END DO
!TVS      ELSE
!TVS        PRINT '(3A,I3,A,I6)',' !! NO VBC executed for instrument,channel,ntvsl: ',&
!TVS                            & tvsname(k),',',tvsch(j,k),',',ntvschan(j,k)
!TVS        vbca(:,j,k) = vbcf(:,j,k)
!TVS      END IF
!TVS    END DO
!TVS  END DO
!TVS!$OMP END DO
!TVS!$OMP WORKSHARE
!TVS  vbcf = vbca
!TVS  ntvschan1 = ntvschan
!TVS!$OMP END WORKSHARE
!TVS!$OMP END PARALLEL
!TVS  DEALLOCATE(hx,pred)
!TVS  n = maxvbc*maxtvsch*ninstrument
!TVS  CALL MPI_BARRIER(MPI_COMM_WORLD,j)
!TVS  CALL MPI_ALLREDUCE(vbcf,vbca,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,j)
!TVS  n = maxtvsch*ninstrument
!TVS  CALL MPI_BARRIER(MPI_COMM_WORLD,j)
!TVS  CALL MPI_ALLREDUCE(ntvschan1,ntvschan,n,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,j)
!TVS
!TVS  RETURN
!TVSEND SUBROUTINE das_vbc
!TVS!-----------------------------------------------------------------------
!TVS!  (in) ichan: channnel
!TVS!  (in) iinst: sensor
!TVS!  (out) a = B_beta^-1 + p R^-1 p^T
!TVS!  (out) b = p R^-1 d
!TVS!-----------------------------------------------------------------------
!TVSSUBROUTINE vbc_local(ichan,iinst,ntvsl,hx,pred,a,b)
!TVS  IMPLICIT NONE
!TVS  INTEGER,PARAMETER :: msw=1
!TVS  INTEGER,PARAMETER :: nmin=400
!TVS  INTEGER,INTENT(IN) :: ichan,iinst,ntvsl
!TVS  REAL(r_size),INTENT(IN) :: hx(maxtvsch,maxtvsprof,ninstrument,nslots)
!TVS  REAL(r_size),INTENT(IN) :: pred(maxvbc,maxtvsch,maxtvsprof,ninstrument,nslots)
!TVS  REAL(r_size),INTENT(OUT) :: a(maxvbc,maxvbc)
!TVS  REAL(r_size),INTENT(OUT) :: b(maxvbc)
!TVS  REAL(r_size) :: dep,dep0
!TVS  REAL(r_size) :: bias,bias0
!TVS  REAL(r_size) :: r,tmp
!TVS  INTEGER:: islot, iprof, i,j,n
!TVS
!TVS  a = 0.0d0
!TVS  b = 0.0d0
!TVS  dep = 0.0d0
!TVS  dep0 = 0.0d0
!TVS  bias = 0.0d0
!TVS  bias0 = 0.0d0
!TVS  n = 0
!TVS  DO islot=1,nslots
!TVS    DO iprof=1,maxtvsprof
!TVS      IF(tvsqc(ichan,iprof,iinst,islot)/=1) CYCLE
!TVS      !
!TVS      ! R
!TVS      !
!TVS      r = tvserr(ichan,iprof,iinst,islot)**2
!TVS      !
!TVS      ! p R^-1 p^T
!TVS      !
!TVS      DO j=1,maxvbc
!TVS        DO i=1,maxvbc
!TVS          a(i,j) = a(i,j) &
!TVS               & + pred(i,ichan,iprof,iinst,islot) &
!TVS               & * pred(j,ichan,iprof,iinst,islot) / r
!TVS        END DO
!TVS      END DO
!TVS      !
!TVS      ! B_beta^-1
!TVS      !
!TVS      IF(msw == 1) THEN ! Y.Sato
!TVS        IF(ntvsl < nmin) THEN
!TVS          tmp = REAL(nmin,r_size) / r
!TVS
!TVS        ELSE
!TVS          tmp = (REAL(ntvsl,r_size) &
!TVS            & / (LOG10(REAL(ntvsl,r_size)/REAL(nmin,r_size))+1.0d0)) / r
!TVS        END IF
!TVS      ELSE IF(msw == 2) THEN ! D.Dee
!TVS        tmp = REAL(ntvsl,r_size) / r
!TVS      ELSE ! Constant
!TVS        tmp = 100.0d0
!TVS      END IF
!TVS      DO i=1,maxvbc
!TVS        a(i,i) = a(i,i) + tmp
!TVS      END DO
!TVS      !
!TVS      ! p R^-1 d
!TVS      !
!TVS      b(:) = b(:) + pred(:,ichan,iprof,iinst,islot) / r &
!TVS                & *(tvsdat(ichan,iprof,iinst,islot)-hx(ichan,iprof,iinst,islot))
!TVS      bias = bias+tvsdat(ichan,iprof,iinst,islot)-hx(ichan,iprof,iinst,islot)
!TVS      dep = dep+(tvsdat(ichan,iprof,iinst,islot)-hx(ichan,iprof,iinst,islot))**2
!TVS      bias0= bias0+tvsdep(ichan,iprof,iinst,islot)
!TVS      dep0= dep0+tvsdep(ichan,iprof,iinst,islot)**2
!TVS      n = n+1
!TVS    END DO
!TVS  END DO
!TVS
!TVS  dep = SQRT(dep / REAL(n,r_size))
!TVS  dep0 = SQRT(dep0 / REAL(n,r_size))
!TVS  bias = bias / REAL(n,r_size)
!TVS  bias0 = bias0 / REAL(n,r_size)
!TVS  PRINT '(2A,I3,4F12.4)',' >> D monit: ',tvsname(iinst),tvsch(ichan,iinst),bias0,bias,dep0,dep
!TVS
!TVS  RETURN
!TVSEND SUBROUTINE vbc_local

END MODULE letkf_tools
