MODULE letkf_tools
!=======================================================================
!
! [PURPOSE:] Module for LETKF with ROMS
!
! [HISTORY:]
!   01/26/2009 Takemasa Miyoshi  created
!   02/03/2009 Takemasa Miyoshi  modified for ROMS
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_obs
  USE common_roms
  USE common_mpi_roms
  USE common_obs_roms
  USE common_letkf

  IMPLICIT NONE

  PRIVATE
  PUBLIC ::  das_letkf

  INTEGER,SAVE :: nobstotal

!  INTEGER,PARAMETER :: nlev_dampinfl = 15 ! number of levels for inflation damp
  REAL(r_size),PARAMETER :: sp_inflation = 0.02d0 !SQRT of cov infl
  REAL(r_size),PARAMETER :: sp_infl_additive = 0.d0 !additive inflation

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
  REAL(r_size) :: parm_infl
  REAL(r_size) :: parm,parm_max,parm_min
  REAL(r_size) :: parm_infl_damp(nlev)
  REAL(r_size) :: pu,pd
  REAL(r_size) :: trans(nbv,nbv)
  INTEGER :: ij,ilev,n,m,i,j,k,nobsl,ierr

  PRINT *,'Hello from das_letkf'
  nobstotal = nobs
  PRINT '(A,I8)','Target observation numbers : NOBS=',nobs
  !
  ! In case of no obs
  !
  IF(nobstotal == 0) THEN
    PRINT '(A)','No observation assimilated'
!$OMP PARALLEL WORKSHARE
    anal3d = gues3d
    anal2d = gues2d
!$OMP END PARALLEL WORKSHARE
    RETURN
  END IF
  !
  ! INFLATION PARAMETER ESTIMATION
  !
  parm_infl = (1.0d0 + sp_inflation)**2 - 1.0d0
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
  ! MAIN ASSIMILATION LOOP
  !
!$OMP PARALLEL PRIVATE(ij,ilev,n,i,k,hdxf,rdiag,dep,trans,parm,nobsl)
  ALLOCATE( hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),dep(1:nobstotal) )
!--- For ILEV = nlev (surface)
  ilev=nlev
!$OMP DO SCHEDULE(DYNAMIC)
    DO ij=1,nij1

      CALL obs_local(ij,ilev,hdxf,rdiag,dep,nobsl)

      IF( nobsl /= 0 ) THEN
!!!!! INFLATION SETTING !!!!!
        parm = parm_infl
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
          IF(n==iv2d_ubar .OR. n==iv2d_vbar .OR. n==iv2d_hbl) THEN
            DO m=1,nbv
              anal2d(ij,m,n)  = mean2d(ij,n) + gues2d(ij,m,n)
            END DO
          ELSE
            DO m=1,nbv
              anal2d(ij,m,n)  = mean2d(ij,n)
              DO k=1,nbv
                anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,k,n) * trans(k,m)
              END DO
            END DO
          END IF
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
!--- For ILEV = 1 - (NLEV-1)
!$OMP DO SCHEDULE(DYNAMIC)
  DO ilev=1,nlev-1
    DO ij=1,nij1

      CALL obs_local(ij,ilev,hdxf,rdiag,dep,nobsl)

      IF( nobsl /= 0 ) THEN
!!!!! INFLATION SETTING !!!!!
        parm = parm_infl * parm_infl_damp(ilev)
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
    PRINT '(A)','===== Additive covariance inflation ====='
    PRINT '(A,F10.4)','  parameter:',sp_infl_additive
    PRINT '(A)','========================================='
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
        anal3d(ij,ilev,m,iv3d_s) = MAX(anal3d(ij,ilev,m,iv3d_s),0.0d0)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  DEALLOCATE(mean3d,mean2d)
  RETURN
END SUBROUTINE das_letkf
!-----------------------------------------------------------------------
! Project global observations to local
!     (hdxf_g,dep_g,rdiag_g) -> (hdxf,dep,rdiag)
!-----------------------------------------------------------------------
SUBROUTINE obs_local(ij,ilev,hdxf,rdiag,dep,nobsl)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ij,ilev
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  INTEGER,INTENT(OUT) :: nobsl
  REAL(r_size) :: minlon,maxlon,minlat,maxlat,dist,dlev
  INTEGER,ALLOCATABLE:: nobs_use(:)
  INTEGER :: imin,imax,jmin,jmax,im,ichan
  INTEGER :: n,nn
!
! INITIALIZE
!
  IF( nobs > 0 ) THEN
    ALLOCATE(nobs_use(nobs))
  END IF
!
! data search
!
  imin = MAX(NINT(ri1(ij) - dist_zero),1)
  imax = MIN(NINT(ri1(ij) + dist_zero),nlon)
  jmin = MAX(NINT(rj1(ij) - dist_zero),1)
  jmax = MIN(NINT(rj1(ij) + dist_zero),nlat)

  nn = 1
  IF( nobs > 0 ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
  nn = nn-1
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
      IF(NINT(obselm(nobs_use(n))) == id_z_obs .AND. ilev < nlev) THEN
        dlev = ABS(nlev - ilev)
      ELSE IF(NINT(obselm(nobs_use(n))) /= id_z_obs) THEN
        dlev = ABS(NINT(obsk(nobs_use(n))) - ilev)
      ELSE
        dlev = 0.0d0
      END IF
      IF(dlev > dist_zerov) CYCLE

      dist = SQRT((obsi(nobs_use(n))-ri1(ij))**2 &
        & + (obsj(nobs_use(n))-rj1(ij))**2)
      IF(dist > dist_zero ) CYCLE

      nobsl = nobsl + 1
      hdxf(nobsl,:) = obshdxf(nobs_use(n),:)
      dep(nobsl)    = obsdep(nobs_use(n))
      !
      ! Observational localization
      !
      rdiag(nobsl) = obserr(nobs_use(n))**2 &
        & * exp(0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2))
    END DO
  END IF
!
! DEBUG
! IF( ILEV == 1 .AND. ILON == 1 ) &
! & WRITE(6,*) 'ILEV,ILON,ILAT,NN,TVNN,NOBSL=',ilev,ij,nn,tvnn,nobsl
!
  IF( nobsl > nobstotal ) THEN
    WRITE(6,'(A,I5,A,I5)') 'FATAL ERROR, NOBSL=',nobsl,' > NOBSTOTAL=',nobstotal
    WRITE(6,*) 'IJ,NN=',ij,nn
    STOP 99
  END IF
!
  IF( nobs > 0 ) THEN
    DEALLOCATE(nobs_use)
  END IF
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

END MODULE letkf_tools
