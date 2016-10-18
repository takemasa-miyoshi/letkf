MODULE efso_tools
!=======================================================================
!
! [PURPOSE:] Module for observation sensitivity calculation
!
! [HISTORY:]
!   07/27/2011 Yoichiro Ohta  created
!   09/29/2011 Yoichiro Ohta  adapted new formulation
!   07/01/2013 Daisuke Hotta  ported to GFS-LETKF system
!   12/19/2013 Guo-Yuan Lien  merged to GFS-LETKF main development
!   01/02/2013 Guo-Yuan Lien  modify output format
!
!=======================================================================
  USE common
  USE common_gfs
  USE common_obs_gfs
  USE common_mpi
  USE common_mpi_gfs
  USE letkf_obs
  USE efso_nml
  USE sigio_module

  IMPLICIT NONE
  PRIVATE
  PUBLIC init_obsense,destroy_obsense,print_obsense,lnorm,loc_advection
  PUBLIC obsense,lon2,lat2,nterm

  REAL(r_size),ALLOCATABLE :: obsense(:,:)
  REAL(r_size),ALLOCATABLE :: lon2(:,:)
  REAL(r_size),ALLOCATABLE :: lat2(:,:)
  INTEGER,PARAMETER :: nterm = 3

CONTAINS

SUBROUTINE init_obsense
  IMPLICIT NONE
  ALLOCATE(obsense(nterm,nobs))
  RETURN
END SUBROUTINE init_obsense

!-----------------------------------------------------------------------
! Compute norm
! [ref: Eq.(6,7,9), Ota et al. 2013]
!-----------------------------------------------------------------------
! [INPUT]
!  fcst3d,fcst2d: (xmean+X)^f_t  (total field)
!  fcer3d,fcer2d: [1/2(K-1)](e^f_t+e^g_t)
! [OUTPUT]
!  fcst3d,fcst2d: C^(1/2)*X^f_t                    [(J/kg)^(1/2)]
!  fcer3d,fcer2d: C^(1/2)*[1/2(K-1)](e^f_t+e^g_t)  [(J/kg)^(1/2)]
!-----------------------------------------------------------------------
SUBROUTINE lnorm(fcst3d,fcst2d,fcer3d,fcer2d)
  IMPLICIT NONE
  REAL(r_size), INTENT(INOUT) :: fcst3d(nij1,nlev,nbv,nv3d)
  REAL(r_size), INTENT(INOUT) :: fcst2d(nij1,nbv,nv2d)
  REAL(r_size), INTENT(INOUT) :: fcer3d(nij1,nlev,nv3d)
  REAL(r_size), INTENT(INOUT) :: fcer2d(nij1,nv2d)
  REAL(r_size),PARAMETER :: tref = 280.0_r_size
  REAL(r_size),PARAMETER :: pref = 1000.0e+2_r_size
  REAL(r_size) :: ensmn3d(nij1,nlev,nv3d)
  REAL(r_size) :: ensmn2d(nij1,nv2d)
  REAL(r_size) :: tmptv(nij1,nlev)
  REAL(r_size) :: pdelta(nij1,nlev)
  REAL(r_size) :: sigweight(nij1,nlev)
  REAL(r_size) :: rinbv, cptr, qweight, rdtrpr
  INTEGER :: iret
  INTEGER :: i,j,k

  ! Calculate ensemble mean of forecast
  ensmn3d(:,:,:) = 0.0_r_size
  ensmn2d(:,:) = 0.0_r_size
  DO i=1,nbv
    ensmn3d(:,:,:) = ensmn3d(:,:,:) + fcst3d(:,:,i,:)
    ensmn2d(:,:) = ensmn2d(:,:) + fcst2d(:,i,:)
  END DO
  rinbv = 1.0_r_size / real(nbv,r_size)
  ensmn3d(:,:,:) = ensmn3d(:,:,:) * rinbv
  ensmn2d(:,:) = ensmn2d(:,:) * rinbv
  ! Calculate ensemble forecast perturbations
  DO i=1,nbv
    fcst3d(:,:,i,:) = fcst3d(:,:,i,:) - ensmn3d(:,:,:)
    fcst2d(:,i,:) = fcst2d(:,i,:) - ensmn2d(:,:)
  END DO
  ! Calculate delta-level pressure
  tmptv = ensmn3d(:,:,iv3d_t) * (1.0d0 + fvirt * ensmn3d(:,:,iv3d_q))
  call sigio_modprd(nij1,nij1,nlev,gfs_nvcoord,gfs_idvc,gfs_idsl,&
                    gfs_vcoord,iret,ensmn2d(:,iv2d_ps),tmptv,pd=pdelta)
  ! Compute combined weights
  DO k=1,nlev
    sigweight(:,k) = sqrt(pdelta(:,k) / ensmn2d(:,iv2d_ps)) * wg1(:)
  END DO
  ! Constants
  cptr = sqrt(cp/tref)
  qweight = sqrt(wmoist/(cp*tref))*hvap
  rdtrpr = sqrt(rd*tref)/pref
  ! For surface variables
  IF(tar_minlev <= 1) THEN
    DO i=1,nv2d
      IF(i == iv2d_ps) THEN
        !!! [(Rd*Tr)(dS/4pi)]^(1/2) * (ps'/Pr)
        fcer2d(:,i) = rdtrpr * wg1(:) * fcer2d(:,i)
        DO j=1,nbv
          fcst2d(:,j,i) = rdtrpr * wg1(:) * fcst2d(:,j,i)
        END DO
      ELSE
        fcer2d(:,i) = 0.0_r_size
        fcst2d(:,:,i) = 0.0_r_size
      END IF
    END DO
  ELSE
    fcer2d(:,:) = 0.0_r_size
    fcst2d(:,:,:) = 0.0_r_size
  END IF
  DO k=1,nlev
    IF(k > tar_maxlev .or. k < tar_minlev) THEN
      fcst3d(:,k,:,:) = 0.0_r_size
      fcer3d(:,k,:) = 0.0_r_size
      CYCLE
    END IF
    DO i=1,nv3d
      IF(i == iv3d_u .or. i == iv3d_v) THEN
        !!! [(dsigma)(dS/4pi)]^(1/2) * u'
        !!! [(dsigma)(dS/4pi)]^(1/2) * v'
        fcer3d(:,k,i) = sigweight(:,k) * fcer3d(:,k,i)
        DO j=1,nbv
          fcst3d(:,k,j,i) = sigweight(:,k) * fcst3d(:,k,j,i)
        END DO
      ELSE IF(i == iv3d_t) THEN
        !!! [(Cp/Tr)(dsigma)(dS/4pi)]^(1/2) * t'
        fcer3d(:,k,i) = cptr * sigweight(:,k) * fcer3d(:,k,i)
        DO j=1,nbv
          fcst3d(:,k,j,i) = cptr * sigweight(:,k) * fcst3d(:,k,j,i)
        END DO
      ELSE IF(i == iv3d_q) THEN
        !!! [(wg*L^2/Cp/Tr)(dsigma)(dS/4pi)]^(1/2) * q'
        fcer3d(:,k,i) = qweight * sigweight(:,k) * fcer3d(:,k,i)
        DO j=1,nbv
          fcst3d(:,k,j,i) = qweight * sigweight(:,k) * fcst3d(:,k,j,i)
        END DO
      ELSE
        fcer3d(:,k,i) = 0.0_r_size
        fcst3d(:,k,:,i) = 0.0_r_size
      END IF
    END DO
  END DO
  DO i=1,nij1
    IF(lon1(i) < tar_minlon .or. lon1(i) > tar_maxlon .or. &
         & lat1(i) < tar_minlat .or. lat1(i) > tar_maxlat) THEN
      fcer2d(i,:) = 0.0_r_size
      fcst2d(i,:,:) = 0.0_r_size
      fcer3d(i,:,:) = 0.0_r_size
      fcst3d(i,:,:,:) = 0.0_r_size
    END IF
  END DO
  RETURN
END SUBROUTINE lnorm

SUBROUTINE loc_advection(ua,va,uf,vf)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: ua(nij1,nlev)
  REAL(r_size),INTENT(IN) :: va(nij1,nlev)
  REAL(r_size),INTENT(IN) :: uf(nij1,nlev)
  REAL(r_size),INTENT(IN) :: vf(nij1,nlev)
  REAL(r_size) :: rad2deg, deg2rad
  REAL(r_size) :: coslat(nij1)
  INTEGER :: i,k
  ALLOCATE(lon2(nij1,nlev))
  ALLOCATE(lat2(nij1,nlev))
  deg2rad = pi/180.0_r_size
  rad2deg = locadv_rate*eft*3600.0_r_size*180.0_r_size/(pi*re)
  DO i=1,nij1
    coslat(i) = 1.0_r_size/cos(lat1(i)*deg2rad)
  END DO
  DO k=1,nlev
    DO i=1,nij1
      lon2(i,k) = lon1(i) - 0.5_r_size * (ua(i,k) + uf(i,k)) &
           & * coslat(i) * rad2deg
      lat2(i,k) = lat1(i) - 0.5_r_size * (va(i,k) + vf(i,k)) &
           & * rad2deg
      IF(lat2(i,k) > 90.0_r_size) THEN
        lat2(i,k) = 180.0_r_size - lat2(i,k)
        lon2(i,k) = lon2(i,k) + 180.0_r_size
      ELSE IF(lat2(i,k) < -90.0_r_size) THEN
        lat2(i,k) = -180.0_r_size - lat2(i,k)
        lon2(i,k) = lon2(i,k) + 180.0_r_size
      END IF
      IF(lon2(i,k) > 360.0_r_size) THEN
        lon2(i,k) = MOD(lon2(i,k),360.0_r_size)
      ELSE IF(lon2(i,k) < 0.0_r_size) THEN
        lon2(i,k) = MOD(lon2(i,k),360.0_r_size) + 360.0_r_size
      END IF
    END DO
  END DO
  RETURN
END SUBROUTINE loc_advection

SUBROUTINE print_obsense
  IMPLICIT NONE
  INTEGER,PARAMETER :: nreg = 3
  INTEGER :: regnh=1, regtr=2, regsh=3          ! Indices for the regions
  REAL(r_size),PARAMETER :: latbound=20._r_size ! Boundary latitude of TR
  INTEGER :: nobs_sense(nid_obs,nobtype+1,nreg)
  REAL(r_size) :: sumsense(nid_obs,nobtype+1,nreg)
  REAL(r_size) :: rate(nid_obs,nobtype+1,nreg)
  INTEGER :: nobs_t
  REAL(r_size) :: sumsense_t,rate_t
  INTEGER :: nob,oid,otype,ireg,iterm
  CHARACTER(len=2) :: charreg(nreg)
  CHARACTER(len=6) :: charotype
  CHARACTER(len=12) :: ofile(nterm)

  IF(nobs == 0) RETURN
  nobs_sense = 0
  sumsense = 0._r_size
  rate = 0._r_size
  charreg(regnh) ='NH'
  charreg(regtr) ='TR'
  charreg(regsh) ='SH'
  ofile(1)='osenseKE.dat'
  ofile(2)='osensePE.dat'
  ofile(3)='osenseME.dat'

  ! Binary output (in obs2 format)
  DO iterm = 1, 3
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',ofile(iterm)
    CALL write_obs2(ofile(iterm),nobs,obselm,obslon,obslat,obslev, &
                    obsdat,obserr,obstyp,obsdif,obsense(iterm,:),obsqc,0)
  END DO

  ! Loop over each observations
  iterm = 1
  DO nob=1,nobs
    ! Select observation elements
    oid = uid_obs(NINT(obselm(nob)))
    IF(oid <= 0 .OR. oid > nid_obs) CYCLE
    otype = NINT(obstyp(nob))
    ! Select observation types
    IF(otype <= 0 .OR. otype > nobtype+1) CYCLE
    ! Select observation regions
    IF(obslat(nob) > latbound) THEN
      ireg=regnh
    ELSE IF(obslat(nob) < -latbound) THEN
      ireg=regsh
    ELSE
      ireg=regtr
    END IF
    ! Sum up
    nobs_sense(oid,otype,ireg) = nobs_sense(oid,otype,ireg) + 1
    sumsense(oid,otype,ireg) = sumsense(oid,otype,ireg) + obsense(iterm,nob)
    IF(obsense(iterm,nob) < 0._r_size) THEN
      rate(oid,otype,ireg) = rate(oid,otype,ireg) + 1._r_size
    END IF
  END DO

  WRITE (6, '(A)') '============================================'
  WRITE (6, '(A,I10)') ' TOTAL NUMBER OF OBSERVATIONS:', nobs
  WRITE (6, '(A)') '============================================'
  WRITE (6, '(A)') '              nobs     dJ(KE)       +rate[%]'
  DO otype = 1,nobtype+1
    IF(otype <= nobtype) THEN
      charotype = obtypelist(otype)
    ELSE
      charotype = 'OTHERS'
    END IF
    nobs_t = SUM(nobs_sense(:,otype,:))
    IF(nobs_t > 0) THEN
      sumsense_t = SUM(sumsense(:,otype,:))
      rate_t = SUM(rate(:,otype,:)) / REAL(nobs_t,r_size) * 100._r_size
      WRITE (6, '(A)') '--------------------------------------------'
      WRITE (6,'(A,1x,A,1x,I8,1x,E12.5,1x,F8.2)') &
           & charotype,' TOTAL',nobs_t,sumsense_t,rate_t
    END IF
    DO ireg = 1,nreg
      DO oid = 1,nid_obs
        IF(nobs_sense(oid,otype,ireg) > 0) THEN
          rate_t = rate(oid,otype,ireg) &
               & / REAL(nobs_sense(oid,otype,ireg),r_size) * 100._r_size
          WRITE (6,'(A,1x,A,1x,A,1x,I8,1x,E12.5,1x,F8.2)') &
               & charotype,charreg(ireg),obelmlist(oid), &
               & nobs_sense(oid,otype,ireg), &
               & sumsense(oid,otype,ireg),   &
               & rate_t
        END IF
      END DO
    END DO
  END DO
  WRITE (6, '(A)') '============================================'

  RETURN
END SUBROUTINE print_obsense

SUBROUTINE destroy_obsense
  IMPLICIT NONE
  IF(allocated(obsense)) DEALLOCATE(obsense)
  IF(allocated(lon2)) DEALLOCATE(lon2)
  IF(allocated(lat2)) DEALLOCATE(lat2)
  RETURN
END SUBROUTINE destroy_obsense

END MODULE efso_tools
