MODULE common_obs_scale
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!   10/04/2012 Guo-Yuan Lien     modified for GFS model
!   07/25/2014 Guo-Yuan Lien     modified for SCALE model
!
!=======================================================================
!
! [LETKF observation format]
!   (In file, all stored in single precision float)
!
!  column  description
!     (1)  variable type (1 to nid_obs; see 'id_*_obs' parameters)
!     (2)  longitude (degree)
!     (3)  latitude (degree)
!     (4)  level/height
!            u,v,t,tv,q,rh: level (hPa)
!            ps: station elevation (m)
!     (5)  observation value
!            wind (m/s)
!            temperature (K)
!            specific humidity (kg/kg)
!            relative humidity (%)
!            surface pressure (hPa)
!     (6)  observation error
!            unit same as observation value
!     (7)  observation platform type (1 to nobtype+1; see 'obtypelist' array)
!
!  --- columns below only exist in obs2 (after the observation operator processing)
!     (8)  observation time relative to analysis time (hour)
!     (9)  h(x) observation in model background
!            unit same as observation value except
!            surface pressure (Pa)
!    (10)  quality control mark (1=pass; others=do not pass)
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_nml
  USE common_scale

  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: nid_obs=15 !H08
  INTEGER,PARAMETER :: nid_obs_varlocal=9 !H08
!
! conventional observations
!
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_tv_obs=3074
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331
!
! surface observations codes > 9999
!
  INTEGER,PARAMETER :: id_ps_obs=14593
  INTEGER,PARAMETER :: id_rain_obs=19999
  INTEGER,PARAMETER :: id_tclon_obs=99991  ! TC vital
  INTEGER,PARAMETER :: id_tclat_obs=99992  ! TC vital
  INTEGER,PARAMETER :: id_tcmip_obs=99993  ! TC vital
!
! radar observations
!
  INTEGER,PARAMETER :: id_radar_ref_obs=4001
  INTEGER,PARAMETER :: id_radar_vr_obs=4002
  INTEGER,PARAMETER :: id_radar_prh_obs=4003
!
! Himawari-8 (H08) observations
!
  INTEGER,PARAMETER :: id_H08IR_obs=8800

  INTEGER,PARAMETER :: elem_uid(nid_obs)= &
     (/id_u_obs, id_v_obs, id_t_obs, id_tv_obs, id_q_obs, id_rh_obs, &
       id_ps_obs, id_rain_obs, id_radar_ref_obs, id_radar_vr_obs, id_radar_prh_obs, &
       id_H08IR_obs, id_tclon_obs, id_tclat_obs, id_tcmip_obs/)

  CHARACTER(3),PARAMETER :: obelmlist(nid_obs)= &
     (/'  U', '  V', '  T', ' Tv', '  Q', ' RH', ' PS', 'PRC', 'REF', ' Vr', 'PRH',&
       'H08', 'TCX', 'TCY', 'TCP'/)

  CHARACTER(3),PARAMETER :: obelmlist_varlocal(nid_obs_varlocal)= &
     (/'WND', '  T', 'MOI', ' PS', 'PRC', 'TCV', 'REF', ' Vr', 'H08'/)

!  INTEGER,PARAMETER :: nobtype = 22
  INTEGER,PARAMETER :: nobtype = 24 ! H08
  CHARACTER(6),PARAMETER :: obtypelist(nobtype)= &
     (/'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR', &
       'VADWND', 'SATEMP', 'ADPSFC', 'SFCSHP', 'SFCBOG', &
       'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND', 'QKSWND', &
       'MSONET', 'GPSIPW', 'RASSDA', 'WDSATR', 'ASCATW', &
       'TMPAPR', 'PHARAD', 'H08IRB', 'TCVITL'/) ! H08

  INTEGER,PARAMETER :: max_obs_info_meta = 3 ! maximum array size for type(obs_info)%meta

  TYPE obs_info
    INTEGER :: nobs = 0
    INTEGER,ALLOCATABLE :: elm(:)
    REAL(r_size),ALLOCATABLE :: lon(:)
    REAL(r_size),ALLOCATABLE :: lat(:)
    REAL(r_size),ALLOCATABLE :: lev(:)
    REAL(r_size),ALLOCATABLE :: dat(:)
    REAL(r_size),ALLOCATABLE :: err(:)
    INTEGER,ALLOCATABLE :: typ(:)
    REAL(r_size),ALLOCATABLE :: dif(:)
    REAL(r_size) :: meta(max_obs_info_meta) = undef
  END TYPE obs_info

  !!!!!!
  TYPE obs_da_value
    INTEGER :: nobs = 0
    INTEGER,ALLOCATABLE :: set(:)
    INTEGER,ALLOCATABLE :: idx(:)
    REAL(r_size),ALLOCATABLE :: val(:)
    !
    ! obsda%lev array is used only for Himawari-8 assimilation.
    ! This array preserves the most sensitive height derived from transmittance outputs from RTTOV.
    ! For Himawari-8 assimilation, LETKF uses obsda%lev instead of obs%lev.
    ! 
#ifdef H08
    REAL(r_size),ALLOCATABLE :: lev(:) ! H08
    REAL(r_size),ALLOCATABLE :: val2(:) ! H08 ! clear sky BT (gues)
#endif
    REAL(r_size),ALLOCATABLE :: ensval(:,:)
    INTEGER,ALLOCATABLE :: qc(:)
    REAL(r_size),ALLOCATABLE :: ri(:)
    REAL(r_size),ALLOCATABLE :: rj(:)
  END TYPE obs_da_value
  !!!!!! need to add %err and %dat for obsda2 if they can be determined in letkf_obs.f90 !!!!!!

  INTEGER,PARAMETER :: nobsformats=3 ! H08
  CHARACTER(30) :: obsformat_name(nobsformats) = &
    (/'CONVENTIONAL ', 'RADAR        ', 'Himawari-8-IR'/)

  INTEGER,PARAMETER :: iqc_good=0
  INTEGER,PARAMETER :: iqc_gross_err=5
  INTEGER,PARAMETER :: iqc_ps_ter=10
  INTEGER,PARAMETER :: iqc_ref_low=11
  INTEGER,PARAMETER :: iqc_ref_mem=12
  INTEGER,PARAMETER :: iqc_radar_vhi=19
  INTEGER,PARAMETER :: iqc_out_vhi=20
  INTEGER,PARAMETER :: iqc_out_vlo=21
  INTEGER,PARAMETER :: iqc_out_h=22
  INTEGER,PARAMETER :: iqc_obs_bad=50
  INTEGER,PARAMETER :: iqc_otype=90
  INTEGER,PARAMETER :: iqc_time=91

  REAL(r_size),SAVE :: MIN_RADAR_REF
  REAL(r_size),SAVE :: RADAR_REF_THRES

!  REAL(r_size),PARAMETER :: UNDEF_OBS = 9.99d9           !Code that will be assigned to obs outside the domain.(so we don't need qc0 array)

!  TYPE obs_csort
!    INTEGER :: nobs = 0
!    LOGICAL :: sorted = .false.
!    INTEGER,ALLOCATABLE :: idx(:)
!    INTEGER,ALLOCATABLE :: nobsgrd(:,:)
!  END TYPE obs_csort

CONTAINS

!-----------------------------------------------------------------------
! Convert a raw obsID to a sequential obsID (1 - nid_obs)
!-----------------------------------------------------------------------
function uid_obs(id_obs)
  implicit none
  integer :: id_obs
  integer :: uid_obs

  select case(id_obs)
  case(id_u_obs)
    uid_obs = 1
  case(id_v_obs)
    uid_obs = 2
  case(id_t_obs)
    uid_obs = 3
  case(id_tv_obs)
    uid_obs = 4
  case(id_q_obs)
    uid_obs = 5
  case(id_rh_obs)
    uid_obs = 6
  case(id_ps_obs)
    uid_obs = 7
  case(id_rain_obs)
    uid_obs = 8
  case(id_radar_ref_obs)
    uid_obs = 9
  case(id_radar_vr_obs)
    uid_obs = 10
  case(id_radar_prh_obs)
    uid_obs = 11
  case(id_h08ir_obs) ! H08
    uid_obs = 12     ! H08
  case(id_tclon_obs)
    uid_obs = 13
  case(id_tclat_obs)
    uid_obs = 14
  case(id_tcmip_obs)
    uid_obs = 15
  case default
    uid_obs = -1     ! error
  end select
end function uid_obs
!-----------------------------------------------------------------------
! Convert a raw obsID to a sequential obsID for variable localization (1 - nid_obs_verlocal)
!-----------------------------------------------------------------------
function uid_obs_varlocal(id_obs)
  implicit none
  integer :: id_obs
  integer :: uid_obs_varlocal

  select case(id_obs)
  case(id_u_obs, id_v_obs)
    uid_obs_varlocal = 1
  case(id_t_obs, id_tv_obs)
    uid_obs_varlocal = 2
  case(id_q_obs, id_rh_obs)
    uid_obs_varlocal = 3
  case(id_ps_obs)
    uid_obs_varlocal = 4
  case(id_rain_obs)
    uid_obs_varlocal = 5
  case(id_tclon_obs, id_tclat_obs, id_tcmip_obs)
    uid_obs_varlocal = 6
  case(id_radar_ref_obs, id_radar_prh_obs)
    uid_obs_varlocal = 7
  case(id_radar_vr_obs)
    uid_obs_varlocal = 8
  case(id_h08ir_obs)      ! H08
    uid_obs_varlocal = 9  ! H08
  case default
    uid_obs_varlocal = -1 ! error
  end select
end function uid_obs_varlocal

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine set_common_obs_scale
  implicit none

  MIN_RADAR_REF = 10.0d0 ** (MIN_RADAR_REF_DBZ/10.0d0)
  RADAR_REF_THRES = 10.0d0 ** (RADAR_REF_THRES_DBZ/10.0d0)

  return
end subroutine set_common_obs_scale

!-----------------------------------------------------------------------
! Transformation from model variables to an observation
!
! stggrd: grid type of u and v
!  0: non-staggered grid
!  1: staggered grid
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY(elm,ri,rj,rk,lon,lat,v3d,v2d,yobs,qc,stggrd)
  use scale_mapproj, only: &
      MPRJ_rotcoef
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rk
  REAL(r_size),INTENT(IN) :: lon,lat
  REAL(r_size),INTENT(IN) :: v3d(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size),INTENT(IN) :: v2d(nlonh,nlath,nv2dd)
  REAL(r_size),INTENT(OUT) :: yobs
  INTEGER,INTENT(OUT) :: qc
  INTEGER,INTENT(IN),OPTIONAL :: stggrd
  REAL(r_size) :: u,v,t,q,topo
  REAL(RP) :: rotc(2)

  INTEGER :: stggrd_ = 0
  if (present(stggrd)) stggrd_ = stggrd

  yobs = undef
  qc = iqc_good

  SELECT CASE (elm)
  CASE(id_u_obs,id_v_obs)  ! U,V
    if (stggrd_ == 1) then
      CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri-0.5,rj,u)  !###### should modity itpl_3d to prevent '1.0' problem....??
      CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj-0.5,v)  !######
    else
      CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri,rj,u)
      CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj,v)
    end if
    call MPRJ_rotcoef(rotc,lon*deg2rad,lat*deg2rad)
    if (elm == id_u_obs) then
      yobs = u * rotc(1) - v * rotc(2)
    else
      yobs = u * rotc(2) + v * rotc(1)
    end if
  CASE(id_t_obs)  ! T
    CALL itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,yobs)
  CASE(id_tv_obs)  ! Tv
    CALL itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,yobs)
    CALL itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,q)
    yobs = yobs * (1.0d0 + fvirt * q)
  CASE(id_q_obs)  ! Q
    CALL itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,yobs)
  CASE(id_ps_obs) ! PS
    CALL itpl_2d(v2d(:,:,iv2dd_t2m),ri,rj,t)
    CALL itpl_2d(v2d(:,:,iv2dd_q2m),ri,rj,q)
    CALL itpl_2d(v2d(:,:,iv2dd_topo),ri,rj,topo)
    CALL itpl_2d(v2d(:,:,iv2dd_ps),ri,rj,yobs)
    call prsadj(yobs,rk-topo,t,q)
    if (abs(rk-topo) > PS_ADJUST_THRES) then
      write (6,'(A,F6.1)') 'warning: PS observation height adjustment exceeds the threshold. dz=', abs(rk-topo)
      qc = iqc_ps_ter
    end if
!  CASE(id_rain_obs) ! RAIN                        ############# (not finished)
!    CALL itpl_2d(v2d(:,:,iv2dd_rain),ri,rj,yobs) !#############
  CASE(id_rh_obs) ! RH
    CALL itpl_3d(v3d(:,:,:,iv3dd_rh),rk,ri,rj,yobs)
!  CASE(id_tclon_obs)
!    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
!    yobs = dummy(1)
!  CASE(id_tclat_obs)
!    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
!    yobs = dummy(2)
!  CASE(id_tcmip_obs)
!    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
!    yobs = dummy(3)
  CASE DEFAULT
    qc = iqc_otype
  END SELECT

  RETURN
END SUBROUTINE Trans_XtoY
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY_radar(elm,radar_lon,radar_lat,radar_z,ri,rj,rk,lon,lat,lev,v3d,v2d,yobs,qc,stggrd)
!  USE common_mpi
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rk,radar_lon,radar_lat,radar_z !!!!! Use only, ri, rj, rk eventually... (radar_lon,lat,z in ri,rj,rk)
  REAL(r_size),INTENT(IN) :: lon,lat,lev
  REAL(r_size),INTENT(IN) :: v3d(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size),INTENT(IN) :: v2d(nlonh,nlath,nv2dd)
  REAL(r_size),INTENT(OUT) :: yobs
  INTEGER,INTENT(OUT) :: qc
  INTEGER,INTENT(IN),OPTIONAL :: stggrd
  INTEGER :: stggrd_ = 0

  REAL(r_size) :: qvr,qcr,qrr,qir,qsr,qgr,ur,vr,wr,tr,pr,rhr
  REAL(r_size) :: dist , dlon , dlat , az , elev , radar_ref,radar_rv

!  integer :: ierr
!  REAL(r_dble) :: rrtimer00,rrtimer
!  rrtimer00 = MPI_WTIME()


  if (present(stggrd)) stggrd_ = stggrd


  yobs = undef
  qc = iqc_good

  if (stggrd_ == 1) then
    CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri-0.5,rj,ur)  !###### should modity itpl_3d to prevent '1.0' problem....??
    CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj-0.5,vr)  !######
  else
    CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri,rj,ur)
    CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj,vr)
  end if
  CALL itpl_3d(v3d(:,:,:,iv3dd_w),rk,ri,rj,wr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,tr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_p),rk,ri,rj,pr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,qvr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qc),rk,ri,rj,qcr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qr),rk,ri,rj,qrr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qi),rk,ri,rj,qir)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qs),rk,ri,rj,qsr)
  CALL itpl_3d(v3d(:,:,:,iv3dd_qg),rk,ri,rj,qgr)


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### Trans_XtoY_radar:itpl_3d:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !Compute az and elevation for the current observation.
  !Simple approach (TODO: implement a more robust computation)

  !Azimuth
  dlon=lon-radar_lon
  dlat=lat-radar_lat
  IF ( dlon == 0.0d0 .and. dlat == 0.0d0  )THEN
!      WRITE(6,*)'OOPS',dlon,dlat,lon,lat,radar_lon,radar_lat
    qc = iqc_out_h
    RETURN
  ELSE
    az = rad2deg*atan2(dlon*cos(radar_lat*deg2rad),dlat)
  ENDIF
  !WRITE(6,*)dlon,dlat,lon,lat,radar_lon(ityp),radar_lat(ityp),dlon*cos(radar_lat(ityp))
  !IF( abs(dlon) > maxdlon )maxdlon=abs(dlon)
  !IF( abs(dlat) > maxdlat )maxdlat=abs(dlat)
  !WRITE(6,*)maxdlon,maxdlat
  IF( az < 0) az = 360.0d0 + az
  !elevation
  CALL com_distll_1(lon,lat,radar_lon,radar_lat,dist)
  elev=rad2deg*atan2(lev-radar_z,dist)


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### Trans_XtoY_radar:radar_coord:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !Check that the azimuth and elevation angles are within the expected range.
  !Some grid points may be at the radar location.
  !IF( .NOT. ( az .GT. 0.0d0 .AND. az .LT. 360.0d0 ) )RETURN
  !IF( .NOT. ( elev .GT. 0.0d0 .AND. elev .LT. 90.0d0 ))RETURN



  !DEBUG---------------------------------------------------------------
  !WRITE(6,*)'RADAR lat ',radar_lat(ityp),' RADAR lon ',radar_lon(ityp)
  !WRITE(6,*)'OBS   lat ',dlat      ,' OBS   lon ',dlon
  !WRITE(6,*)'AZIMUTH   ',az
  !WRITE(6,*)'RADAR z   ',radar_z ,'  OBS   z   ',lev
  !WRITE(6,*)'ELEVATION ',elev
  !DEGUB---------------------------------------------------------------

  !WRITE(6,*)'BCRV',dlon,dlat,az,elev

  CALL calc_ref_vr(qvr,qcr,qrr,qir,qsr,qgr,ur,vr,wr,tr,pr,az,elev,radar_ref,radar_rv)


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### Trans_XtoY_radar:calc_ref_vr:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  SELECT CASE (elm)
  CASE(id_radar_ref_obs)
!!!!    if (radar_ref < MIN_RADAR_REF) then
!!!!      !In this case we will replace the observation by -RH
!!!!      !This allows us to use pseudo rh observations in some cases.
!!!!      !Later we will take the decision on what to do with these cases...

!!!!      CALL itpl_3d(v3d(:,:,:,iv3dd_rh),rk,ri,rj,rhr)
!!!!!      qc = 
!!!!      yobs = -rhr

!!!!    else                      !!!!!! --------- Pesudo RH: TO BE DONE...
    if (radar_ref < MIN_RADAR_REF) then
      qc = iqc_ref_low
      yobs = MIN_RADAR_REF_DBZ + LOW_REF_SHIFT  !!! even if the above qc is bad, still return the value
    else
      yobs = 10.0d0 * log10(radar_ref)
    end if
!!!!    end if
  CASE(id_radar_vr_obs)
    if (radar_ref < MIN_RADAR_REF) then
      qc = iqc_ref_low
    end if
    yobs = radar_rv  !!! even if the above qc is bad, still return the value
  CASE DEFAULT
    qc = iqc_otype
  END SELECT


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### Trans_XtoY_radar:conversion:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  RETURN
END SUBROUTINE Trans_XtoY_radar
!-----------------------------------------------------------------------
! TC center search
!  [AUTHORS:] T. Miyoshi and M. Kunii
!-----------------------------------------------------------------------
!SUBROUTINE tctrk(ps,t2,ri,rj,trk)
!  IMPLICIT NONE
!  INTEGER,PARAMETER :: isearch = 10 !search radius [grid points]
!  REAL(r_size),INTENT(IN) :: ps(nlon,nlat)
!  REAL(r_size),INTENT(IN) :: t2(nlon,nlat)
!  REAL(r_size),INTENT(IN) :: ri,rj
!  REAL(r_size),INTENT(OUT) :: trk(3) !1:lon, 2:lat, 3:minp
!  REAL(r_size) :: wk(100,3)
!  REAL(r_size) :: slp(nlon,nlat)
!  REAL(r_size) :: p1,p2,p3,p4,p5,a,c,d,xx,yy
!  INTEGER :: i,j,i0,i1,j0,j1,n

!  i0 = MAX(1,FLOOR(ri)-isearch)
!  i1 = MIN(nlon,CEILING(ri)+isearch)
!  j0 = MAX(1,FLOOR(rj)-isearch)
!  j1 = MIN(nlat,CEILING(rj)+isearch)
!  trk = undef

!  DO j=j0,j1
!    DO i=i0,i1
!      slp(i,j) = ps(i,j) * (1.0d0 - 0.0065d0 * phi0(i,j) / &
!        & (t2(i,j) + 0.0065d0 * phi0(i,j))) ** -5.257d0
!    END DO
!  END DO

!  n=0
!  DO j=j0+1,j1-1
!    DO i=i0+1,i1-1
!      IF(slp(i,j) > slp(i  ,j-1)) CYCLE
!      IF(slp(i,j) > slp(i  ,j+1)) CYCLE
!      IF(slp(i,j) > slp(i-1,j  )) CYCLE
!      IF(slp(i,j) > slp(i+1,j  )) CYCLE
!      IF(slp(i,j) > slp(i-1,j-1)) CYCLE
!      IF(slp(i,j) > slp(i+1,j-1)) CYCLE
!      IF(slp(i,j) > slp(i-1,j+1)) CYCLE
!      IF(slp(i,j) > slp(i+1,j+1)) CYCLE
!      p1 = slp(i,j)
!      p2 = slp(i-1,j)
!      p3 = slp(i+1,j)
!      p4 = slp(i,j-1)
!      p5 = slp(i,j+1)
!      c = (p3-p2)*0.5d0
!      d = (p5-p4)*0.5d0
!      a = (p2+p3+p4+p5)*0.25d0 - p1
!      IF(a == 0.0d0) CYCLE
!      xx = -0.5d0 * c / a
!      yy = -0.5d0 * d / a
!      n = n+1
!      wk(n,3) = p1 - a*(xx*xx + yy*yy)
!      wk(n,2) = lat(i,j) * (1.0d0 - yy) + lat(i,j+1) * yy
!      wk(n,1) = lon(i,j) * (1.0d0 - xx) + lon(i+1,j) * xx
!    END DO
!  END DO

!  j=1
!  IF(n > 1) THEN
!    a = wk(1,3)
!    DO i=2,n
!      IF(wk(i,3) < a) THEN
!        a = wk(i,3)
!        j = i
!      END IF
!    END DO
!  END IF
!  trk = wk(j,:)

!END SUBROUTINE tctrk
!-----------------------------------------------------------------------
! Compute relative humidity (RH)
!-----------------------------------------------------------------------
SUBROUTINE calc_rh(t,q,p,rh)
  IMPLICIT NONE
  REAL(r_size),PARAMETER :: t0=273.15d0
  REAL(r_size),PARAMETER :: e0c=6.11d0
  REAL(r_size),PARAMETER :: al=17.3d0
  REAL(r_size),PARAMETER :: bl=237.3d0
  REAL(r_size),PARAMETER :: e0i=6.1121d0
  REAL(r_size),PARAMETER :: ai=22.587d0
  REAL(r_size),PARAMETER :: bi=273.86d0
  REAL(r_size),INTENT(IN) :: t,q,p
  REAL(r_size),INTENT(OUT) :: rh
  REAL(r_size) :: e,es,tc

  e = q * p * 0.01d0 / (0.378d0 * q + 0.622d0)

  tc = t-t0
  IF(tc >= 0.0d0) THEN
    es = e0c * exp(al*tc/(bl+tc))
  ELSE IF(tc <= -15.d0) THEN
    es = e0i * exp(ai*tc/(bi+tc))
  ELSE
    es = e0c * exp(al*tc/(bl+tc)) * (15.0d0+tc)/15.0d0 &
       + e0i * exp(ai*tc/(bi+tc)) * (-tc) / 15.0d0
  END IF

  rh = e/es

  RETURN
END SUBROUTINE calc_rh
!-----------------------------------------------------------------------
! Pressure adjustment for a different height level
!-----------------------------------------------------------------------
SUBROUTINE prsadj(p,dz,t,q)
  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: p
  REAL(r_size),INTENT(IN) :: dz ! height difference (target - original) [m]
  REAL(r_size),INTENT(IN) :: t  ! temperature [K] at original level
  REAL(r_size),INTENT(IN) :: q  ! humidity [kg/kg] at original level
  REAL(r_size),PARAMETER :: gamma=5.0d-3 ! lapse rate [K/m]
  REAL(r_size) :: tv

  IF(dz /= 0) THEN
    tv = t * (1.0d0 + 0.608d0 * q)
    p = p * ((-gamma*dz+tv)/tv)**(gg/(gamma*rd)) !tv is at original level
!    p = p * (tv/(tv+gamma*dz))**(gg/(gamma*rd)) !tv is at target level
  END IF

  RETURN
END SUBROUTINE prsadj



!-----------------------------------------------------------------------
! Compute radar reflectivity and radial wind.
! Radial wind computations for certain methods depend on model reflectivity
! so both functions has been merged into a single one.
! First reflectivity is computed, and the the radial velocity is computed.
!-----------------------------------------------------------------------
SUBROUTINE calc_ref_vr(qv,qc,qr,qci,qs,qg,u,v,w,t,p,az,elev,ref,vr)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: qv        !Water vapor
  REAL(r_size), INTENT(IN) :: qc,qr     !Cloud and rain water
  REAL(r_size), INTENT(IN) :: qci,qs,qg !Cloud ice, snow and graupel
  REAL(r_size), INTENT(IN) :: t,p       !Temperature and pressure.
  REAL(r_size), INTENT(IN) :: u,v,w     !velocities with respecto to earth.
  REAL(r_size), INTENT(INOUT) :: ref      !Reflectivity.
  REAL(r_size)              :: ro
  REAL(r_size), INTENT(IN) :: az     !Azimuth respect to the radar.
  REAL(r_size), INTENT(IN) :: elev   !Elevation angle respect to the surface.
  REAL(r_size), INTENT(INOUT) :: vr    !Radial velocity.
  REAL(r_size)  :: qms , qmg !Melting species concentration (METHOD_REF_CALC 3)
  REAL(r_size)  :: qt        !Total condensate mixing ratio (METHOD_REF_CALC 1)
  REAL(r_size)  :: zr , zs , zg !Rain, snow and graupel's reflectivities.
  REAL(r_size)  :: wr , ws , wg !Rain, snow and graupel's mean terminal velocities.
  REAL(r_size)  :: wt           !Total mean terminal velocity.
  REAL(r_size)  :: nor, nos, nog !Rain, snow and graupel's intercepting parameters.
  REAL(r_size)  :: ror, ros, rog , roi !Rain, snow and graupel, ice densities.
  REAL(r_size)  :: a,b,c,d,Cd    !Constant for fall speed computations.
  REAL(r_size)  :: cf, pip , roo
  REAL(r_size)  :: ki2 , kr2
  REAL(r_size)  :: lr , ls , lg
  REAL(r_size)  :: tmp_factor , rofactor
  REAL(r_size)  :: p0
  REAL(r_size)  :: Fs, Fg , zms , zmg , fws , fwg !METHOD_REF_CALC 3
  REAL(r_size)  :: qrp , qsp , qgp
  REAL(r_size)  :: maxf                     !Maximum mixture relative concentration. (METHOD_REF_CALC 3)

  !Note: While equivalent reflectivity is assumed to be independent of the radar, in 
  !practice short wavelengths as those associated with K band radars migh frequently
  !experience Mie scattering. In that case, the equivalent reflectivity is not longer
  !radar independent and an appropiate relationship between the forecasted concentrations
  !and the reflectivity should be used.


  !REAL(r_size)  :: trqr,trqs,trqg
  ![P] Pa
  ![T] K
  ![q...] dimensionless
  ![ref] mm^6/m^3
  ![refdb] refdb
  ![U,V,W] m/s
  ![az, elev] degree
  ![vr] m/s

  !This model reflectivity won't be lower than this value.

  !Initialize reflectivities
  zr=0.0d0
  zs=0.0d0
  zg=0.0d0
  zms=0.0d0
  zmg=0.0d0
  ref=0.0d0

  !Compute air density (all methods use this)

  ro =  p / (rd * t)

  !Begin computation of reflectivity and vr

  if (METHOD_REF_CALC == 1) then

    !WRF method: See for example Sugimoto et al. Evaluation of the Performance of Ra
    !dial Velocity Assimilation with WRF-3DVAR System and Simulated Multiple-Doppler
    !Radar Data
    !Only rain is used to estimate the terminal velocity of hidrometeors.
    !Only rain is used to compute equivalent reflectivity.
    !Marshall-Palmer distributions are assumed to find the relationship bestween
    !concentration and equivalent reflectivity.
    ! Sun and Crook 1997 , 1998.
    !Derived for C-band radars.
    !Attenuation is not computed.

    !Reflectivity
    nor=8.0d6      ![m^-4]
    ror=1000.0d0   ![Kg/m3]
    pip=pi ** 1.75 !factor
    cf =10.0d18 * 72 !factor
    p0=1.0d5            !Reference pressure.

    qt=qr + qs + qg  !Assume that the total condensate is rain water
                     !But ignore cloud ice and cloud water

    IF( qt .GT. 0.0d0 )THEN
    ref = cf * ( ( ro * qt )**1.75 )
    ref = ref / ( pip * ( nor ** 0.75 ) * ( ror ** 1.75 ) )
   !ref= 2.04d4 *( ( ro * qt * 1.0d3 ) ** 1.75 ) !Original Sun and Crook expresion.
    ELSE
    ref=0.0d0
    ENDIF

    !Radial wind

    IF ( qt .GT. 0.0d0 )THEN
    a=(p0/p)**0.4
    wt = 5.40d0 * a * ( qt ** 0.125 )
    ELSE
    wt=0d0
    ENDIF
   !WRITE(6,*)qr,qs,qg


  else if (METHOD_REF_CALC == 2) then
    !Observation operator from Tong and Xue 2006, 2008 a and b.
    !Based on Smith et al 1975.
    !It includes reflectivity contribution by all the microphisical species.
    !is assumes Marshall and Palmer distributions.
    !Based on C band radars.
    !Attenuation is not computed.
    nor=8.0d6      ![m^-4]
    nos=3.0d6      ![m^-4]
    nog=4.0d4      ![m^-4]
    ror=1000.0d0   ![Kg/m3]
    ros=100.0d0    ![Kg/m3]
    rog=913.0d0    ![Kg/m3] 
    roi=917.0d0    ![Kg/m3]
    roo=1.0d0      ![Kg/m3] Surface air density.
    ki2=0.176d0    !Dielectric factor for ice.
    kr2=0.930d0    !Dielectric factor for water.
    pip=pi ** 1.75 !factor
    cf =1.0d18 * 720 !factor 

    IF( qr .GT. 0.0d0 )THEN
    zr= cf * ( ( ro * qr )**1.75 )
    zr= zr / ( pip * ( nor ** 0.75 ) * ( ror ** 1.75 ) )
    ENDIF
    !The contribution of snow depends on temperature (bright band effect)
    IF( qs .GT. 0.0d0 )THEN
    IF ( t <= 273.16 )THEN
     zs = cf * ki2 * ( ros ** 0.25 ) * ( ( ro * qs ) ** 1.75 )
     zs = zs / ( pip * kr2 * ( nos ** 0.75  ) * ( roi ** 2 ) )
    ELSE
    !WARNING: This formulation has to be checked the paper says that 
     !ros instead of roi should be used in thes expresion, but that 
     !leads to unrealistic jumps in the model derived reflectivity.
     zs = cf * ( ( ro * qs ) ** 1.75 )
     zs = zs / ( pip * ( nos ** 0.75 ) * ( roi ** 1.75 ) )
    ENDIF
    ENDIF

    !Only dry graupel contribution is ussed.
    IF( qg .GT. 0.0d0 )THEN
    zg= ( cf / ( pip * ( nog ** 0.75) * ( rog ** 1.75 ) ) ) ** 0.95
    zg= zg * ( ( ro * qg ) ** 1.6625 )

    !zg=(ki2/kr2)*( cf / ( pip * ( nog ** 0.75) * ( rog ** 1.75 ) ) ) 
    !zg= zg * ( ( ro * qg ) ** 1.75 )

    !if( qg * ro * 1000 .GT. 1.0d0 )THEN
    !  WRITE(6,*)qg, zg, ro, rd , t , p, 10*log10(zg)
    !endif
    ENDIF

    ref = zr + zs + zg

    !Compute reflectivity weigthed terminal velocity.
    !Lin et al 1983.
    IF( ref > 0.0d0 )THEN
    !There are hidrometeors, compute their terminal velocity.
    !Change units to be consistent with Lin et al 1983 and
    !to obtain wt in m/s
    nor=nor*1e-3      ![cm^-4]
    nos=nos*1e-3      ![cm^-4]
    nog=nog*1e-3      ![cm^-4]
    ror=ror*1e-3        ![g/cm3]
    ros=ros*1e-3        ![g/cm3]
    rog=rog*1e-3      ![g/cm3] 
    roo=roo*1e-3      ![g/cm3] Surface air density.
    ro= ro*1e-3

    a=2115d0   ![cm**1-b / s]
    b=0.8d0
    c=152.93d0 ![cm**1-b / s]
    d=0.25d0
    Cd=0.6d0

    lr= ( pi * ror * nor / ( ro * qr ) ) ** 0.25
    ls= ( pi * ros * nos / ( ro * qs ) ) ** 0.25
    lg= ( pi * rog * nog / ( ro * qg ) ) ** 0.25

    rofactor= ( roo / ro  ) ** 0.25
    CALL com_gamma( 4.0_r_size + b , tmp_factor )
    wr= a * tmp_factor / ( 6.0d0 * ( lr ** b ) )
    wr= 1.0d-2*wr * rofactor
    CALL com_gamma( 4.0_r_size + d , tmp_factor )
    ws= c * tmp_factor / ( 6.0d0 * ( ls ** d ) )
    ws= 1.0d-2*ws * rofactor
    CALL com_gamma( 4.5_r_size , tmp_factor )
    wg= tmp_factor * ( ( ( 4.0d0 * gg * 100.0d0 * rog )/( 3.0d0 * Cd * ro ) ) ** 0.5 )
    wg= 1.0d-2*wg / ( 6.0d0 * ( lg ** 0.5 ) )

    !Reflectivity weighted terminal velocity. 
    wt = ( wr * zr + ws * zs + wg * zg )/ ( zr + zs + zg )

    ELSE

    wt=0.0d0

    ENDIF

  else if (METHOD_REF_CALC == 3) then
    !Observation operator from Xue et al 2007
    !Asumes power law between Z and mass concentration of different 
    !hydrometeor categories. 
    !Derived for X-Band radars tacking into account Mie scattering.
    !Includes a computation of the attenuation. (k)

    MAXF=0.5d0

    !First we need to compute the mixtures between rain, snow and hail
    !Following Jung et al 2007 eq (2) and (3)
    Fg=0.0d0
    Fs=0.0d0
    fwg=0.0d0
    fws=0.0d0
    IF( qr .GT. 0.0d0 .AND. qg .GT. 0.0d0)THEN
      Fg=MAXF * ( min( qr/qg , qg/qr ) )**(1.0d0/3.0d0)
      fwg= qr / ( qr + qg )
    ENDIF
    IF( qr .GT. 0.0d0 .AND. qs .GT. 0.0d0)THEN
      Fs=MAXF * ( min( qr/qs , qs/qr ) )**(1.0d0/3.0d0)
      fws= qr / ( qr + qs )
    ENDIF


    !Correct the rain, snow and hail mixing ratios assuming
    !that we have a mixture due to melting.

    qrp=(1.0d0-Fs-Fg)*qr

    qsp=(1.0d0-Fs)*qs

    qgp=(1.0d0-Fg)*qg

    !Compute the new species resulting from melting.

    qms=Fs * (qr + qs) !Melting snow concentration.

    qmg=Fg * (qr + qg) !Melting hail concentration.

    !Compute reflectivities for each species including the melting species.

    IF( qrp .GT. 0.0d0)THEN
    zr= 2.53d4 * ( ro * qrp * 1.0d3 )**1.84
    ENDIF
    IF( qsp .GT. 0.0d0)THEN
    zs= 3.48d3 * ( ro * qsp * 1.0d3 )**1.66
    ENDIF
    IF( qgp .GT. 0.0d0)THEN
    zg= 8.18d4 * ( ro * qgp * 1.0d3 )**1.50
    ENDIF
    IF( qms .GT. 0.0d0 )THEN
    zms=( 0.00491 + 5.75*fws - 5.588*(fws**2) )*1.0d5
    zms= zms * ( ro * qms * 1.0d3 )**( 1.67 - 0.202*fws + 0.398*(fws**2) )

    ENDIF
    IF( qmg .GT. 0.0d0 )THEN
    zmg=( 0.809 + 10.13*fwg -5.98*(fwg**2) )*1.0d5
    zmg= zmg * ( ro * qmg * 1.0d3 )**( 1.48 + 0.0448*fwg - 0.0313*(fwg**2) )
    ENDIF

    ref = zr +  zg  + zs + zms + zmg

    !Compute reflectivity weigthed terminal velocity.
    !Lin et al 1983. (The distribution parameters are 
    !consistent with the work of Jung et al 2007)

    IF( ref > 0.0d0 )THEN
      !There are hidrometeors, compute their terminal velocity.
      !Units according to Lin et al 1983.
      nor=8.0d-2      ![cm^-4]
      nos=3.0d-2      ![cm^-4]
      nog=4.0d-4      ![cm^-4]
      ror=1.0d0        ![g/cm3]
      ros=0.1d0        ![g/cm3]
      rog=0.917d0      ![g/cm3] 
      roo=0.001d0      ![g/cm3] Surface air density.
      ro=1.0d-3 * ro
      a=2115d0   ![cm**1-b / s]
      b=0.8d0
      c=152.93d0 ![cm**1-b / s]
      d=0.25d0
      Cd=0.6d0

      rofactor= ( roo / ro  ) ** 0.5

      IF ( qr .GT. 0.0d0 )THEN
      CALL com_gamma( 4.0_r_size + b , tmp_factor )
      lr= ( pi * ror * nor / ( ro * qr ) ) ** 0.25
      wr= a * tmp_factor / ( 6.0d0 * ( lr ** b ) )
      wr= 1.0d-2 * wr * rofactor
      ELSE
      wr=0.0d0
      ENDIF

      IF( qs .GT. 0.0d0 )THEN
      ls= ( pi * ros * nos / ( ro * qs ) ) ** 0.25
      CALL com_gamma( 4.0_r_size + d , tmp_factor )
      ws= c * tmp_factor / ( 6.0d0 * ( ls ** d ) )
      ws= 1.0d-2 * ws * rofactor
      ELSE
      ws=0.0d0
      ENDIF

      IF ( qg .GT. 0.0d0 )THEN
      lg= ( pi * rog * nog / ( ro * qg ) ) ** 0.25
      CALL com_gamma( 4.5_r_size , tmp_factor )
      wg= tmp_factor * ( ( ( 4.0d0 * gg * 100.0d0 * rog )/( 3.0d0 * Cd * ro ) ) ** 0.5 )
      wg= 1.0d-2 * wg / ( 6.0d0 * ( lg ** 0.5 ) )
      ELSE
      wg=0.0d0
      ENDIF

      !Reflectivity weighted terminal velocity. 
      !The melting species are assumed to fail as their non-melting counterparts.
      !however this might not be the case for melting snow.
      wt = ( wr * zr + ws *  zs +  ws * zms + wg *  zg + wg * zmg ) / ( zr + zs + zg + zms + zmg )

    ELSE

      wt=0.0d0

    ENDIF

  else  !IF OVER DIFFERENT OPTIONS

    WRITE(6,*)'ERROR: Not recognized method for radar reflectivity and wind computation'
    STOP

  end if ! [METHOD_REF_CALC == ?]


  !Compute radial velocity
  !WRITE(6,*)'ICRV',u,v,w,wt,az,elev
  vr = u * cos(elev*deg2rad) * sin(az*deg2rad)
  vr = vr + v * cos(elev*deg2rad) * cos(az*deg2rad)
  IF( USE_TERMINAL_VELOCITY )THEN
    vr = vr + (w - wt)*sin(elev*deg2rad)
  ELSE
    vr = vr + (w)*sin(elev*deg2rad)
  ENDIF

  !WRITE(6,*) u , v , w
  !WRITE(6,*) wt , vr
  !WRITE(6,*) elev , az , deg2rad
  !STOP


  RETURN
END SUBROUTINE calc_ref_vr



!-----------------------------------------------------------------------
! Coordinate conversion (find rk in pressure levels)
!
! rk = 0.0d0  : surface observation
!-----------------------------------------------------------------------
SUBROUTINE phys2ijk(p_full,elem,ri,rj,rlev,rk,qc)
  use scale_grid_index, only: &
      KHALO
  IMPLICIT NONE

  REAL(r_size),INTENT(IN) :: p_full(nlevh,nlonh,nlath)
  INTEGER,INTENT(IN) :: elem
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rlev ! pressure levels (for 3D variable only)
  REAL(r_size),INTENT(OUT) :: rk
  INTEGER,INTENT(OUT) :: qc
  REAL(r_size) :: ak
!  REAL(r_size) :: lnps(nlonh,nlath)
  REAL(r_size) :: lnps(nlevh,nlonh,nlath)
  REAL(r_size) :: plev(nlevh)
  REAL(r_size) :: ptmp
  INTEGER :: i,j,k, ii, jj, ks

  qc = iqc_good
!
! rlev -> rk
!
  if (ri < 1.0d0 .or. ri > nlonh .or. rj < 1.0d0 .or. rj > nlath) then
    write (6,'(A)') 'warning: observation is outside of the horizontal domain'
    rk = undef
    qc = iqc_out_h
    return
  end if
  !
  IF(elem > 9999) THEN ! surface observation
    rk = rlev
  ELSE
    !
    ! horizontal interpolation
    !
    i = CEILING(ri)
    j = CEILING(rj)
    !
    ! Find the lowest valid level
    !
    ks = 1+KHALO
    do jj = j-1, j
      do ii = i-1, i


!print *, p_full(:,ii,jj)

        DO k=1+KHALO,nlev+KHALO
          if (p_full(k,ii,jj) >= 0.0d0) exit
        END DO
        if (k > ks) ks = k
      end do
    end do
!    DO k=1+KHALO,nlev+KHALO
!!      IF(i <= nlon+IHALO) THEN
!        lnps(i-1:i,j-1:j) = LOG(p_full(k,i-1:i,j-1:j))
!!      ELSE
!!        lnps(i-1,j-1:j) = LOG(p_full(k,i-1,j-1:j))
!!        lnps(1,j-1:j) = LOG(p_full(k,1,j-1:j))
!!      END IF
!      CALL itpl_2d(lnps,ri,rj,plev(k))
!    END DO

    lnps(:,i-1:i,j-1:j) = LOG(p_full(:,i-1:i,j-1:j))
    call itpl_2d_column(lnps,ri,rj,plev)

    !
    ! Log pressure
    !
    rk = LOG(rlev)
    !
    ! determine if rk is within bound.
    !
    IF(rk < plev(nlev+KHALO)) THEN
      call itpl_2d(p_full(nlev+KHALO,:,:),ri,rj,ptmp)
      write(6,'(A,F8.1,A,F8.1,A,I5)') 'warning: observation is too high: ptop=', ptmp, ', lev=', rlev, ', elem=', elem
      rk = undef
      qc = iqc_out_vhi
      RETURN
    END IF
    IF(rk > plev(ks)) THEN
      call itpl_2d(p_full(ks,:,:),ri,rj,ptmp)
!print *, ks, rk, plev(ks)
      write(6,'(A,F8.1,A,F8.1,A,I5)') 'warning: observation is too low: pbottom=', ptmp, ', lev=', rlev, ', elem=', elem
      rk = undef
      qc = iqc_out_vlo

!print *, plev
!stop
!print *, elem, ri, rj, rlev, rk, qc


      RETURN
    END IF
    !
    ! find rk
    !
    DO k=ks+1,nlev+KHALO
      IF(plev(k) < rk) EXIT ! assuming descending order of plev
    END DO
    ak = (rk - plev(k-1)) / (plev(k) - plev(k-1))
    rk = REAL(k-1,r_size) + ak
  END IF

  RETURN
END SUBROUTINE phys2ijk
!-----------------------------------------------------------------------
! Coordinate conversion (find rk in height levels)
!
! rk = 0.0d0  : surface observation
!-----------------------------------------------------------------------
SUBROUTINE phys2ijkz(z_full,ri,rj,rlev,rk,qc)
  use scale_grid_index, only: &
      KHALO
!  use common_mpi
  IMPLICIT NONE

  REAL(r_size),INTENT(IN) :: z_full(nlevh,nlonh,nlath)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rlev ! height levels
  REAL(r_size),INTENT(OUT) :: rk
  INTEGER,INTENT(OUT) :: qc
  REAL(r_size) :: ak
  REAL(r_size) :: zlev(nlevh)
  REAL(r_size) :: ztmp
  INTEGER :: i,j,k, ii, jj, ks

!  integer :: ierr
!  REAL(r_dble) :: rrtimer00,rrtimer
!  rrtimer00 = MPI_WTIME()


  qc = iqc_good
!
! rlev -> rk
!
  if (ri < 1.0d0 .or. ri > nlonh .or. rj < 1.0d0 .or. rj > nlath) then
    write (6,'(A)') 'warning: observation is outside of the horizontal domain'
    rk = undef
    qc = iqc_out_h
    return
  end if


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:check_domain:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !
  ! horizontal interpolation
  !
  i = CEILING(ri)
  j = CEILING(rj)
  !
  ! Find the lowest valid level
  !
  ks = 1+KHALO
  do jj = j-1, j
    do ii = i-1, i
      DO k=1+KHALO,nlev+KHALO
        if (z_full(k,ii,jj) > -300.0d0 .and. z_full(k,ii,jj) < 10000.0d0) exit
      END DO
      if (k > ks) ks = k
    end do
  end do


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:find_lowesr_level:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


!  DO k=1+KHALO,nlev+KHALO
!    CALL itpl_2d(z_full(k,:,:),ri,rj,zlev(k))
!  END DO

  call itpl_2d_column(z_full,ri,rj,zlev)


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:itpl_2d:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !
  ! determine if rlev is within bound.
  !
  IF(rlev > zlev(nlev+KHALO)) THEN
    call itpl_2d(z_full(nlev+KHALO,:,:),ri,rj,ztmp)
    write(6,'(A,F8.1,A,F8.1)') 'warning: observation is too high: ztop=', ztmp, ', lev=', rlev
    rk = undef
    qc = iqc_out_vhi
    RETURN
  END IF
  IF(rlev < zlev(ks)) THEN
    call itpl_2d(z_full(ks,:,:),ri,rj,ztmp)
    write(6,'(A,F8.1,A,F8.1)') 'warning: observation is too low: zbottom=', ztmp, ', lev=', rlev
    rk = undef
    qc = iqc_out_vlo
    RETURN
  END IF


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:check_height:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  !
  ! find rk
  !
  DO k=ks+1,nlev+KHALO
    IF(zlev(k) > rlev) EXIT ! assuming ascending order of zlev
  END DO
  ak = (rlev - zlev(k-1)) / (zlev(k) - zlev(k-1))
  rk = REAL(k-1,r_size) + ak


!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,F18.10)') '###### phys2ijkz:calc_rk:',rrtimer-rrtimer00
!  rrtimer00=rrtimer


  RETURN
END SUBROUTINE phys2ijkz
!-----------------------------------------------------------------------
! Coordinate conversion
!-----------------------------------------------------------------------
SUBROUTINE phys2ij(rlon,rlat,rig,rjg)
  use scale_grid, only: &
      GRID_CXG, &
      GRID_CYG, &
      DX, &
      DY
  use scale_mapproj, only: &
      MPRJ_lonlat2xy
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: rlon
  REAL(r_size),INTENT(IN) :: rlat
  REAL(r_size),INTENT(OUT) :: rig
  REAL(r_size),INTENT(OUT) :: rjg
!
! rlon,rlat -> ri,rj
!
  call MPRJ_lonlat2xy(rlon*pi/180.0_r_size,rlat*pi/180.0_r_size,rig,rjg)
  rig = (rig - GRID_CXG(1)) / DX + 1.0d0
  rjg = (rjg - GRID_CYG(1)) / DY + 1.0d0

  RETURN
END SUBROUTINE phys2ij

SUBROUTINE ij2phys(rig,rjg,rlon,rlat)
  use scale_grid, only: &
      GRID_CXG, &
      GRID_CYG, &
      DX, &
      DY
  use scale_mapproj, only: &
      MPRJ_xy2lonlat
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: rig
  REAL(r_size),INTENT(IN) :: rjg
  REAL(r_size),INTENT(OUT) :: rlon ! (deg)
  REAL(r_size),INTENT(OUT) :: rlat ! (deg)
  REAL(r_size) :: x, y ! (m)
!
! ri,rj -> rlon,rlat
!
  x = (rig - 1.0d0) * DX + GRID_CXG(1) 
  y = (rjg - 1.0d0) * DY + GRID_CYG(1) 

  call MPRJ_xy2lonlat(x,y,rlon,rlat)

  rlon = rlon * rad2deg
  rlat = rlat * rad2deg

  RETURN
END SUBROUTINE ij2phys
!
!-----------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------
SUBROUTINE itpl_2d(var,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlonh,nlath)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
     & + var(i  ,j-1) *    ai  * (1-aj) &
     & + var(i-1,j  ) * (1-ai) *    aj  &
     & + var(i  ,j  ) *    ai  *    aj

  RETURN
END SUBROUTINE itpl_2d

SUBROUTINE itpl_2d_column(var,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlevh,nlonh,nlath)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5(nlevh)
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  var5(:) = var(:,i-1,j-1) * (1-ai) * (1-aj) &
        & + var(:,i  ,j-1) *    ai  * (1-aj) &
        & + var(:,i-1,j  ) * (1-ai) *    aj  &
        & + var(:,i  ,j  ) *    ai  *    aj

  RETURN
END SUBROUTINE itpl_2d_column

SUBROUTINE itpl_3d(var,ri,rj,rk,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlevh,nlonh,nlath)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rk
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj,ak
  INTEGER :: i,j,k

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)
  k = CEILING(rk)
  ak = rk - REAL(k-1,r_size)

  var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
     & + var(i  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
     & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
     & + var(i  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
     & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
     & + var(i  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
     & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
     & + var(i  ,j  ,k  ) *    ai  *    aj  *    ak

  RETURN
END SUBROUTINE itpl_3d
!-----------------------------------------------------------------------
! Monitor observation departure by giving the v3dg,v2dg data
!-----------------------------------------------------------------------
subroutine monit_obs(v3dg,v2dg,obs,obsda,topo,nobs,bias,rmse,monit_type)
  use scale_process, only: &
      PRC_myrank
  use scale_grid_index, only: &
      IHALO, JHALO, KHALO, &
      IS, IE, JS, JE, KS, KE, KA
  use scale_grid, only: &
      GRID_CZ, &
      GRID_FZ
  use scale_comm, only: &
      COMM_vars8, &
      COMM_wait



  implicit none

  REAL(RP),intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),intent(in) :: v2dg(nlon,nlat,nv2d)
  type(obs_info),intent(in) :: obs(OBS_IN_NUM)
  type(obs_da_value),intent(in) :: obsda
  real(r_size),intent(in) :: topo(nlon,nlat)
  INTEGER,INTENT(OUT) :: nobs(nid_obs)
  REAL(r_size),INTENT(OUT) :: bias(nid_obs)
  REAL(r_size),INTENT(OUT) :: rmse(nid_obs)
  LOGICAL,INTENT(OUT) :: monit_type(nid_obs)

  REAL(r_size) :: v3dgh(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size) :: v2dgh(nlonh,nlath,nv2dd)
  integer :: n,i,j,k,proc,iv3d,iv2d
  real(r_size) :: ri,rj,rk

  real(r_size) :: oelm(obsda%nobs)
  real(r_size) :: ohx(obsda%nobs)
  integer :: oqc(obsda%nobs)

!  real(r_size),allocatable :: oelm(:)
!  real(r_size),allocatable :: ohx(:)
!  integer,allocatable :: oqc(:)

!  REAL(r_size) :: timer
!  INTEGER :: ierr

  real(r_size) :: ztop

#ifdef H08
! -- for Himawari-8 obs --
  INTEGER :: nprof_H08 ! num of H08 obs
  REAL(r_size),ALLOCATABLE :: ri_H08(:),rj_H08(:)
  REAL(r_size),ALLOCATABLE :: lon_H08(:),lat_H08(:)
  REAL(r_size),ALLOCATABLE :: tmp_ri_H08(:),tmp_rj_H08(:)
  REAL(r_size),ALLOCATABLE :: tmp_lon_H08(:),tmp_lat_H08(:)
  INTEGER,ALLOCATABLE :: n2prof(:) ! obs num 2 prof num

  REAL(r_size),ALLOCATABLE :: yobs_H08(:),plev_obs_H08(:)
  REAL(r_size),ALLOCATABLE :: yobs_H08_clr(:)
  REAL(r_size),ALLOCATABLE :: CA(:) ! (Okamoto et al., 2014QJRMS)
  INTEGER :: ns
  INTEGER,ALLOCATABLE :: qc_H08(:)
#endif

! -- for TC vital assimilation --
!  INTEGER :: obs_idx_TCX, obs_idx_TCY, obs_idx_TCP ! obs index
!  INTEGER :: bTC_proc ! the process where the background TC is located.
! bTC: background TC in each subdomain
! bTC(1,:) : tcx (m), bTC(2,:): tcy (m), bTC(3,:): mslp (Pa)
!  REAL(r_size),ALLOCATABLE :: bTC(:,:)
!  REAL(r_size),ALLOCATABLE :: bufr(:,:)
!  REAL(r_size) :: bTC_mslp

!CALL MPI_BARRIER(MPI_COMM_a,ierr)
!CALL CPU_TIME(timer)
!if (myrank == 0) print *, '######', timer

  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_u) = v3dg(:,:,:,iv3d_u)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_v) = v3dg(:,:,:,iv3d_v)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_w) = v3dg(:,:,:,iv3d_w)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_t) = v3dg(:,:,:,iv3d_t)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_p) = v3dg(:,:,:,iv3d_p)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_q) = v3dg(:,:,:,iv3d_q)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qc) = v3dg(:,:,:,iv3d_qc)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qr) = v3dg(:,:,:,iv3d_qr)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qi) = v3dg(:,:,:,iv3d_qi)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qs) = v3dg(:,:,:,iv3d_qs)
  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_qg) = v3dg(:,:,:,iv3d_qg)
!  v3dgh(1+KHALO:nlev+KHALO,1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv3dd_rh) =

  ztop = GRID_FZ(KE) - GRID_FZ(KS-1)
!$OMP PARALLEL DO PRIVATE(j,i,k)
  do j = 1, nlat
    do i = 1, nlon
      do k = 1, nlev
        v3dgh(k+KHALO,i+IHALO,j+JHALO,iv3dd_hgt) = (ztop - topo(i,j)) / ztop * GRID_CZ(k+KHALO) + topo(i,j)
      end do
    enddo
  enddo
!$OMP END PARALLEL DO

  !!! use the 1st level as the surface (although it is not)
  v2dgh(:,:,iv2dd_topo) = v3dgh(1+KHALO,:,:,iv3dd_hgt)

  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_ps) = v3dg(1,:,:,iv3d_p)
!  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_rain) =
  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_u10m) = v3dg(1,:,:,iv3d_u)
  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_v10m) = v3dg(1,:,:,iv3d_v)
  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_t2m) = v3dg(1,:,:,iv3d_t)
  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_q2m) = v3dg(1,:,:,iv3d_q)

#ifdef H08
  v2dgh(1+IHALO:nlon+IHALO,1+JHALO:nlat+JHALO,iv2dd_skint) = v3dg(1,:,:,iv3d_t)

  !!! assume the point where terrain height is less than 10 m is the ocean. T.Honda (02/09/2016)
!$OMP PARALLEL DO PRIVATE(j,i)
  do j = 1, nlat
    do i = 1, nlon
      v2dgh(i+IHALO,j+JHALO,iv2dd_lsmask) = min(max(topo(i,j) - 10.0d0, 0.0d0), 1.0d0)
    enddo
  enddo
!$OMP END PARALLEL DO
#endif

  do iv3d = 1, nv3dd
!!!!!$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j  = JS, JE
      do i  = IS, IE
        v3dgh(   1:KS-1,i,j,iv3d) = v3dgh(KS,i,j,iv3d)
        v3dgh(KE+1:KA,  i,j,iv3d) = v3dgh(KE,i,j,iv3d)
      end do
    end do
  end do

  do iv3d = 1, nv3dd
    call COMM_vars8( v3dgh(:,:,:,iv3d), iv3d )
  end do
  do iv3d = 1, nv3dd
    call COMM_wait ( v3dgh(:,:,:,iv3d), iv3d )
  end do

  do iv2d = 1, nv2dd
    call COMM_vars8( v2dgh(:,:,iv2d), iv2d )
  end do
  do iv2d = 1, nv2dd
    call COMM_wait ( v2dgh(:,:,iv2d), iv2d )
  end do



!    allocate (oelm(obsda%nobs))
!    allocate (ohx(obsda%nobs))
!    allocate (oqc(obsda%nobs))

  oqc = -1

!  obs_idx_TCX = -1
!  obs_idx_TCY = -1
!  obs_idx_TCP = -1

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,ri,rj,rk)
  do n = 1, obsda%nobs

    if (obsda%qc(n) /= iqc_good) write(6, *) '############', obsda%qc(n)

    oelm(n) = obs(obsda%set(n))%elm(obsda%idx(n))

!    select case (int(oelm(n)))
!    case (id_tclon_obs)
!      obs_idx_TCX = n
!      cycle
!    case (id_tclat_obs)
!      obs_idx_TCY = n
!      cycle
!    case (id_tcmip_obs)
!      obs_idx_TCP = n
!      cycle
!    end select

#ifdef H08
    if(int(oelm(n)) == id_H08IR_obs)cycle
#endif

    call rij_g2l_auto(proc,obsda%ri(n),obsda%rj(n),ri,rj)
    if (PRC_myrank /= proc) then
      write(6, *) '############ Error!', PRC_myrank,proc,obsda%ri(n),obsda%rj(n),ri,rj
      cycle
!      stop
    end if

    if (DEPARTURE_STAT_T_RANGE <= 0.0d0 .or. &
        abs(obs(obsda%set(n))%dif(obsda%idx(n))) <= DEPARTURE_STAT_T_RANGE) then

      oqc(n) = iqc_otype

      select case (obs(obsda%set(n))%elm(obsda%idx(n)))
      case(id_u_obs,id_v_obs,id_t_obs,id_tv_obs,id_q_obs,id_ps_obs) !,id_rh_obs)
        call phys2ijk(v3dgh(:,:,:,iv3dd_p),obs(obsda%set(n))%elm(obsda%idx(n)), &
                      ri,rj,obs(obsda%set(n))%lev(obsda%idx(n)),rk,oqc(n))
        if (oqc(n) == iqc_good) then
          call Trans_XtoY(obs(obsda%set(n))%elm(obsda%idx(n)),ri,rj,rk, &
                          obs(obsda%set(n))%lon(obsda%idx(n)),obs(obsda%set(n))%lat(obsda%idx(n)), &
                          v3dgh,v2dgh,ohx(n),oqc(n),stggrd=1)
        end if

      case(id_radar_ref_obs,id_radar_vr_obs,id_radar_prh_obs)
        if (DEPARTURE_STAT_RADAR) then
          call phys2ijkz(v3dgh(:,:,:,iv3dd_hgt),ri,rj,obs(obsda%set(n))%lev(obsda%idx(n)),rk,oqc(n))
          if (oqc(n) == iqc_good) then
            call Trans_XtoY_radar(obs(obsda%set(n))%elm(obsda%idx(n)),obs(obsda%set(n))%meta(1), &
                                  obs(obsda%set(n))%meta(2),obs(obsda%set(n))%meta(3),ri,rj,rk, &
                                  obs(obsda%set(n))%lon(obsda%idx(n)),obs(obsda%set(n))%lat(obsda%idx(n)), &
                                  obs(obsda%set(n))%lev(obsda%idx(n)),v3dgh,v2dgh,ohx(n),oqc(n),stggrd=1)
            if (oqc(n) == iqc_ref_low) oqc(n) = iqc_good ! when process the observation operator, we don't care if reflectivity is too small
          end if
        end if
      end select

      if (oqc(n) == iqc_good) then
        ohx(n) = obs(obsda%set(n))%dat(obsda%idx(n)) - ohx(n)
      end if
!write (6, '(2I6,2F8.2,4F12.4,I3)') obs(obsda%set(n))%elm(obsda%idx(n)), obs(obsda%set(n))%typ(obsda%idx(n)), obs(obsda%set(n))%lon(obsda%idx(n)), obs(obsda%set(n))%lat(obsda%idx(n)), obs(obsda%set(n))%lev(obsda%idx(n)), obs(obsda%set(n))%dat(obsda%idx(n)), obs(obsda%set(n))%err(obsda%idx(n)), ohx(n), oqc(n)

    end if

  end do ! [ n = 1, obsda%nobs ]
!$OMP END PARALLEL DO


#ifdef H08
!
! -- Count the number of the Himawari-8 obs "location" (=num of prof).
!    Then, Trans_XtoY_H08 will be called without openMP.
!

  if (DEPARTURE_STAT_H08) then !-- [DEPARTURE_STAT_H08]

    ALLOCATE(tmp_ri_H08(obsda%nobs))
    ALLOCATE(tmp_rj_H08(obsda%nobs))
    ALLOCATE(tmp_lon_H08(obsda%nobs))
    ALLOCATE(tmp_lat_H08(obsda%nobs))
    ALLOCATE(n2prof(obsda%nobs))

    n2prof = 0
    nprof_H08 = 0
    do n = 1, obsda%nobs
      oelm(n) = obs(obsda%set(n))%elm(obsda%idx(n))
      if(oelm(n) /= id_H08IR_obs)cycle

      call rij_g2l_auto(proc,obsda%ri(n),obsda%rj(n),ri,rj)
      if (PRC_myrank /= proc) then
        write(6, *) '############ Error from H08 monitor!', PRC_myrank,proc
        cycle
      end if

      if(nprof_H08 > 1)then
        if((tmp_ri_H08(nprof_H08)==ri) .and. (tmp_ri_H08(nprof_H08)==ri))then
          n2prof(n) = nprof_H08
          cycle
        else
          nprof_H08 = nprof_H08 + 1
          tmp_ri_H08(nprof_H08) = ri
          tmp_rj_H08(nprof_H08) = rj
          tmp_lon_H08(nprof_H08) = obs(obsda%set(n))%lon(obsda%idx(n))
          tmp_lat_H08(nprof_H08) = obs(obsda%set(n))%lat(obsda%idx(n))
          n2prof(n) = nprof_H08
        endif
      else ! nprof_H08 <= 1
        nprof_H08 = nprof_H08 + 1
        tmp_ri_H08(nprof_H08) = ri
        tmp_rj_H08(nprof_H08) = rj
        tmp_lon_H08(nprof_H08) = obs(obsda%set(n))%lon(obsda%idx(n))
        tmp_lat_H08(nprof_H08) = obs(obsda%set(n))%lat(obsda%idx(n))
        n2prof(n) = nprof_H08
      endif
    end do ! [ n = 1, obsda%nobs ]

    IF(nprof_H08 >=1)THEN ! [nprof_H08 >=1]
      ALLOCATE(ri_H08(nprof_H08))
      ALLOCATE(rj_H08(nprof_H08))
      ALLOCATE(lon_H08(nprof_H08))
      ALLOCATE(lat_H08(nprof_H08))

      ri_H08 = tmp_ri_H08(1:nprof_H08)
      rj_H08 = tmp_rj_H08(1:nprof_H08)
      lon_H08 = tmp_lon_H08(1:nprof_H08)
      lat_H08 = tmp_lat_H08(1:nprof_H08)

      DEALLOCATE(tmp_ri_H08, tmp_rj_H08)
      DEALLOCATE(tmp_lon_H08, tmp_lat_H08)

      ALLOCATE(yobs_H08(nprof_H08*nch))
      ALLOCATE(yobs_H08_clr(nprof_H08*nch))
      ALLOCATE(CA(nprof_H08*nch))
      ALLOCATE(plev_obs_H08(nprof_H08*nch))
      ALLOCATE(qc_H08(nprof_H08*nch))


      CALL Trans_XtoY_H08(nprof_H08,ri_H08,rj_H08,&
                          lon_H08,lat_H08,v3dgh,v2dgh,&
                          yobs_H08,plev_obs_H08,&
                          qc_H08,stggrd=1,yobs_H08_clr=yobs_H08_clr)

      ! yobs here should be positive!!
      yobs_H08 = abs(yobs_H08)

      write (6, '(A)')"MEAN-HIMAWARI-8-STATISTICS"

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,ns)
      do n = 1, obsda%nobs
        oelm(n) = obs(obsda%set(n))%elm(obsda%idx(n))
        if(oelm(n) /= id_H08IR_obs)cycle
  
        ns = (n2prof(n) - 1) * nch + nint(obsda%lev(n) - 6.0) 

        if (DEPARTURE_STAT_T_RANGE <= 0.0d0 .or. & 
          abs(obs(obsda%set(n))%dif(obsda%idx(n))) <= DEPARTURE_STAT_T_RANGE) then
!          oqc(n) = iqc_otype

          oqc(n) = qc_H08(ns)
          ohx(n) = yobs_H08(ns)
!          if(plev_obs_H08(ns) < H08_LIMIT_LEV) oqc(n) = iqc_obs_bad

          CA(n) =  (abs(yobs_H08(ns) - yobs_H08_clr(ns)) & ! CM
                   +  abs(obs(obsda%set(n))%dat(obsda%idx(n)) - yobs_H08_clr(ns)) & ! CO
                   &) * 0.5d0 
                   

          if(oqc(n) == iqc_good) then
            ohx(n) = obs(obsda%set(n))%dat(obsda%idx(n)) - ohx(n) 
            write (6, '(A,2I6,2F8.2,5F11.4,I6,F10.4)')"H08-O-A-B",&
                  obs(obsda%set(n))%elm(obsda%idx(n)), &
                  nint(obsda%lev(n)), & ! obsda%lev includes the band num.
                  obs(obsda%set(n))%lon(obsda%idx(n)), &
                  obs(obsda%set(n))%lat(obsda%idx(n)), &
                  ohx(n), &! O-A
                  obsda%val(n), &! O-B
                  plev_obs_H08(ns), &
                  obs(obsda%set(n))%dat(obsda%idx(n)), &
                  obs(obsda%set(n))%err(obsda%idx(n)), &
                  oqc(n),&
                  CA(n) 
          endif

        endif ! [DEPARTURE_STAT_T_RANGE]
      end do ! [ n = 1, obsda%nobs ]
!$OMP END PARALLEL DO

      DEALLOCATE(yobs_H08, yobs_H08_clr, plev_obs_H08, qc_H08, CA)
    ENDIF ! [nprof_H08 >=1]

    IF(ALLOCATED(n2prof)) DEALLOCATE(n2prof)
    IF(ALLOCATED(tmp_ri_H08)) DEALLOCATE(tmp_ri_H08)
    IF(ALLOCATED(tmp_rj_H08)) DEALLOCATE(tmp_rj_H08)
    IF(ALLOCATED(tmp_lon_H08)) DEALLOCATE(tmp_lon_H08)
    IF(ALLOCATED(tmp_lat_H08)) DEALLOCATE(tmp_lat_H08)
  endif !-- [DEPARTURE_STAT_H08]

#endif

! ###  -- TC vital assimilation -- ###
!  if (obs_idx_TCX > 0 .and. obs_idx_TCY > 0 .and. obs_idx_TCP > 0 .and.&
!    obs(obsda%set(obs_idx_TCX))%dif(obsda%idx(obs_idx_TCX)) == &
!    obs(obsda%set(obs_idx_TCY))%dif(obsda%idx(obs_idx_TCY)) .and. &
!    obs(obsda%set(obs_idx_TCY))%dif(obsda%idx(obs_idx_TCY)) == &
!    obs(obsda%set(obs_idx_TCP))%dif(obsda%idx(obs_idx_TCP)) .and. & 
!    (DEPARTURE_STAT_T_RANGE <= 0.0d0 .or. &
!    abs(obs(obsda%set(obs_idx_TCX))%dif(obsda%idx(obs_idx_TCX))) <= DEPARTURE_STAT_T_RANGE))then
!
!    allocate(bTC(3,0:MEM_NP-1))
!    allocate(bufr(3,0:MEM_NP-1))
!
!    bTC = 9.99d33
!    bufr = 9.99d33
!
!    call search_tc_subdom(obsda%ri(obs_idx_TCX),obsda%rj(obs_idx_TCX),v2dg,bTC(1,PRC_myrank),bTC(2,PRC_myrank),bTC(3,PRC_myrank))
!
!    CALL MPI_BARRIER(MPI_COMM_d,ierr)
!    CALL MPI_ALLREDUCE(bTC,bufr,3*MEM_NP,MPI_r_size,MPI_MIN,MPI_COMM_d,ierr)
!    bTC = bufr
!
!    deallocate(bufr)
!
!
!    ! Assume MSLP of background TC is lower than 1100 (hPa). 
!    bTC_mslp = 1100.0d2
!    do n = 0, MEM_NP - 1
!      write(6,'(3e20.5)')bTC(1,n),bTC(2,n),bTC(3,n) ! debug
!      if (bTC(3,n) < bTC_mslp ) then
!        bTC_mslp = bTC(3,n)
!        bTC_proc = n
!      endif
!    enddo ! [ n = 0, MEM_NP - 1]
!
!    do n = 1, 3
!      if(n==1) i = obs_idx_TCX
!      if(n==2) i = obs_idx_TCY
!      if(n==3) i = obs_idx_TCP
!
!      ohx(i) = bTC(n,bTC_proc)
!      oqc(i) = iqc_otype
!      if(bTC_MSLP < 1100.0d2) oqc(i) = iqc_good
!
!      if (oqc(i) == iqc_good) then
!        ohx(i) = obs(obsda%set(i))%dat(obsda%idx(i)) - ohx(i)
!      end if
!    enddo
!
!  endif ! [DEPARTURE_STAT_T_RANGE]


  call monit_dep(obsda%nobs,oelm,ohx,oqc,nobs,bias,rmse)

  monit_type = .false.
  monit_type(uid_obs(id_u_obs)) = .true.
  monit_type(uid_obs(id_v_obs)) = .true.
  monit_type(uid_obs(id_t_obs)) = .true.
  monit_type(uid_obs(id_tv_obs)) = .true.
  monit_type(uid_obs(id_q_obs)) = .true.
!  monit_type(uid_obs(id_rh_obs)) = .true.
  monit_type(uid_obs(id_ps_obs)) = .true.
  if (DEPARTURE_STAT_RADAR) then
    monit_type(uid_obs(id_radar_ref_obs)) = .true.
    monit_type(uid_obs(id_radar_vr_obs)) = .true.
!    monit_type(uid_obs(id_radar_prh_obs)) = .true.
  end if
  if (DEPARTURE_STAT_H08) then
    monit_type(uid_obs(id_H08IR_obs)) = .true.
  end if


!    deallocate (oelm)
!    deallocate (ohx)
!    deallocate (oqc)


!CALL MPI_BARRIER(MPI_COMM_a,ierr)
!CALL CPU_TIME(timer)
!if (myrank == 0) print *, '######', timer



end subroutine monit_obs
!-----------------------------------------------------------------------
! Monitor departure
!  ofmt: output format
!    0: U,V,T(Tv),Q,RH,PS (default)
!    1: U,V,T(Tv),Q,RH,PS,RAIN
!-----------------------------------------------------------------------
SUBROUTINE monit_dep(nn,elm,dep,qc,nobs,bias,rmse)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elm(nn)
  REAL(r_size),INTENT(IN) :: dep(nn)
  INTEGER,INTENT(IN) :: qc(nn)
  INTEGER,INTENT(OUT) :: nobs(nid_obs)
  REAL(r_size),INTENT(OUT) :: bias(nid_obs)
  REAL(r_size),INTENT(OUT) :: rmse(nid_obs)
  INTEGER :: n,i,ielm

  nobs = 0
  rmse = 0.0d0
  bias = 0.0d0

  DO n=1,nn
    IF(qc(n) /= iqc_good) CYCLE

    ielm = NINT(elm(n))
    if (ielm == id_tv_obs) then ! compute Tv as T
      ielm = id_t_obs
    end if
    i = uid_obs(ielm)

    nobs(i) = nobs(i) + 1
    bias(i) = bias(i) + dep(n)
    rmse(i) = rmse(i) + dep(n)**2
  END DO

  DO i = 1, nid_obs
    IF(nobs(i) == 0) THEN
      bias(i) = undef
      rmse(i) = undef
    ELSE
      bias(i) = bias(i) / REAL(nobs(i),r_size)
      rmse(i) = SQRT(rmse(i) / REAL(nobs(i),r_size))
    END IF
  END DO

  RETURN
END SUBROUTINE monit_dep
!-----------------------------------------------------------------------
! Monitor departure
!-----------------------------------------------------------------------
SUBROUTINE monit_print(nobs,bias,rmse,monit_type)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nobs(nid_obs)
  REAL(r_size),INTENT(IN) :: bias(nid_obs)
  REAL(r_size),INTENT(IN) :: rmse(nid_obs)
  LOGICAL,INTENT(IN),OPTIONAL :: monit_type(nid_obs)

  character(12) :: var_show(nid_obs)
  character(12) :: nobs_show(nid_obs)
  character(12) :: bias_show(nid_obs)
  character(12) :: rmse_show(nid_obs)

  integer :: i, itv, n
  character(4) :: nstr
  character(12) :: tmpstr(nid_obs)
  character(12) :: tmpstr2(nid_obs)

  logical :: monit_type_(nid_obs)

  monit_type_ = .true.
  if (present(monit_type)) monit_type_ = monit_type

  n = 0
  itv = uid_obs(id_tv_obs)
  do i = 1, nid_obs
    if (monit_type_(i) .and. i /= itv) then
      n = n + 1
      write(var_show(n),'(A12)') obelmlist(i)
      write(nobs_show(n),'(I12)') nobs(i)
      if (nobs(i) > 0) then
        write(bias_show(n),'(ES12.3)') bias(i)
        write(rmse_show(n),'(ES12.3)') rmse(i)
      else
        write(bias_show(n),'(A12)') 'N/A'
        write(rmse_show(n),'(A12)') 'N/A'
      end if
    end if
  end do
  write(nstr, '(I4)') n
  tmpstr(1:n) = '============'
  tmpstr2(1:n) = '------------'

  WRITE(6,'(A,' // trim(nstr) // 'A)') '======', tmpstr(1:n)
  WRITE(6,'(6x,' // trim(nstr) // 'A)')          var_show(1:n)
  WRITE(6,'(A,' // trim(nstr) // 'A)') '------', tmpstr2(1:n)
  WRITE(6,'(A,' // trim(nstr) // 'A)') 'BIAS  ', bias_show(1:n)
  WRITE(6,'(A,' // trim(nstr) // 'A)') 'RMSE  ', rmse_show(1:n)
  WRITE(6,'(A,' // trim(nstr) // 'A)') 'NUMBER', nobs_show(1:n)
  WRITE(6,'(A,' // trim(nstr) // 'A)') '======', tmpstr(1:n)

  RETURN
END SUBROUTINE monit_print
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_info_allocate(obs)
  IMPLICIT NONE
  TYPE(obs_info),INTENT(INOUT) :: obs

  call obs_info_deallocate(obs)

  ALLOCATE( obs%elm (obs%nobs) )
  ALLOCATE( obs%lon (obs%nobs) )
  ALLOCATE( obs%lat (obs%nobs) )
  ALLOCATE( obs%lev (obs%nobs) )
  ALLOCATE( obs%dat (obs%nobs) )
  ALLOCATE( obs%err (obs%nobs) )
  ALLOCATE( obs%typ (obs%nobs) )
  ALLOCATE( obs%dif (obs%nobs) )

  obs%elm = 0
  obs%lon = 0.0d0
  obs%lat = 0.0d0
  obs%lev = 0.0d0
  obs%dat = 0.0d0
  obs%err = 0.0d0
  obs%typ = 0
  obs%dif = 0.0d0

  RETURN
END SUBROUTINE obs_info_allocate
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_info_deallocate(obs)
  IMPLICIT NONE
  TYPE(obs_info),INTENT(INOUT) :: obs

  IF(ALLOCATED(obs%elm)) DEALLOCATE(obs%elm)
  IF(ALLOCATED(obs%lon)) DEALLOCATE(obs%lon)
  IF(ALLOCATED(obs%lat)) DEALLOCATE(obs%lat)
  IF(ALLOCATED(obs%lev)) DEALLOCATE(obs%lev)
  IF(ALLOCATED(obs%dat)) DEALLOCATE(obs%dat)
  IF(ALLOCATED(obs%err)) DEALLOCATE(obs%err)
  IF(ALLOCATED(obs%typ)) DEALLOCATE(obs%typ)
  IF(ALLOCATED(obs%dif)) DEALLOCATE(obs%dif)

  RETURN
END SUBROUTINE obs_info_deallocate
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_da_value_allocate(obs,member)
  IMPLICIT NONE
  TYPE(obs_da_value),INTENT(INOUT) :: obs
  INTEGER,INTENT(IN) :: member

  call obs_da_value_deallocate(obs)

  ALLOCATE( obs%set    (obs%nobs) )
  ALLOCATE( obs%idx    (obs%nobs) )
  ALLOCATE( obs%val    (obs%nobs) )
#ifdef H08
  ALLOCATE( obs%lev    (obs%nobs) ) ! H08
  ALLOCATE( obs%val2    (obs%nobs) ) ! H08
#endif
  ALLOCATE( obs%qc     (obs%nobs) )
  ALLOCATE( obs%ri     (obs%nobs) )
  ALLOCATE( obs%rj     (obs%nobs) )

  obs%idx = 0
  obs%val = 0.0d0
#ifdef H08
  obs%lev = 0.0d0 ! H08
  obs%val2 = 0.0d0 ! H08
#endif
  obs%qc = 0
  obs%ri = 0.0d0
  obs%rj = 0.0d0

  if (member > 0) then
    ALLOCATE( obs%ensval (member,obs%nobs) )
    obs%ensval = 0.0d0
  end if

  RETURN
END SUBROUTINE obs_da_value_allocate
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE obs_da_value_deallocate(obs)
  IMPLICIT NONE
  TYPE(obs_da_value),INTENT(INOUT) :: obs

  IF(ALLOCATED(obs%set    )) DEALLOCATE(obs%set    )
  IF(ALLOCATED(obs%idx    )) DEALLOCATE(obs%idx    )
  IF(ALLOCATED(obs%val    )) DEALLOCATE(obs%val    )
#ifdef H08
  IF(ALLOCATED(obs%lev    )) DEALLOCATE(obs%lev    ) ! H08
  IF(ALLOCATED(obs%val2   )) DEALLOCATE(obs%val2   ) ! H08
#endif
  IF(ALLOCATED(obs%ensval )) DEALLOCATE(obs%ensval )
  IF(ALLOCATED(obs%qc     )) DEALLOCATE(obs%qc     )
  IF(ALLOCATED(obs%ri     )) DEALLOCATE(obs%ri     )
  IF(ALLOCATED(obs%rj     )) DEALLOCATE(obs%rj     )

  RETURN
END SUBROUTINE obs_da_value_deallocate
!-----------------------------------------------------------------------
! Basic modules for observation input
!-----------------------------------------------------------------------
SUBROUTINE get_nobs(cfile,nrec,nn)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nrec
  INTEGER,INTENT(OUT) :: nn
  REAL(r_sngl),ALLOCATABLE :: wk(:)
  INTEGER :: ios
  INTEGER :: iu,iv,it,iq,irh,ips,itc
  INTEGER :: iunit
  LOGICAL :: ex

  ALLOCATE(wk(nrec))
  nn = 0
  iu = 0
  iv = 0
  it = 0
  iq = 0
  irh = 0
  ips = 0
  itc = 0
  iunit=91
  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
    DO
      READ(iunit,IOSTAT=ios) wk
      IF(ios /= 0) EXIT
!      SELECT CASE(NINT(wk(1)))
!      CASE(id_u_obs)
!        iu = iu + 1
!      CASE(id_v_obs)
!        iv = iv + 1
!      CASE(id_t_obs,id_tv_obs)
!        it = it + 1
!      CASE(id_q_obs)
!        iq = iq + 1
!      CASE(id_rh_obs)
!        irh = irh + 1
!      CASE(id_ps_obs)
!        ips = ips + 1
!      CASE(id_tclon_obs)
!        itc = itc + 1
!      END SELECT
      nn = nn + 1
    END DO
!    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
!    WRITE(6,'(A12,I10)') '          U:',iu
!    WRITE(6,'(A12,I10)') '          V:',iv
!    WRITE(6,'(A12,I10)') '      T(Tv):',it
!    WRITE(6,'(A12,I10)') '          Q:',iq
!    WRITE(6,'(A12,I10)') '         RH:',irh
!    WRITE(6,'(A12,I10)') '         Ps:',ips
!    WRITE(6,'(A12,I10)') '   TC TRACK:',itc
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF
  DEALLOCATE(wk)

  RETURN
END SUBROUTINE get_nobs

SUBROUTINE read_obs(cfile,obs)
  use scale_mapproj, only: &
      MPRJ_lonlat2xy
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(INOUT) :: obs
  REAL(r_sngl) :: wk(8)
  REAL(r_size) :: x, y
  INTEGER :: n,iunit

!  call obs_info_allocate(obs)

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,obs%nobs
    READ(iunit) wk
    SELECT CASE(NINT(wk(1)))
    CASE(id_u_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_v_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_t_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_tv_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_q_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_ps_obs)
      wk(5) = wk(5) * 100.0 ! hPa -> Pa
      wk(6) = wk(6) * 100.0 ! hPa -> Pa
    CASE(id_rh_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = wk(5) * 0.01 ! percent input
      wk(6) = wk(6) * 0.01 ! percent input
    CASE(id_tcmip_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = wk(5) * 100.0 ! hPa -> Pa
      wk(6) = real(OBSERR_TCP,kind=r_sngl)
    CASE(id_tclon_obs)
      call MPRJ_lonlat2xy(REAL(wk(2),kind=r_size)*pi/180.0_r_size,&
                          REAL(wk(3),kind=r_size)*pi/180.0_r_size,&
                          x,y)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = real(x,kind=r_sngl)
      wk(6) = real(OBSERR_TCX,kind=r_sngl)
    CASE(id_tclat_obs)
      call MPRJ_lonlat2xy(REAL(wk(2),kind=r_size)*pi/180.0_r_size,&
                          REAL(wk(3),kind=r_size)*pi/180.0_r_size,&
                          x,y)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = real(y,kind=r_sngl)
      wk(6) = real(OBSERR_TCY,kind=r_sngl)
    END SELECT
    obs%elm(n) = NINT(wk(1))
    obs%lon(n) = REAL(wk(2),r_size)
    obs%lat(n) = REAL(wk(3),r_size)
    obs%lev(n) = REAL(wk(4),r_size)
    obs%dat(n) = REAL(wk(5),r_size)
    obs%err(n) = REAL(wk(6),r_size)
    obs%typ(n) = NINT(wk(7))
    obs%dif(n) = REAL(wk(8),r_size)
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs

SUBROUTINE write_obs(cfile,obs,append,missing)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(IN) :: obs
  LOGICAL,INTENT(IN),OPTIONAL :: append
  LOGICAL,INTENT(IN),OPTIONAL :: missing
  LOGICAL :: append_
  LOGICAL :: missing_
  REAL(r_sngl) :: wk(8)
  INTEGER :: n,iunit

  iunit=92
  append_ = .false.
  IF(present(append)) append_ = append
  missing_ = .true.
  IF(present(missing)) missing_ = missing

  IF(append_) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='append')
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  END IF
  DO n=1,obs%nobs
    if (missing_ .or. abs(obs%dat(n) - undef) > tiny(wk(5))) then
      wk(1) = REAL(obs%elm(n),r_sngl)
      wk(2) = REAL(obs%lon(n),r_sngl)
      wk(3) = REAL(obs%lat(n),r_sngl)
      wk(4) = REAL(obs%lev(n),r_sngl)
      wk(5) = REAL(obs%dat(n),r_sngl)
      wk(6) = REAL(obs%err(n),r_sngl)
      wk(7) = REAL(obs%typ(n),r_sngl)
      wk(8) = REAL(obs%dif(n),r_sngl)
      SELECT CASE(NINT(wk(1)))
      CASE(id_u_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_v_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_t_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_tv_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_q_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
      CASE(id_ps_obs)
        wk(5) = wk(5) * 0.01 ! Pa -> hPa
        wk(6) = wk(6) * 0.01 ! Pa -> hPa
      CASE(id_rh_obs)
        wk(4) = wk(4) * 0.01 ! Pa -> hPa
        wk(5) = wk(5) * 100.0 ! percent output
        wk(6) = wk(6) * 100.0 ! percent output
      CASE(id_tcmip_obs)
        wk(5) = wk(5) * 0.01 ! Pa -> hPa
        wk(6) = wk(6) * 0.01 ! Pa -> hPa
      END SELECT
      WRITE(iunit) wk
    end if
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs

SUBROUTINE read_obs_da(cfile,obs,im)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_da_value),INTENT(INOUT) :: obs
  INTEGER,INTENT(IN) :: im
#ifdef H08
!  REAL(r_sngl) :: wk(7) ! H08
  REAL(r_sngl) :: wk(8) ! H08
#else
  REAL(r_sngl) :: wk(6) ! H08
#endif
  INTEGER :: n,iunit

!  call obs_da_value_allocate(obs)

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,obs%nobs
    READ(iunit) wk
    obs%set(n) = NINT(wk(1))
    obs%idx(n) = NINT(wk(2))  !!!!!! will overflow......
    if (im == 0) then
      obs%val(n) = REAL(wk(3),r_size)
    else
      obs%ensval(im,n) = REAL(wk(3),r_size)
    end if
    obs%qc(n) = NINT(wk(4))
    obs%ri(n) = REAL(wk(5),r_size)
    obs%rj(n) = REAL(wk(6),r_size)
#ifdef H08
    obs%lev(n) = REAL(wk(7),r_size) ! H08
    obs%val2(n) = REAL(wk(8),r_size) ! H08
#endif
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs_da

SUBROUTINE write_obs_da(cfile,obs,im,append)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_da_value),INTENT(IN) :: obs
  INTEGER,INTENT(IN) :: im
  LOGICAL,INTENT(IN),OPTIONAL :: append
  LOGICAL :: append_
#ifdef H08
!  REAL(r_sngl) :: wk(7) ! H08
  REAL(r_sngl) :: wk(8) ! H08
#else
  REAL(r_sngl) :: wk(6) 
#endif
  INTEGER :: n,iunit

  iunit=92
  append_ = .false.
  IF(present(append)) append_ = append
  IF(append_) THEN
    IF(obs%nobs <= 0) RETURN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='append',STATUS='replace')
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential',STATUS='replace')
  END IF
  DO n=1,obs%nobs
    wk(1) = REAL(obs%set(n),r_sngl)
    wk(2) = REAL(obs%idx(n),r_sngl)  !!!!!! will overflow......
    if (im == 0) then
      wk(3) = REAL(obs%val(n),r_sngl)
    else
      wk(3) = REAL(obs%ensval(im,n),r_sngl)
    end if
    wk(4) = REAL(obs%qc(n),r_sngl)
    wk(5) = REAL(obs%ri(n),r_sngl)
    wk(6) = REAL(obs%rj(n),r_sngl)
#ifdef H08
    wk(7) = REAL(obs%lev(n),r_sngl) ! H08
    wk(8) = REAL(obs%val2(n),r_sngl) ! H08
#endif
    WRITE(iunit) wk
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs_da

SUBROUTINE get_nobs_radar(cfile,nn,radarlon,radarlat,radarz)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(OUT) :: nn
  REAL(r_sngl) :: wk(7),tmp
  INTEGER :: ios
  INTEGER :: ir,iv
  INTEGER :: iunit
  LOGICAL :: ex
  REAL(r_size),INTENT(OUT) :: radarlon,radarlat,radarz

  nn = 0
  iv = 0
  ir = 0
  iunit=91

  radarlon=0.0d0
  radarlat=0.0d0
  radarz  =0.0d0

  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
    READ(iunit,IOSTAT=ios)tmp
    IF(ios /= 0) THEN
      WRITE(6,'(2A)') cfile,': Reading error -- skipped'
      RETURN
    END IF
    radarlon=REAL(tmp,r_size)
    READ(iunit,IOSTAT=ios)tmp
    IF(ios /= 0) THEN 
      WRITE(6,'(2A)') cfile,': Reading error -- skipped'
      RETURN
    END IF
    radarlat=REAL(tmp,r_size)
    READ(iunit,IOSTAT=ios)tmp
    IF(ios /= 0) THEN 
      WRITE(6,'(2A)') cfile,': Reading error -- skipped'
      RETURN
    END IF
    radarz=REAL(tmp,r_size)
    DO
      READ(iunit,IOSTAT=ios) wk
      IF(ios /= 0) EXIT
      SELECT CASE(NINT(wk(1)))
      CASE(id_radar_ref_obs)
        ir = ir + 1
      CASE(id_radar_vr_obs)
        iv = iv + 1
      END SELECT
      nn = nn + 1
    END DO
    WRITE(6,*)' RADAR FILE ', cfile
    WRITE(6,*)' RADAR LON = ',radarlon
    WRITE(6,*)' RADAR LAT = ',radarlat
    WRITE(6,*)' RADAR Z   = ',radarz
    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
    WRITE(6,'(A12,I10)') '   REFLECTIVITY:',ir
    WRITE(6,'(A12,I10)') 'RADIAL VELOCITY:',iv
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF

  RETURN
END SUBROUTINE get_nobs_radar

SUBROUTINE read_obs_radar(cfile,obs)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(INOUT) :: obs
  REAL(r_sngl) :: wk(7)
  REAL(r_sngl) :: tmp
  INTEGER :: n,iunit,ios

!  call obs_info_allocate(obs)

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  READ(iunit, iostat=ios)tmp
  IF(ios /= 0) RETURN
  READ(iunit, iostat=ios)tmp
  IF(ios /= 0) RETURN
  READ(iunit, iostat=ios)tmp
  IF(ios /= 0) RETURN
  DO n=1,obs%nobs
    READ(iunit) wk
    obs%elm(n) = NINT(wk(1))
    obs%lon(n) = REAL(wk(2),r_size)
    obs%lat(n) = REAL(wk(3),r_size)
    obs%lev(n) = REAL(wk(4),r_size)
    obs%dat(n) = REAL(wk(5),r_size)
    obs%err(n) = REAL(wk(6),r_size)
!    obs%typ(n) = NINT(wk(7))
    obs%typ(n) = 22
    obs%dif(n) = 0.0d0
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs_radar

SUBROUTINE write_obs_radar(cfile,obs,append,missing)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(IN) :: obs
  LOGICAL,INTENT(IN),OPTIONAL :: append
  LOGICAL,INTENT(IN),OPTIONAL :: missing
  LOGICAL :: append_
  LOGICAL :: missing_
  REAL(r_sngl) :: wk(7)
  INTEGER :: n,iunit

  iunit=92
  append_ = .false.
  IF(present(append)) append_ = append
  missing_ = .true.
  IF(present(missing)) missing_ = missing

  IF(append_) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='append')
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  END IF
  WRITE(iunit) REAL(obs%meta(1),r_sngl)
  WRITE(iunit) REAL(obs%meta(2),r_sngl)
  WRITE(iunit) REAL(obs%meta(3),r_sngl)
  DO n=1,obs%nobs
    if (missing_ .or. abs(obs%dat(n) - undef) > tiny(wk(5))) then
      wk(1) = REAL(obs%elm(n),r_sngl)
      wk(2) = REAL(obs%lon(n),r_sngl)
      wk(3) = REAL(obs%lat(n),r_sngl)
      wk(4) = REAL(obs%lev(n),r_sngl)
      wk(5) = REAL(obs%dat(n),r_sngl)
      wk(6) = REAL(obs%err(n),r_sngl)
      wk(7) = REAL(obs%typ(n),r_sngl)
      WRITE(iunit) wk
    end if
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs_radar

subroutine read_obs_all(obs)
  implicit none

  type(obs_info), intent(out) :: obs(OBS_IN_NUM)
  integer :: iof
  logical :: ex

  do iof = 1, OBS_IN_NUM
    inquire (file=trim(OBS_IN_NAME(iof)), exist=ex)
    if (.not. ex) then
      write(6,*) 'WARNING: FILE ',trim(OBS_IN_NAME(iof)),' NOT FOUND'


      obs(iof)%nobs = 0
      call obs_info_allocate(obs(iof))


      cycle
    end if

    select case (OBS_IN_FORMAT(iof))
    case (1)
      call get_nobs(trim(OBS_IN_NAME(iof)),8,obs(iof)%nobs)
    case (2)
      call get_nobs_radar(trim(OBS_IN_NAME(iof)), obs(iof)%nobs, obs(iof)%meta(1), obs(iof)%meta(2), obs(iof)%meta(3))
    case (3) !H08 
      call get_nobs_H08(trim(OBS_IN_NAME(iof)),obs(iof)%nobs) ! H08
    case default
      write(6,*) 'Error: Unsupported observation file format!'
      stop
    end select

    write(6,'(5A,I9,A)') 'OBS FILE [', trim(OBS_IN_NAME(iof)), '] (FORMAT ', &
                         trim(obsformat_name(OBS_IN_FORMAT(iof))), '): TOTAL ', &
                         obs(iof)%nobs, ' OBSERVATIONS'

    call obs_info_allocate(obs(iof))

    select case (OBS_IN_FORMAT(iof))
    case (1)
      call read_obs(trim(OBS_IN_NAME(iof)),obs(iof))
    case (2)
      call read_obs_radar(trim(OBS_IN_NAME(iof)),obs(iof))
    case (3) ! H08 
      call read_obs_H08(trim(OBS_IN_NAME(iof)),obs(iof)) ! H08
    end select
  end do ! [ iof = 1, OBS_IN_NUM ]

  return
end subroutine read_obs_all

subroutine write_obs_all(obs, missing, file_suffix)
  implicit none

  type(obs_info), intent(in) :: obs(OBS_IN_NUM)
  logical, intent(in), optional :: missing
  character(len=*), intent(in), optional :: file_suffix
  logical :: missing_
  integer :: iof, strlen1, strlen2
  character(200) :: filestr

  missing_ = .true.
  IF(present(missing)) missing_ = missing

  do iof = 1, OBS_IN_NUM
    if (present(file_suffix)) then
      strlen1 = len(trim(OBS_IN_NAME(iof)))
      strlen2 = len(trim(file_suffix))
      write (filestr(1:strlen1),'(A)') trim(OBS_IN_NAME(iof))
      write (filestr(strlen1+1:strlen1+strlen2),'(A)') trim(file_suffix)
    else
      filestr = OBS_IN_NAME(iof)
    end if
    select case (OBS_IN_FORMAT(iof))
    case (1)
      call write_obs(trim(filestr),obs(iof),missing=missing_)
    case (2)
      call write_obs_radar(trim(filestr),obs(iof),missing=missing_)
    case (3) ! H08 
      call write_obs_H08(trim(filestr),obs(iof),missing=missing_) ! H08
    end select
  end do ! [ iof = 1, OBS_IN_NUM ]

  return
end subroutine write_obs_all
!
!-----------------------------------------------------------------------
!   TC vital obs subroutines by T. Honda (03/28/2016)
!-----------------------------------------------------------------------
!
SUBROUTINE search_tc_subdom(ritc,rjtc,v2d,yobs_tcx,yobs_tcy,yobs_mslp)
  use scale_grid, only: &
      GRID_CXG, &
      GRID_CYG, &
      DX, &
      DY
  use scale_grid_index, only: &
      IHALO, JHALO
  use scale_process, only: &
      PRC_myrank

  IMPLICIT NONE
  INTEGER :: il, jl, ig, jg
  REAL(r_size) :: xdis, ydis, rdis
  REAL(r_size),INTENT(IN) :: ritc, rjtc
  REAL(r_size),INTENT(IN) :: v2d(nlonh,nlath,nv2dd)
  REAL(r_size),INTENT(OUT) :: yobs_mslp !(Pa)
!  REAL(r_size),INTENT(OUT) :: yobs_lon, yobs_lat !(deg)
  REAL(r_size),INTENT(OUT) :: yobs_tcx, yobs_tcy !(m)

  REAL(r_size) :: slp2d(nlonh,nlath)
  REAL(r_size) :: dz, t, q, var5

  yobs_mslp = 9.99d33
  yobs_tcx = 9.99d33
  yobs_tcy = 9.99d33

  DO jl = 1, nlat 
  DO il = 1, nlon 
    t = v2d(il,jl,iv2dd_t2m)
    q = v2d(il,jl,iv2dd_q2m)
    dz = -1.0d0 * v2d(il,jl,iv2dd_topo)
    slp2d(il,jl) = v2d(il,jl,iv2dd_ps)
    call prsadj(slp2d(il,jl),dz,t,q)
  ENDDO
  ENDDO

  DO jl = JHALO + 1, nlat - JHALO
  DO il = IHALO + 1, nlon - IHALO
    call ij_l2g(PRC_myrank, il, jl, ig, jg)
    xdis = abs(real(ig,kind=r_size) - ritc) * DX
    ydis = abs(real(jg,kind=r_size) - rjtc) * DY
    rdis = sqrt(xdis*xdis + ydis*ydis)

    IF(rdis > TC_SEARCH_DIS)CYCLE

    IF(IHALO >= 2 .and. JHALO >= 2)THEN
      call wgt_ave2d(slp2d(:,:),il,jl,var5)
    ELSE
      var5 = slp2d(il,jl)
    ENDIF

    if(var5 < yobs_mslp)then
      yobs_mslp = var5
      yobs_tcx = (real(ig,kind=r_size) - 1.0d0) * DX + GRID_CXG(1)
      yobs_tcy = (real(jg,kind=r_size) - 1.0d0) * DY + GRID_CYG(1)
    endif
  ENDDO
  ENDDO

  RETURN
END SUBROUTINE search_tc_subdom

!-- 25 points weighted average (tentative)--
! 2D weight is...
!     1 1 1 1 1
!     1 3 3 3 1
!     1 3 5 3 1
!     1 3 3 3 1
!     1 1 1 1 1
!
SUBROUTINE wgt_ave2d(var,i,j,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlonh,nlath)
  INTEGER,INTENT(IN) :: i,j
  REAL(r_size),INTENT(OUT) :: var5

  var5 = ((var(i,j) * 5.0d0 + &
         (sum(var(i-1:i+1,j-1:j+1)) - var(i,j)) * 3.0d0 + &
         (sum(var(i-2:i+2,j-2:j+2)) - sum(var(i-1:i+1,j-1:j+1))) * 1.0d0)) / 45.0d0

  RETURN
END SUBROUTINE wgt_ave2d


!
!-----------------------------------------------------------------------
!   Himawari-8 obs subroutines by T. Honda (10/29/2015)
!-----------------------------------------------------------------------
#ifdef H08
!
SUBROUTINE Trans_XtoY_H08(nprof,ri,rj,lon,lat,v3d,v2d,yobs,plev_obs,qc,stggrd,yobs_H08_clr)
  use scale_mapproj, only: &
      MPRJ_rotcoef
  use scale_H08_fwd
  use scale_grid_index, only: &
    KHALO

  IMPLICIT NONE
  INTEGER :: n, np, k, ch
  INTEGER,INTENT(IN) :: nprof ! Num of Brightness Temp "Loc" observed by Himawari-8
                              ! NOTE: multiple channels (obs) on each grid point !!
  REAL(r_size),INTENT(IN) :: ri(nprof),rj(nprof)
  REAL(r_size),INTENT(IN) :: lon(nprof),lat(nprof)
  REAL(r_size),INTENT(IN) :: v3d(nlevh,nlonh,nlath,nv3dd)
  REAL(r_size),INTENT(IN) :: v2d(nlonh,nlath,nv2dd)
  INTEGER,INTENT(IN),OPTIONAL :: stggrd
  REAL(RP) :: rotc(2)

  INTEGER :: stggrd_ = 0

! -- 2D (nlevh,nbtobs) or 1D (nbtobs) profiles for RTTOV --  
  REAL(r_size) :: prs2d(nlevh,nprof)
  REAL(r_size) :: tk2d(nlevh,nprof)
  REAL(r_size) :: qv2d(nlevh,nprof)
  REAL(r_size) :: qliq2d(nlevh,nprof)
  REAL(r_size) :: qice2d(nlevh,nprof)

  REAL(r_size) :: tsfc1d(nprof)
  REAL(r_size) :: qsfc1d(nprof)
  REAL(r_size) :: psfc1d(nprof)
  REAL(r_size) :: usfc1d(nprof)
  REAL(r_size) :: vsfc1d(nprof)
  REAL(r_size) :: lon1d(nprof)
  REAL(r_size) :: lat1d(nprof)
  REAL(r_size) :: topo1d(nprof)
  REAL(r_size) :: lsmask1d(nprof)

! -- brightness temp from RTTOV
  REAL(r_size) :: btall_out(nch,nprof) ! NOTE: RTTOV always calculates all (10) channels!!
  REAL(r_size) :: btclr_out(nch,nprof) ! NOTE: RTTOV always calculates all (10) channels!!
! -- transmittance from RTTOV
  REAL(r_size) :: trans_out(nlev,nch,nprof)
 
  REAL(r_size) :: max_weight(nch,nprof)
  REAL(r_size) :: tmp_weight

  REAL(r_size),INTENT(OUT) :: yobs(nprof*nch)
  REAL(r_size),OPTIONAL,INTENT(OUT) :: yobs_H08_clr(nprof*nch)
  INTEGER,INTENT(OUT) :: qc(nprof*nch)
  REAL(r_size),INTENT(OUT) :: plev_obs(nch*nprof)

  REAL(r_size) :: rdp ! delta p
  INTEGER :: slev, elev

  REAL(r_size) :: utmp, vtmp ! U10m & V10m tmp for rotation

  if (present(stggrd)) stggrd_ = stggrd

  yobs = undef
  qc = iqc_good

  lon1d(:) = lon(:)
  lat1d(:) = lat(:)

! -- make profile arrays for RTTOV --
  DO np = 1, nprof ! -- make profiles

    CALL itpl_2d(v2d(:,:,iv2dd_skint),ri(np),rj(np),tsfc1d(np)) ! T2 is better??
!    CALL itpl_2d(v2d(:,:,iv2dd_t2m),ri(np),rj(np),tsfc1d(np))
    CALL itpl_2d(v2d(:,:,iv2dd_q2m),ri(np),rj(np),qsfc1d(np))
    CALL itpl_2d(v2d(:,:,iv2dd_topo),ri(np),rj(np),topo1d(np))
    CALL itpl_2d(v2d(:,:,iv2dd_lsmask),ri(np),rj(np),lsmask1d(np))
    CALL itpl_2d(v2d(:,:,iv2dd_ps),ri(np),rj(np),psfc1d(np))
!    call prsadj(yobs,rk-topo,t,q)
!    if (abs(rk-topo) > PS_ADJUST_THRES) then
!      write (6,'(A,F6.1)') 'warning: PS observation height adjustment exceeds the threshold. dz=', abs(rk-topo)
!      qc = iqc_ps_ter
!    end if

    if (stggrd_ == 1) then
      CALL itpl_2d(v2d(:,:,iv2dd_u10m),ri(np)-0.5,rj(np),utmp)  !###### should modity itpl_3d to prevent '1.0' problem....??
      CALL itpl_2d(v2d(:,:,iv2dd_v10m),ri(np),rj(np)-0.5,vtmp)  !######
    else
      CALL itpl_2d(v2d(:,:,iv2dd_u10m),ri(np),rj(np),utmp)
      CALL itpl_2d(v2d(:,:,iv2dd_v10m),ri(np),rj(np),vtmp)
    end if
    call MPRJ_rotcoef(rotc,lon(np)*deg2rad,lat(np)*deg2rad)
    usfc1d(np) = utmp * rotc(1) - vtmp * rotc(2)
    vsfc1d(np) = utmp * rotc(2) + vtmp * rotc(1)

    CALL itpl_2d_column(v3d(:,:,:,iv3dd_p),ri(np),rj(np),prs2d(:,np))
    CALL itpl_2d_column(v3d(:,:,:,iv3dd_t),ri(np),rj(np),tk2d(:,np))
    CALL itpl_2d_column(v3d(:,:,:,iv3dd_q),ri(np),rj(np),qv2d(:,np))
    CALL itpl_2d_column(v3d(:,:,:,iv3dd_qc),ri(np),rj(np),qliq2d(:,np))
    CALL itpl_2d_column((v3d(:,:,:,iv3dd_qi) &
                       + v3d(:,:,:,iv3dd_qs) &
                       + v3d(:,:,:,iv3dd_qg)),ri(np),rj(np),qice2d(:,np))

  ENDDO ! -- make profiles


!
! -- NOTE: The channel number for RTTOV is always 10, because it should be the same
!          with that in Himawari-8 RTTOV coef files.
!
!        : Satellite zenith angles are computed within SCALE_RTTOV_fwd using (lon,lat).
!

  slev = 1 + KHALO
  elev = nlevh - KHALO

  CALL SCALE_RTTOV_fwd(nch, & ! num of channels
                       nlev,& ! num of levels
                       nprof,& ! num of profs
                       prs2d(elev:slev:-1,1:nprof),& ! (Pa)
                       tk2d(elev:slev:-1,1:nprof),& ! (K)
                       qv2d(elev:slev:-1,1:nprof),& ! (kg/kg)
                       qliq2d(elev:slev:-1,1:nprof),& ! (kg/kg)
                       qice2d(elev:slev:-1,1:nprof),& ! (kg/kg)
                       tsfc1d(1:nprof),& ! (K)
                       qsfc1d(1:nprof),& ! (kg/kg)
                       psfc1d(1:nprof),& ! (Pa)
                       usfc1d(1:nprof),& ! (m/s)
                       vsfc1d(1:nprof),& ! (m/s)
                       topo1d(1:nprof),& ! (m)
                       lon1d(1:nprof),& ! (deg)
                       lat1d(1:nprof),& ! (deg)
                       lsmask1d(1:nprof),& ! (0-1)
                       btall_out(1:nch,1:nprof),& ! (K)
                       btclr_out(1:nch,1:nprof),& ! (K)
                       trans_out(nlev:1:-1,1:nch,1:nprof))
!
! -- Compute max weight level using trans_out 
! -- (Transmittance from each user pressure level to Top Of the Atmosphere)
! -- btall_out is substituted into yobs

  n = 0
  DO np = 1, nprof
  DO ch = 1, nch
    n = n + 1

    rdp = 1.0d0 / (prs2d(slev+1,np) - prs2d(slev,np))
    max_weight(ch,np) = abs((trans_out(2,ch,np) - trans_out(1,ch,np)) * rdp )


    plev_obs(n) = (prs2d(slev+1,np) + prs2d(slev,np)) * 0.5d0 ! (Pa)

    DO k = 2, nlev-1

      rdp = 1.0d0 / abs(prs2d(slev+k,np) - prs2d(slev+k-1,np))
      tmp_weight = (trans_out(k+1,ch,np) - trans_out(k,ch,np)) * rdp 
      if(tmp_weight > max_weight(ch,np))then
        max_weight(ch,np) = tmp_weight
        plev_obs(n) = (prs2d(slev+k,np) + prs2d(slev+k-1,np)) * 0.5d0 ! (Pa)
      endif
    ENDDO

    yobs(n) = btall_out(ch,np)
!
! ## comment out by T.Honda (02/09/2016)
! -- tentative QC here --
!    IF(plev_obs(n) >= H08_LIMIT_LEV)THEN
!      qc(n) = iqc_good
!    ELSE
!      qc(n) = iqc_obs_bad
!    ENDIF

    SELECT CASE(H08_CH_USE(ch))
    CASE(1)
      qc(n) = iqc_good
    CASE DEFAULT
      qc(n) = iqc_obs_bad
    END SELECT

    IF(H08_REJECT_LAND .and. (lsmask1d(np) > 0.5d0))THEN
      qc(n) = iqc_obs_bad
    ENDIF

    IF(abs(btall_out(ch,np) - btclr_out(ch,np)) > H08_CLDSKY_THRS)THEN
! Cloudy sky
      yobs(n) = yobs(n) * (-1.0d0)
    ELSE
! Clear sky
      yobs(n) = yobs(n) * 1.0d0
    ENDIF
   
    yobs_H08_clr(n) = btclr_out(ch,np)

  ENDDO ! ch
  ENDDO ! np

  RETURN
END SUBROUTINE Trans_XtoY_H08
#endif

SUBROUTINE get_nobs_H08(cfile,nn)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(OUT) :: nn ! num of all H08 obs
  REAL(r_sngl) :: wk(4+nch)
  INTEGER :: ios 
  INTEGER :: iprof
  INTEGER :: iunit
  LOGICAL :: ex

  nn = 0 
  iprof = 0
  iunit=91
      
      
  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
    
!    READ(iunit,IOSTAT=ios)wk
!    IF(ios /= 0) THEN 
!      WRITE(6,'(2A)') cfile,': Reading error -- skipped'
!      RETURN
!    END IF
    
    DO
      READ(iunit,IOSTAT=ios) wk
      IF(ios /= 0) EXIT
      iprof = iprof + 1
      nn = nn + nch
    END DO
    WRITE(6,*)' H08 FILE ', cfile
    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
    WRITE(6,'(A12,I10)') '   num of prof:',iprof
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF

  RETURN
END SUBROUTINE get_nobs_H08

SUBROUTINE read_obs_H08(cfile,obs)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(INOUT) :: obs
  REAL(r_sngl) :: wk(4+nch)
!  REAL(r_sngl) :: tmp
  INTEGER :: n,iunit

  INTEGER :: nprof, np, ch

  nprof = obs%nobs / nch
!  call obs_info_allocate(obs)

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')

  n = 0
  DO np=1,nprof
    READ(iunit) wk

    DO ch = 1, nch
      n = n + 1

      obs%elm(n) = NINT(wk(1))
      obs%typ(n) = NINT(wk(2))
      obs%lon(n) = REAL(wk(3),r_size)
      obs%lat(n) = REAL(wk(4),r_size)
      obs%dat(n) = REAL(wk(4+ch),r_size)
      obs%dif(n) = 0.0d0
      obs%lev(n) = ch + 6.0 ! substitute channnel number instead of the obs level
      obs%err(n) = REAL(OBSERR_H08(ch),r_size)
    END DO
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs_H08

SUBROUTINE write_obs_H08(cfile,obs,append,missing)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  TYPE(obs_info),INTENT(IN) :: obs
  LOGICAL,INTENT(IN),OPTIONAL :: append
  LOGICAL,INTENT(IN),OPTIONAL :: missing
  LOGICAL :: append_
  LOGICAL :: missing_
  REAL(r_sngl) :: wk(4+nch)
  INTEGER :: n,iunit
  INTEGER :: iprof, ns, ne

  iunit=92
  append_ = .false.
  IF(present(append)) append_ = append
  missing_ = .true.
  IF(present(missing)) missing_ = missing

  iprof = obs%nobs / nch


  IF(append_) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='append')
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  END IF

  DO n=1,iprof
    ns = (n-1)*nch + 1
    ne = n*nch

    wk(1) = REAL(obs%elm(ns),r_sngl)
    wk(2) = REAL(obs%typ(ns),r_sngl)
    wk(3) = REAL(obs%lon(ns),r_sngl)
    wk(4) = REAL(obs%lat(ns),r_sngl)
    wk(5:5+nch-1) = REAL(obs%dat(ns:ne),r_size)
    WRITE(iunit) wk

  ENDDO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs_H08

END MODULE common_obs_scale
