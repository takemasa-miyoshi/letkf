PROGRAM obsmake
!=======================================================================
!
! [PURPOSE:] For generating observations
!
! [HISTORY:]
!   01/26/2009 Takemasa Miyoshi  created
!   05/06/2009 Takemasa Miyoshi  modified for station.tbl format change
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_roms
  USE common_obs_roms

  IMPLICIT NONE

  LOGICAL,PARAMETER :: verbose=.FALSE.
  INTEGER,PARAMETER :: id_hrf   =1
  INTEGER,PARAMETER :: id_prof  =2
  INTEGER,PARAMETER :: id_sst   =3
  INTEGER,PARAMETER :: kprof = 17
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: error(:)
  REAL(r_sngl) :: wk(6)
  REAL(r_size) :: alon,depth(nlev),surfd
  REAL(r_size) :: profdep(kprof)
  REAL(r_size) :: rk,wk1,wk2
  REAL(r_size) :: wk1prof(kprof),wk2prof(kprof)
  INTEGER :: id,i,j,k,kk,elm,nn,ios
  INTEGER :: nobkind,nu,nv,nt,ns,nz,obelm(5)
  REAL(r_size) :: a,oberr(5)
  CHARACTER(100) :: cdummy
  CHARACTER(7) :: truef='true.nc'
  LOGICAL :: ex
  !
  ! Initial settings
  !
  CALL set_common_roms
  profdep=(/ 0.0, -5.0, -10.0, -15.0, -20.0, -30.0, -40.0, -50.0, -60.0, -75.0, -100.0, -125.0, -150.0, -200.0, -250.0, -300.0, -400.0 /)
  !
  ! read table obserr.tbl
  !
  OPEN(11,FILE='obserr.tbl')
  READ(11,'(A)') cdummy
  READ(11,'(A)') cdummy
  nobkind = 0
  DO
    READ(11,'(L1,X,I5,X,ES9.2)',IOSTAT=ios) ex,i,a
    IF(ios /= 0) EXIT
    IF(ex) THEN
      nobkind = nobkind+1
      obelm(nobkind) = i
      oberr(nobkind) = a
    END IF
  END DO
  CLOSE(11)
  !
  ! find n
  !
  DO nu=1,nobkind
    IF(obelm(nu) == id_u_obs) EXIT
  END DO
  DO nv=1,nobkind
    IF(obelm(nv) == id_v_obs) EXIT
  END DO
  DO nt=1,nobkind
    IF(obelm(nt) == id_t_obs) EXIT
  END DO
  DO ns=1,nobkind
    IF(obelm(ns) == id_s_obs) EXIT
  END DO
  DO nz=1,nobkind
    IF(obelm(nz) == id_z_obs) EXIT
  END DO
  !
  ! Read true (nature run)
  !
  ALLOCATE(v3d(nlon,nlat,nlev,nv3d))
  ALLOCATE(v2d(nlon,nlat,nv2d))
  CALL read_grd(truef,v3d,v2d)
  !
  ! Read station (station.tbl)
  !
  OPEN(10,FILE='station.tbl')
  READ(10,'(A)') cdummy
  READ(10,'(A)') cdummy
  nn = 0
  DO
    READ(10,'(3I4)',IOSTAT=ios) id,i,j
    IF(ios /= 0) EXIT
    IF(id == id_hrf) nn = nn + 2
    IF(id == id_prof) nn = nn + kprof*2
    IF(id == id_sst) nn = nn + 1
  END DO
  CLOSE(10)
  IF(verbose) PRINT '(A,I)','nn = ',nn
  !
  ! Random number
  !
  ALLOCATE(error(nn))
  CALL com_randn(nn,error)
!  error = 0.0d0 !! for testing purpose
  !
  ! Output OBS
  !
  nn = 0
  OPEN(10,FILE='station.tbl')
  READ(10,'(A)') cdummy
  READ(10,'(A)') cdummy
  ! hfradar_uv
  OPEN(91,FILE='hrf.letkf',FORM='unformatted',ACCESS='sequential')
  OPEN(81,FILE='hrf.ascii',FORM='formatted')
  WRITE(81,'(A)') '%   lon       lat       u        v       speed      Uerr      Verr'
  WRITE(81,'(A)') '%  (deg)     (deg)    (cm/s)   (cm/s)    (cm/s)    (cm/s)    (cm/s))'
  ! prof
  OPEN(92,FILE='prof.letkf',FORM='unformatted',ACCESS='sequential')
  OPEN(82,FILE='prof.ascii',FORM='formatted')
  ! ship_sst
  OPEN(93,FILE='sst.letkf',FORM='unformatted',ACCESS='sequential')
  OPEN(83,FILE='sst.ascii',FORM='formatted')
  DO
    READ(10,'(3I4)',IOSTAT=ios) id,i,j
    alon = lon(i,j)
    IF(alon > 180.0d0) alon = alon - 360.0d0
    IF(ios /= 0) EXIT
    SELECT CASE(id)
    !
    ! hfradar_uv
    !
    CASE(id_hrf)
      wk1=v3d(i,j,nlev,iv3d_u)+error(nn+1)*oberr(nu)
      wk2=v3d(i,j,nlev,iv3d_v)+error(nn+2)*oberr(nv)
      WRITE(81,'(7(F9.4X))') lat(i,j),alon,wk1*100.0,wk2*100.0,SQRT(wk1**2+wk2**2)*100.0,oberr(nu)*100.0,oberr(nv)*100.0
      wk(1)=REAL(id_u_obs,r_sngl)
      wk(2)=REAL(i,r_sngl)
      wk(3)=REAL(j,r_sngl)
      wk(4)=0.0d0
      wk(5)=REAL(wk1,r_sngl)
      wk(6)=REAL(oberr(nu),r_sngl)
      WRITE(91) wk
      wk(1)=REAL(id_v_obs,r_sngl)
      wk(5)=REAL(wk2,r_sngl)
      wk(6)=REAL(oberr(nv),r_sngl)
      WRITE(91) wk
      nn = nn+2
    !
    ! profile
    !
    CASE(id_prof)
      CALL calc_depth(v2d(i,j,iv2d_z),phi0(i,j),depth)
      CALL com_interp_spline(nlev,depth,v3d(i,j,:,iv3d_t),kprof,profdep,wk1prof)
      CALL com_interp_spline(nlev,depth,v3d(i,j,:,iv3d_s),kprof,profdep,wk2prof)
      DO k=1,kprof
        IF(k == 1) THEN
          wk1 =-999.0d0
          wk2 =-999.0d0
        ELSE
          wk1 = wk1prof(k) + error(nn+1)*oberr(nt)
          wk2 = wk2prof(k) + error(nn+2)*oberr(ns)
        END IF
        WRITE(82,'(5(X,F9.4))') alon,lat(i,j),profdep(k),wk1,wk2
        IF(k /= 1) THEN
          wk(1)=REAL(id_t_obs,r_sngl)
          wk(2)=REAL(i,r_sngl)
          wk(3)=REAL(j,r_sngl)
          wk(4)=profdep(k)
          wk(5)=REAL(wk1,r_sngl)
          wk(6)=REAL(oberr(nt),r_sngl)
          WRITE(92) wk
          wk(1)=REAL(id_s_obs,r_sngl)
          wk(5)=REAL(wk2,r_sngl)
          wk(6)=REAL(oberr(ns),r_sngl)
          WRITE(92) wk
        END IF
        nn = nn+2
      END DO
    !
    ! SST
    !
    CASE(id_sst)
      wk1=v3d(i,j,nlev,iv3d_t)+error(nn+1)*oberr(nt)
      WRITE(83,'(4(X,F9.4))') alon,lat(i,j),wk1,-999.0
      wk(1)=REAL(id_t_obs,r_sngl)
      wk(2)=REAL(i,r_sngl)
      wk(3)=REAL(j,r_sngl)
      wk(4)=0.0d0
      wk(5)=REAL(wk1,r_sngl)
      wk(6)=REAL(oberr(nt),r_sngl)
      WRITE(93) wk
      nn = nn+1
    END SELECT
  END DO
  CLOSE(81)
  CLOSE(82)
  CLOSE(83)
  CLOSE(91)
  CLOSE(92)
  CLOSE(93)

  DEALLOCATE(v3d,v2d,error)

  STOP
END PROGRAM obsmake
