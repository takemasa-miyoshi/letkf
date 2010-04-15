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
  USE common_obs
  USE common_roms

  IMPLICIT NONE

  LOGICAL,PARAMETER :: verbose=.FALSE.
  INTEGER,PARAMETER :: id_hrf   =1
  INTEGER,PARAMETER :: id_prof  =2
  INTEGER,PARAMETER :: id_sst   =3
  INTEGER,PARAMETER :: id_ssh   =4
  INTEGER,PARAMETER :: id_glider=5
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: error(:)
  REAL(r_sngl) :: wk(6)
  REAL(r_sngl) :: wk2(6)
  REAL(r_size) :: alon,depth,cs(nlev)
  REAL(r_size) :: err
  INTEGER :: id,i,j,k,elm,nn,ios
  CHARACTER(100) :: cdummy

  CHARACTER(7) :: truef='true.nc'
  !
  ! Initial settings
  !
  CALL set_common_roms
  cs=(/ -0.921962804129276, -0.783682417800398, -0.666141723246888, -0.566230075492245, -0.481303371858989, -0.409114077354825, -0.347751744988103, -0.295592456917874, -0.251255848428496, -0.213568577393686, -0.181533272470966, -0.154302138258862, -0.131154518892926, -0.111477826314426, -0.0947513284896102, -0.0805323685443561, -0.0684446501117523, -0.0581682788712608, -0.0494312967346868, -0.0420024846361897, -0.0356852434564318, -0.0303123911431393, -0.0257417383369477, -0.0218523254140864, -0.0185412213614628, -0.0157207997682656, -0.0133164198448193, -0.0112644510982221, -0.00951058938897752, -0.00800841980405428, -0.0067181883136592, -0.00560574970434628, -0.00464166394613383, -0.00380041707951966, -0.00305974600362987, -0.00240004929641906, -0.00180386847462664, -0.00125542596534104, -0.000740207561862293, -0.000244578313804319 /)

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
  nobs = 0
  DO
    READ(10,'(3I4,I3,I6,ES9.2)',IOSTAT=ios) id,i,j,k,elm,err
    IF(ios /= 0) EXIT
    nobs = nobs + 1
  END DO
  CLOSE(10)
  IF(verbose) PRINT '(A,I)','nobs = ',nobs
  !
  ! Random number
  !
  ALLOCATE(error(nobs))
  CALL com_randn(nobs,error)
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
  ! ship_sst
  OPEN(92,FILE='sst.letkf',FORM='unformatted',ACCESS='sequential')
  OPEN(82,FILE='sst.ascii',FORM='formatted')
  ! ship_ssh
  OPEN(93,FILE='ssh.letkf',FORM='unformatted',ACCESS='sequential')
  OPEN(83,FILE='ssh.ascii',FORM='formatted')
  ! prof
  OPEN(94,FILE='prof.letkf',FORM='unformatted',ACCESS='sequential')
  OPEN(84,FILE='prof.ascii',FORM='formatted')
  ! prof uv
  OPEN(95,FILE='profuv.letkf',FORM='unformatted',ACCESS='sequential')
  OPEN(85,FILE='profuv.ascii',FORM='formatted')
  DO
    READ(10,'(3I4,I3,I6,ES9.2)',IOSTAT=ios) id,i,j,k,elm,err
    IF(ios /= 0) EXIT
    !
    ! LETKF obs format
    !
    wk(1)=REAL(elm,r_sngl)
    wk(2)=REAL(i,r_sngl)
    wk(3)=REAL(j,r_sngl)
    wk(4)=REAL(k,r_sngl)
    nn = nn+1
    SELECT CASE(elm)
    CASE(id_u_obs)
      wk(5)=REAL(v3d(i,j,k,iv3d_u)+error(nn)*err,r_sngl)
    CASE(id_v_obs)
      wk(5)=REAL(v3d(i,j,k,iv3d_v)+error(nn)*err,r_sngl)
    CASE(id_t_obs)
      wk(5)=REAL(v3d(i,j,k,iv3d_t)+error(nn)*err,r_sngl)
    CASE(id_s_obs)
      wk(5)=REAL(v3d(i,j,k,iv3d_s)+error(nn)*err,r_sngl)
    CASE(id_z_obs)
      wk(5)=REAL(v2d(i,j,iv2d_z)+error(nn)*err,r_sngl)
    END SELECT
    wk(6)=REAL(err,r_sngl)
    !
    ! real obs format
    !
    alon = lon(i,j)
    IF(alon > 180.0d0) alon = alon - 360.0d0
    !! hfradar_uv
    IF(id == id_hrf .AND. elm == id_u_obs) wk2 = wk
    IF(id == id_hrf .AND. elm == id_v_obs) THEN
      WRITE(81,'(7(F9.4X))') lat(i,j),alon,wk2(5)*100.0,wk(5)*100.0,SQRT(wk(5)**2+wk2(5)**2)*100.0,wk2(6)*100.0,wk(6)*100.0
      WRITE(91) wk2
      WRITE(91) wk
    END IF
    !! sst
    IF(id == id_sst) THEN
      WRITE(82,'(4(X,F9.4))') alon,lat(i,j),wk(5),-999.0
      WRITE(92) wk
    END IF
    !! ssh
    IF(id == id_ssh) THEN
      WRITE(83,'(4(X,F9.4))') alon,lat(i,j),wk(5),-999.0
      WRITE(93) wk
    END IF
    !! profile
    IF(id == id_prof .AND. elm == id_t_obs) wk2 = wk
    IF(id == id_prof .AND. elm == id_s_obs) THEN
      depth = v2d(i,j,iv2d_z)+(v2d(i,j,iv2d_z)+3500.0d0)* &
        & (10.0d0*(REAL(k,r_size)-REAL(nlev,r_size)-0.5d0)/ &
        & REAL(nlev,r_size)+cs(k)*3500.0d0)/3510.0d0
      WRITE(84,'(5(X,F9.4))') alon,lat(i,j),depth,wk2(5),wk(5)
      WRITE(94) wk2
      WRITE(94) wk
    END IF
    !! profile uv
    IF(id == id_prof .AND. elm == id_u_obs) wk2 = wk
    IF(id == id_prof .AND. elm == id_v_obs) THEN
      depth = v2d(i,j,iv2d_z)+(v2d(i,j,iv2d_z)+3500.0d0)* &
        & (10.0d0*(REAL(k,r_size)-REAL(nlev,r_size)-0.5d0)/ &
        & REAL(nlev,r_size)+cs(k)*3500.0d0)/3510.0d0
      WRITE(85,'(5(X,F9.4))') alon,lat(i,j),depth,wk2(5)*100.0,wk(5)*100.0
      WRITE(95) wk2
      WRITE(95) wk
    END IF
  END DO
  CLOSE(91)
  CLOSE(92)
  CLOSE(93)
  CLOSE(94)
  CLOSE(95)
  CLOSE(81)
  CLOSE(82)
  CLOSE(83)
  CLOSE(84)
  CLOSE(85)

  DEALLOCATE(v3d,v2d,error)

  STOP
END PROGRAM obsmake
