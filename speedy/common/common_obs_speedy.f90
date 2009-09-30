MODULE common_obs_speedy
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_obs
  USE common_speedy
  USE common_mpi_speedy
  USE common_letkf

  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: nslots=1 ! number of time slots for 4D-LETKF
  INTEGER,PARAMETER :: nbslot=1 ! basetime slot
  REAL(r_size),PARAMETER :: sigma_obs=700.0d3
!  REAL(r_size),PARAMETER :: sigma_obsv=0.4d0
  REAL(r_size),PARAMETER :: sigma_obsv=0.1d0
  REAL(r_size),PARAMETER :: sigma_obst=3.0d0
  REAL(r_size),SAVE :: dist_zero
  REAL(r_size),SAVE :: dist_zerov
  REAL(r_size),ALLOCATABLE,SAVE :: dlon_zero(:)
  REAL(r_size),SAVE :: dlat_zero
  REAL(r_size),ALLOCATABLE,SAVE :: obselm(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslon(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslat(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslev(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsdat(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obserr(:)
!  REAL(r_size),ALLOCATABLE,SAVE :: obsk(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsdep(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obshdxf(:,:)
  INTEGER,SAVE :: nobsgrd(nlon,nlat)

CONTAINS
!-----------------------------------------------------------------------
! Initialize
!-----------------------------------------------------------------------
SUBROUTINE set_common_obs_speedy
  IMPLICIT NONE
  REAL(r_size) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size) :: v2d(nlon,nlat,nv2d)
  REAL(r_size) :: p_full(nlon,nlat,nlev)
  REAL(r_size),PARAMETER :: threshold_dz=1000.0d0
  REAL(r_size),PARAMETER :: gross_error=10.0d0
  REAL(r_size) :: dz,tg,qg
  REAL(r_size) :: ri,rj
  REAL(r_size) :: dlon1,dlon2,dlon,dlat
  REAL(r_size),ALLOCATABLE :: wk2d(:,:)
  INTEGER,ALLOCATABLE :: iwk2d(:,:)
  REAL(r_size),ALLOCATABLE :: tmpelm(:)
  REAL(r_size),ALLOCATABLE :: tmplon(:)
  REAL(r_size),ALLOCATABLE :: tmplat(:)
  REAL(r_size),ALLOCATABLE :: tmplev(:)
  REAL(r_size),ALLOCATABLE :: tmpdat(:)
  REAL(r_size),ALLOCATABLE :: tmperr(:)
  REAL(r_size),ALLOCATABLE :: tmpk(:)
  REAL(r_size),ALLOCATABLE :: tmpdep(:)
  REAL(r_size),ALLOCATABLE :: tmphdxf(:,:)
  INTEGER,ALLOCATABLE :: tmpqc0(:,:)
  INTEGER,ALLOCATABLE :: tmpqc(:)
  REAL(r_size),ALLOCATABLE :: tmp2elm(:)
  REAL(r_size),ALLOCATABLE :: tmp2lon(:)
  REAL(r_size),ALLOCATABLE :: tmp2lat(:)
  REAL(r_size),ALLOCATABLE :: tmp2lev(:)
  REAL(r_size),ALLOCATABLE :: tmp2dat(:)
  REAL(r_size),ALLOCATABLE :: tmp2err(:)
!  REAL(r_size),ALLOCATABLE :: tmp2k(:)
  REAL(r_size),ALLOCATABLE :: tmp2dep(:)
  REAL(r_size),ALLOCATABLE :: tmp2hdxf(:,:)
  INTEGER :: nobslots(nslots)
  INTEGER :: n,i,j,ierr,islot,nn,l,im
  INTEGER :: nj(0:nlat-1)
  INTEGER :: njs(1:nlat-1)
  CHARACTER(9) :: obsfile='obsTT.dat'
  CHARACTER(11) :: guesfile='gsTTNNN.grd'

  WRITE(6,'(A)') 'Hello from set_common_obs_speedy'

  dist_zero = sigma_obs * SQRT(10.0d0/3.0d0) * 2.0d0
  dist_zerov = sigma_obsv * SQRT(10.0d0/3.0d0) * 2.0d0
  dlat_zero = dist_zero / pi / re * 180.0d0
  ALLOCATE(dlon_zero(nij1))
  DO i=1,nij1
    dlon_zero(i) = dlat_zero / COS(pi*lat1(i)/180.0d0)
  END DO

  DO islot=1,nslots
    WRITE(obsfile(4:5),'(I2.2)') islot
    CALL get_nobs_mpi(obsfile,nobslots(islot))
  END DO
  nobs = SUM(nobslots)
  WRITE(6,'(I10,A)') nobs,' TOTAL OBSERVATIONS INPUT'
!
! INITIALIZE GLOBAL VARIABLES
!
  ALLOCATE( tmpelm(nobs) )
  ALLOCATE( tmplon(nobs) )
  ALLOCATE( tmplat(nobs) )
  ALLOCATE( tmplev(nobs) )
  ALLOCATE( tmpdat(nobs) )
  ALLOCATE( tmperr(nobs) )
  ALLOCATE( tmpk(nobs) )
  ALLOCATE( tmpdep(nobs) )
  ALLOCATE( tmphdxf(nobs,nbv) )
  ALLOCATE( tmpqc0(nobs,nbv) )
  ALLOCATE( tmpqc(nobs) )
  tmpqc0 = 0
  tmphdxf = 0.0d0
!
! LOOP of timeslots
!
  nn=0
  timeslots: DO islot=1,nslots
    IF(nobslots(islot) == 0) CYCLE
    WRITE(obsfile(4:5),'(I2.2)') islot
    CALL read_obs_mpi(obsfile,nobslots(islot),&
      & tmpelm(nn+1:nn+nobslots(islot)),tmplon(nn+1:nn+nobslots(islot)),&
      & tmplat(nn+1:nn+nobslots(islot)),tmplev(nn+1:nn+nobslots(islot)),&
      & tmpdat(nn+1:nn+nobslots(islot)),tmperr(nn+1:nn+nobslots(islot)) )
    l=0
    DO
      im = myrank+1 + nprocs * l
      IF(im > nbv) EXIT
      WRITE(guesfile(3:7),'(I2.2,I3.3)') islot,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',guesfile
      CALL read_grd(guesfile,v3d,v2d)
      CALL calc_pfull(nlon,nlat,v2d(:,:,iv2d_ps),p_full)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,ri,rj,dz,tg,qg)
      DO n=1,nobslots(islot)
        CALL phys2ijk(p_full,tmpelm(nn+n),&
          & tmplon(nn+n),tmplat(nn+n),tmplev(nn+n),ri,rj,tmpk(nn+n))
        IF(CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri)) THEN
!$OMP CRITICAL
          WRITE(6,'(A)') '* X-coordinate out of range'
          WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', rlon=',tmplon(nn+n)
!$OMP END CRITICAL
          CYCLE
        END IF
        IF(CEILING(rj) < 2 .OR. nlat < CEILING(rj)) THEN
!$OMP CRITICAL
          WRITE(6,'(A)') '* Y-coordinate out of range'
          WRITE(6,'(A,F6.2,A,F6.2)') '*   rj=',rj,', rlat=',tmplat(nn+n)
!$OMP END CRITICAL
          CYCLE
        END IF
        IF(CEILING(tmpk(n+nn)) > nlev) THEN
          CALL itpl_2d(phi0,ri,rj,dz)
!$OMP CRITICAL
          WRITE(6,'(A)') '* Z-coordinate out of range'
          WRITE(6,'(A,F6.2,A,F10.2,A,F6.2,A,F6.2,A,F10.2)') &
           & '*   rk=',tmpk(nn+n),', rlev=',tmplev(nn+n),&
           & ', (lon,lat)=(',tmplon(nn+n),',',tmplat(nn+n),'), phi0=',dz
!$OMP END CRITICAL
          CYCLE
        END IF
        IF(CEILING(tmpk(nn+n)) < 2 .AND. NINT(tmpelm(nn+n)) /= id_ps_obs) THEN
          IF(NINT(tmpelm(nn+n)) == id_u_obs .OR.&
           & NINT(tmpelm(nn+n)) == id_v_obs) THEN
            tmpk(nn+n) = 1.00001d0
          ELSE
            CALL itpl_2d(phi0,ri,rj,dz)
!$OMP CRITICAL
            WRITE(6,'(A)') '* Z-coordinate out of range'
            WRITE(6,'(A,F6.2,A,F10.2,A,F6.2,A,F6.2,A,F10.2)') &
             & '*   rk=',tmpk(nn+n),', rlev=',tmplev(nn+n),&
             & ', (lon,lat)=(',tmplon(nn+n),',',tmplat(nn+n),'), phi0=',dz
!$OMP END CRITICAL
            CYCLE
          END IF
        END IF
        IF(NINT(tmpelm(nn+n)) == id_ps_obs .AND. tmpdat(nn+n) < -100.0d0) THEN
          CYCLE
        END IF
!       IF(NINT(tmpelm(nn+n)) == id_ps_obs) THEN
!         CALL itpl_2d(phi0,ri,rj,dz)
!         dz = dz - tmplev(nn+n)
!         IF(ABS(dz) < threshold_dz) THEN ! pressure adjustment threshold
!           CALL itpl_2d(t(:,:,1),ri,rj,tg)
!           CALL itpl_2d(q(:,:,1),ri,rj,qg)
!           CALL prsadj(tmpdat(nn+n),dz,tg,qg)
!         ELSE
!!OMP CRITICAL
!           PRINT '(A)','PS obs vertical adjustment beyond threshold'
!           PRINT '(A,F10.2,A,F6.2,A,F6.2,A)',&
!             & '  dz=',dz,', (lon,lat)=(',tmplon(nn+n),',',tmplat(nn+n),')'
!!OMP END CRITICAL
!           CYCLE
!         END IF
!       END IF
        !
        ! observational operator
        !
        CALL Trans_XtoY(tmpelm(nn+n),&
          & ri,rj,tmpk(nn+n),v3d,v2d,p_full,tmphdxf(nn+n,im))
        tmpqc0(nn+n,im) = 1
      END DO
!$OMP END PARALLEL DO
      l = l+1
    END DO
    nn = nn + nobslots(islot)
  END DO timeslots

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ALLOCATE(wk2d(nobs,nbv))
  wk2d = tmphdxf
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(wk2d,tmphdxf,nobs*nbv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  DEALLOCATE(wk2d)
  ALLOCATE(iwk2d(nobs,nbv))
  iwk2d = tmpqc0
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(iwk2d,tmpqc0,nobs*nbv,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  DEALLOCATE(iwk2d)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i)
  DO n=1,nobs
    tmpqc(n) = MINVAL(tmpqc0(n,:))
    IF(tmpqc(n) /= 1) CYCLE
    tmpdep(n) = tmphdxf(n,1)
    DO i=2,nbv
      tmpdep(n) = tmpdep(n) + tmphdxf(n,i)
    END DO
    tmpdep(n) = tmpdep(n) / REAL(nbv,r_size)
    DO i=1,nbv
      tmphdxf(n,i) = tmphdxf(n,i) - tmpdep(n) ! Hdx
    END DO
    tmpdep(n) = tmpdat(n) - tmpdep(n) ! y-Hx
    IF(ABS(tmpdep(n)) > gross_error*tmperr(n)) THEN !gross error
      tmpqc(n) = 0
    END IF
  END DO
!$OMP END PARALLEL DO
  DEALLOCATE(tmpqc0)

  WRITE(6,'(I10,A)') SUM(tmpqc),' OBSERVATIONS TO BE ASSIMILATED'

  CALL monit_dep(nobs,tmpelm,tmpdep,tmpqc)
!
! temporal observation localization
!
!  nn = 0
!  DO islot=1,nslots
!    tmperr(nn+1:nn+nobslots(islot)) = tmperr(nn+1:nn+nobslots(islot)) &
!      & * exp(0.25d0 * (REAL(islot-nbslot,r_size) / sigma_obst)**2)
!    nn = nn + nobslots(islot)
!  END DO
!
! SELECT OBS IN THE NODE
!
  nn = 0
  DO n=1,nobs
    IF(tmpqc(n) /= 1) CYCLE
!    IF(tmplat(n) < MINVAL(lat1) .OR. MAXVAL(lat1) < tmplat(n)) THEN
!      dlat = MIN( ABS(MINVAL(lat1)-tmplat(n)),ABS(MAXVAL(lat1)-tmplat(n)) )
!      IF(dlat > dlat_zero) CYCLE
!    END IF
!    IF(tmplon(n) < MINVAL(lon1) .OR. MAXVAL(lon1) < tmplon(n)) THEN
!      dlon1 = ABS(MINVAL(lon1) - tmplon(n))
!      dlon1 = MIN(dlon1,360.0d0-dlon1)
!      dlon2 = ABS(MAXVAL(lon1) - tmplon(n))
!      dlon2 = MIN(dlon2,360.0d0-dlon2)
!      dlon =  MIN(dlon1,dlon2) &
!         & * pi*re*COS(tmplat(n)*pi/180.d0)/180.0d0
!      IF(dlon > dist_zero) CYCLE
!    END IF
    nn = nn+1
    tmpelm(nn) = tmpelm(n)
    tmplon(nn) = tmplon(n)
    tmplat(nn) = tmplat(n)
    tmplev(nn) = tmplev(n)
    tmpdat(nn) = tmpdat(n)
    tmperr(nn) = tmperr(n)
    tmpk(nn) = tmpk(n)
    tmpdep(nn) = tmpdep(n)
    tmphdxf(nn,:) = tmphdxf(n,:)
    tmpqc(nn) = tmpqc(n)
  END DO
  nobs = nn
  WRITE(6,'(I10,A,I3.3)') nobs,' OBSERVATIONS TO BE ASSIMILATED IN MYRANK ',myrank
!
! SORT
!
  ALLOCATE( tmp2elm(nobs) )
  ALLOCATE( tmp2lon(nobs) )
  ALLOCATE( tmp2lat(nobs) )
  ALLOCATE( tmp2lev(nobs) )
  ALLOCATE( tmp2dat(nobs) )
  ALLOCATE( tmp2err(nobs) )
!  ALLOCATE( tmp2k(nobs) )
  ALLOCATE( tmp2dep(nobs) )
  ALLOCATE( tmp2hdxf(nobs,nbv) )
  ALLOCATE( obselm(nobs) )
  ALLOCATE( obslon(nobs) )
  ALLOCATE( obslat(nobs) )
  ALLOCATE( obslev(nobs) )
  ALLOCATE( obsdat(nobs) )
  ALLOCATE( obserr(nobs) )
!  ALLOCATE( obsk(nobs) )
  ALLOCATE( obsdep(nobs) )
  ALLOCATE( obshdxf(nobs,nbv) )
  nobsgrd = 0
  nj = 0
!$OMP PARALLEL PRIVATE(i,j,n,nn)
!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    DO n=1,nobs
      IF(tmplat(n) < lat(j) .OR. lat(j+1) <= tmplat(n)) CYCLE
      nj(j) = nj(j) + 1
    END DO
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    njs(j) = SUM(nj(0:j-1))
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    nn = 0
    DO n=1,nobs
      IF(tmplat(n) < lat(j) .OR. lat(j+1) <= tmplat(n)) CYCLE
      nn = nn + 1
      tmp2elm(njs(j)+nn) = tmpelm(n)
      tmp2lon(njs(j)+nn) = tmplon(n)
      tmp2lat(njs(j)+nn) = tmplat(n)
      tmp2lev(njs(j)+nn) = tmplev(n)
      tmp2dat(njs(j)+nn) = tmpdat(n)
      tmp2err(njs(j)+nn) = tmperr(n)
!      tmp2k(njs(j)+nn) = tmpk(n)
      tmp2dep(njs(j)+nn) = tmpdep(n)
      tmp2hdxf(njs(j)+nn,:) = tmphdxf(n,:)
    END DO
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    IF(nj(j) == 0) THEN
      nobsgrd(:,j) = njs(j)
      CYCLE
    END IF
    nn = 0
    DO i=1,nlon
      DO n=njs(j)+1,njs(j)+nj(j)
        IF(i < nlon) THEN
          IF(tmp2lon(n) < lon(i) .OR. lon(i+1) <= tmp2lon(n)) CYCLE
        ELSE
          IF(tmp2lon(n) < lon(nlon) .OR. 360.0d0 <= tmp2lon(n)) CYCLE
        END IF
        nn = nn + 1
        obselm(njs(j)+nn) = tmp2elm(n)
        obslon(njs(j)+nn) = tmp2lon(n)
        obslat(njs(j)+nn) = tmp2lat(n)
        obslev(njs(j)+nn) = tmp2lev(n)
        obsdat(njs(j)+nn) = tmp2dat(n)
        obserr(njs(j)+nn) = tmp2err(n)
!        obsk(njs(j)+nn) = tmp2k(n)
        obsdep(njs(j)+nn) = tmp2dep(n)
        obshdxf(njs(j)+nn,:) = tmp2hdxf(n,:)
      END DO
      nobsgrd(i,j) = njs(j) + nn
    END DO
    IF(nn /= nj(j)) THEN
!$OMP CRITICAL
      WRITE(6,'(A,2I)') 'OBS DATA SORT ERROR: ',nn,nj(j)
      WRITE(6,'(F6.2,A,F6.2)') lat(j),'< LAT <',lat(j+1)
      WRITE(6,'(F6.2,A,F6.2)') MINVAL(tmp2lat(njs(j)+1:njs(j)+nj(j))),'< OBSLAT <',MAXVAL(tmp2lat(njs(j)+1:njs(j)+nj(j)))
!$OMP END CRITICAL
    END IF
  END DO
!$OMP END DO
!$OMP END PARALLEL
  DEALLOCATE( tmp2elm )
  DEALLOCATE( tmp2lon )
  DEALLOCATE( tmp2lat )
  DEALLOCATE( tmp2lev )
  DEALLOCATE( tmp2dat )
  DEALLOCATE( tmp2err )
!  DEALLOCATE( tmp2k )
  DEALLOCATE( tmp2dep )
  DEALLOCATE( tmp2hdxf )
  DEALLOCATE( tmpelm )
  DEALLOCATE( tmplon )
  DEALLOCATE( tmplat )
  DEALLOCATE( tmplev )
  DEALLOCATE( tmpdat )
  DEALLOCATE( tmperr )
  DEALLOCATE( tmpk )
  DEALLOCATE( tmpdep )
  DEALLOCATE( tmphdxf )
  DEALLOCATE( tmpqc )

  RETURN
END SUBROUTINE set_common_obs_speedy
!-----------------------------------------------------------------------
! Transformation from model variables to an observation
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY(elm,ri,rj,rk,v3d,v2d,p_full,yobs)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rk
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_size),INTENT(IN) :: p_full(nlon,nlat,nlev)
  REAL(r_size),INTENT(OUT) :: yobs
  REAL(r_size) :: rh(nlon,nlat,nlev)
  INTEGER :: i,j,k
  INTEGER :: is,ie,js,je,ks,ke
  ie = CEILING( ri )
  is = ie-1
  je = CEILING( rj )
  js = je-1
  ke = CEILING( rk )
  ks = ke-1

  SELECT CASE (NINT(elm))
  CASE(id_u_obs)  ! U
    CALL itpl_3d(v3d(:,:,:,iv3d_u),ri,rj,rk,yobs)
  CASE(id_v_obs)  ! V
    CALL itpl_3d(v3d(:,:,:,iv3d_v),ri,rj,rk,yobs)
  CASE(id_t_obs)  ! T
    CALL itpl_3d(v3d(:,:,:,iv3d_t),ri,rj,rk,yobs)
  CASE(id_q_obs)  ! Q
    CALL itpl_3d(v3d(:,:,:,iv3d_q),ri,rj,rk,yobs)
  CASE(id_ps_obs) ! PS
    CALL itpl_2d(v2d(:,:,iv2d_ps),ri,rj,yobs)
  CASE(id_rh_obs) ! RH
    DO k=ks,ke
      DO j=js,je
        IF(ie <= nlon ) THEN
          CALL calc_rh(v3d(is,j,k,iv3d_t),v3d(is,j,k,iv3d_q),&
            & p_full(is,j,k),rh(is,j,k))
          CALL calc_rh(v3d(ie,j,k,iv3d_t),v3d(ie,j,k,iv3d_q),&
            & p_full(ie,j,k),rh(ie,j,k))
        ELSE
          CALL calc_rh(v3d(is,j,k,iv3d_t),v3d(is,j,k,iv3d_q),&
            & p_full(is,j,k),rh(is,j,k))
          CALL calc_rh(v3d( 1,j,k,iv3d_t),v3d( 1,j,k,iv3d_q),&
            & p_full( 1,j,k),rh( 1,j,k))
        END IF
      END DO
    END DO
    CALL itpl_3d(rh,ri,rj,rk,yobs)
  END SELECT

  RETURN
END SUBROUTINE Trans_XtoY
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
  REAL(r_size),INTENT(IN) :: t  ! temperature [K] at target level
  REAL(r_size),INTENT(IN) :: q  ! humidity [kg/kg] at target level
  REAL(r_size),PARAMETER :: gamma=5.0d-3 ! lapse rate [K/m]
  REAL(r_size) :: tv

  tv = t * (1.0d0 + 0.608d0 * q)
  IF(dz /= 0) THEN
!    p = p * ((-gamma*dz+tv)/tv)**(gg/(gamma*rd)) !tv is at original level
    p = p * (tv/(tv+gamma*dz))**(gg/(gamma*rd)) !tv is at target level
  END IF

  RETURN
END SUBROUTINE prsadj
!-----------------------------------------------------------------------
! Coordinate conversion
!-----------------------------------------------------------------------
SUBROUTINE phys2ijk(p_full,elem,rlon,rlat,rlev,ri,rj,rk)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: p_full(nlon,nlat,nlev)
  REAL(r_size),INTENT(IN) :: elem
  REAL(r_size),INTENT(IN) :: rlon
  REAL(r_size),INTENT(IN) :: rlat
  REAL(r_size),INTENT(IN) :: rlev ! pressure levels
  REAL(r_size),INTENT(OUT) :: ri
  REAL(r_size),INTENT(OUT) :: rj
  REAL(r_size),INTENT(OUT) :: rk
  REAL(r_size) :: aj,ak
  REAL(r_size) :: lnps(nlon,nlat)
  REAL(r_size) :: plev(nlev)
  INTEGER :: i,j,k
!
! rlon -> ri
!
  IF(rlon == 0.0 .OR. rlon == 360.0) THEN
    ri = REAL(nlon+1,r_size)
  ELSE
    ri = rlon / 360.0d0 * REAL(nlon,r_size) + 1.0d0
  END IF
  IF(CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri)) RETURN
!
! rlat -> rj
!
  DO j=1,nlat
    IF(rlat < lat(j)) EXIT
  END DO
  IF(j == 1) THEN
    rj = (rlat + 90.0d0) / (lat(1) + 90.0d0)
  ELSE IF(j == nlat+1) THEN
    aj = (rlat - lat(nlat)) / (90.0d0 - lat(nlat))
    rj = REAL(nlat,r_size) + aj
  ELSE
    aj = (rlat - lat(j-1)) / (lat(j) - lat(j-1))
    rj = REAL(j-1,r_size) + aj
  END IF
  IF(CEILING(rj) < 2 .OR. nlat < CEILING(rj)) RETURN
!
! rlev -> rk
!
  IF(NINT(elem) == id_ps_obs) THEN ! surface pressure observation
    rk = 0.0d0
  ELSE
    !
    ! horizontal interpolation
    !
    i = CEILING(ri)
    j = CEILING(rj)
    DO k=1,nlev
      IF(i <= nlon) THEN
        lnps(i-1:i,j-1:j) = LOG(p_full(i-1:i,j-1:j,k))
      ELSE
        lnps(i-1,j-1:j) = LOG(p_full(i-1,j-1:j,k))
        lnps(1,j-1:j) = LOG(p_full(1,j-1:j,k))
      END IF
      CALL itpl_2d(lnps,ri,rj,plev(k))
    END DO
    !
    ! Log pressure
    !
    rk = LOG(rlev)
    !
    ! find rk
    !
    DO k=2,nlev-1
      IF(plev(k) < rk) EXIT ! assuming descending order of plev
    END DO
    ak = (rk - plev(k-1)) / (plev(k) - plev(k-1))
    rk = REAL(k-1,r_size) + ak
  END IF

  RETURN
END SUBROUTINE phys2ijk
!-----------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------
SUBROUTINE itpl_2d(var,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlon,nlat)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  IF(i <= nlon) THEN
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(i  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(i  ,j  ) *    ai  *    aj
  ELSE
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(1  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(1  ,j  ) *    ai  *    aj
  END IF

  RETURN
END SUBROUTINE itpl_2d

SUBROUTINE itpl_3d(var,ri,rj,rk,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlon,nlat,nlev)
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

  IF(i <= nlon) THEN
    var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(i  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
       & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
       & + var(i  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
       & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
       & + var(i  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
       & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
       & + var(i  ,j  ,k  ) *    ai  *    aj  *    ak
  ELSE
    var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(1  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
       & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
       & + var(1  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
       & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
       & + var(1  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
       & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
       & + var(1  ,j  ,k  ) *    ai  *    aj  *    ak
  END IF

  RETURN
END SUBROUTINE itpl_3d
!-----------------------------------------------------------------------
! Monitor departure
!-----------------------------------------------------------------------
SUBROUTINE monit_dep(nn,elm,dep,qc)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elm(nn)
  REAL(r_size),INTENT(IN) :: dep(nn)
  INTEGER,INTENT(IN) :: qc(nn)
  REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_q,rmse_ps,rmse_rh
  REAL(r_size) :: bias_u,bias_v,bias_t,bias_q,bias_ps,bias_rh
  INTEGER :: n,iu,iv,it,iq,ips,irh

  rmse_u = 0.0d0
  rmse_v = 0.0d0
  rmse_t = 0.0d0
  rmse_q = 0.0d0
  rmse_ps = 0.0d0
  rmse_rh = 0.0d0
  bias_u = 0.0d0
  bias_v = 0.0d0
  bias_t = 0.0d0
  bias_q = 0.0d0
  bias_ps = 0.0d0
  bias_rh = 0.0d0
  iu = 0
  iv = 0
  it = 0
  iq = 0
  ips = 0
  irh = 0
  DO n=1,nn
    IF(qc(n) /= 1) CYCLE
    SELECT CASE(NINT(elm(n)))
    CASE(id_u_obs)
      rmse_u = rmse_u + dep(n)**2
      bias_u = bias_u + dep(n)
      iu = iu + 1
    CASE(id_v_obs)
      rmse_v = rmse_v + dep(n)**2
      bias_v = bias_v + dep(n)
      iv = iv + 1
    CASE(id_t_obs)
      rmse_t = rmse_t + dep(n)**2
      bias_t = bias_t + dep(n)
      it = it + 1
    CASE(id_q_obs)
      rmse_q = rmse_q + dep(n)**2
      bias_q = bias_q + dep(n)
      iq = iq + 1
    CASE(id_ps_obs)
      rmse_ps = rmse_ps + dep(n)**2
      bias_ps = bias_ps + dep(n)
      ips = ips + 1
    CASE(id_rh_obs)
      rmse_rh = rmse_rh + dep(n)**2
      bias_rh = bias_rh + dep(n)
      irh = irh + 1
    END SELECT
  END DO
  IF(iu == 0) THEN
    rmse_u = undef
    bias_u = undef
  ELSE
    rmse_u = SQRT(rmse_u / REAL(iu,r_size))
    bias_u = bias_u / REAL(iu,r_size)
  END IF
  IF(iv == 0) THEN
    rmse_v = undef
    bias_v = undef
  ELSE
    rmse_v = SQRT(rmse_v / REAL(iv,r_size))
    bias_v = bias_v / REAL(iv,r_size)
  END IF
  IF(it == 0) THEN
    rmse_t = undef
    bias_t = undef
  ELSE
    rmse_t = SQRT(rmse_t / REAL(it,r_size))
    bias_t = bias_t / REAL(it,r_size)
  END IF
  IF(iq == 0) THEN
    rmse_q = undef
    bias_q = undef
  ELSE
    rmse_q = SQRT(rmse_q / REAL(iq,r_size))
    bias_q = bias_q / REAL(iq,r_size)
  END IF
  IF(ips == 0) THEN
    rmse_ps = undef
    bias_ps = undef
  ELSE
    rmse_ps = SQRT(rmse_ps / REAL(ips,r_size))
    bias_ps = bias_ps / REAL(ips,r_size)
  END IF
  IF(irh == 0) THEN
    rmse_rh = undef
    bias_rh = undef
  ELSE
    rmse_rh = SQRT(rmse_rh / REAL(irh,r_size))
    bias_rh = bias_rh / REAL(irh,r_size)
  END IF

  WRITE(6,'(A)') '== OBSERVATIONAL DEPARTURE ============================================='
  WRITE(6,'(6A12)') 'U','V','T','Q','PS','RH'
  WRITE(6,'(6ES12.3)') bias_u,bias_v,bias_t,bias_q,bias_ps,bias_rh
  WRITE(6,'(6ES12.3)') rmse_u,rmse_v,rmse_t,rmse_q,rmse_ps,rmse_rh
  WRITE(6,'(A)') '== NUMBER OF OBSERVATIONS TO BE ASSIMILATED ============================'
  WRITE(6,'(6A12)') 'U','V','T','Q','PS','RH'
  WRITE(6,'(6I12)') iu,iv,it,iq,ips,irh
  WRITE(6,'(A)') '========================================================================'

  RETURN
END SUBROUTINE monit_dep
!-----------------------------------------------------------------------
! Monitor departure from gues/anal mean
!-----------------------------------------------------------------------
SUBROUTINE monit_mean(file)
  CHARACTER(4),INTENT(IN) :: file
  REAL(r_size) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size) :: v2d(nlon,nlat,nv2d)
  REAL(r_size) :: p_full(nlon,nlat,nlev)
  REAL(r_size) :: elem
  REAL(r_size) :: bias_u,bias_v,bias_t,bias_ps,bias_q,bias_rh
  REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_ps,rmse_q,rmse_rh
  REAL(r_size) :: hdxf,dep,ri,rj,rk
  INTEGER :: n,iu,iv,it,iq,ips,irh
  CHARACTER(11) :: filename='filexxx.grd'

  rmse_u  = 0.0d0
  rmse_v  = 0.0d0
  rmse_t  = 0.0d0
  rmse_q  = 0.0d0
  rmse_ps = 0.0d0
  rmse_rh = 0.0d0
  bias_u = 0.0d0
  bias_v = 0.0d0
  bias_t = 0.0d0
  bias_q = 0.0d0
  bias_ps = 0.0d0
  bias_rh = 0.0d0
  iu  = 0
  iv  = 0
  it  = 0
  iq  = 0
  ips = 0
  irh = 0

  WRITE(filename(1:7),'(A4,A3)') file,'_me'
  CALL read_grd(filename,v3d,v2d)
  CALL calc_pfull(nlon,nlat,v2d(:,:,iv2d_ps),p_full)

  DO n=1,nobs
    CALL phys2ijk(p_full,obselm(n),obslon(n),obslat(n),obslev(n),ri,rj,rk)
    IF(CEILING(rk) > nlev) CYCLE
    IF(CEILING(rk) < 2 .AND. NINT(obselm(n)) /= id_ps_obs) THEN
      IF(NINT(obselm(n)) == id_u_obs .OR. NINT(obselm(n)) == id_v_obs) THEN
        rk = 1.00001d0
      ELSE
        CYCLE
      END IF
    END IF
    CALL Trans_XtoY(obselm(n),ri,rj,rk,v3d,v2d,p_full,hdxf)
    dep = obsdat(n) - hdxf
    SELECT CASE(NINT(obselm(n)))
    CASE(id_u_obs)
      rmse_u = rmse_u + dep**2
      bias_u = bias_u + dep
      iu = iu + 1
    CASE(id_v_obs)
      rmse_v = rmse_v + dep**2
      bias_v = bias_v + dep
      iv = iv + 1
    CASE(id_t_obs)
      rmse_t = rmse_t + dep**2
      bias_t = bias_t + dep
      it = it + 1
    CASE(id_q_obs)
      rmse_q = rmse_q + dep**2
      bias_q = bias_q + dep
      iq = iq + 1
    CASE(id_ps_obs)
      rmse_ps = rmse_ps + dep**2
      bias_ps = bias_ps + dep
      ips = ips + 1
    CASE(id_rh_obs)
      rmse_rh = rmse_rh + dep**2
      bias_rh = bias_rh + dep
      irh = irh + 1
    END SELECT
  END DO

  IF(iu == 0) THEN
    rmse_u = undef
    bias_u = undef
  ELSE
    rmse_u = SQRT(rmse_u / REAL(iu,r_size))
    bias_u = bias_u / REAL(iu,r_size)
  END IF
  IF(iv == 0) THEN
    rmse_v = undef
    bias_v = undef
  ELSE
    rmse_v = SQRT(rmse_v / REAL(iv,r_size))
    bias_v = bias_v / REAL(iv,r_size)
  END IF
  IF(it == 0) THEN
    rmse_t = undef
    bias_t = undef
  ELSE
    rmse_t = SQRT(rmse_t / REAL(it,r_size))
    bias_t = bias_t / REAL(it,r_size)
  END IF
  IF(iq == 0) THEN
    rmse_q = undef
    bias_q = undef
  ELSE
    rmse_q = SQRT(rmse_q / REAL(iq,r_size))
    bias_q = bias_q / REAL(iq,r_size)
  END IF
  IF(ips == 0) THEN
    rmse_ps = undef
    bias_ps = undef
  ELSE
    rmse_ps = SQRT(rmse_ps / REAL(ips,r_size))
    bias_ps = bias_ps / REAL(ips,r_size)
  END IF
  IF(irh == 0) THEN
    rmse_rh = undef
    bias_rh = undef
  ELSE
    rmse_rh = SQRT(rmse_rh / REAL(irh,r_size))
    bias_rh = bias_rh / REAL(irh,r_size)
  END IF

  WRITE(6,'(3A)') '== PARTIAL OBSERVATIONAL DEPARTURE (',file,') =============================='
  WRITE(6,'(6A12)') 'U','V','T','Q','PS','RH'
  WRITE(6,'(6ES12.3)') bias_u,bias_v,bias_t,bias_q,bias_ps,bias_rh
  WRITE(6,'(6ES12.3)') rmse_u,rmse_v,rmse_t,rmse_q,rmse_ps,rmse_rh
  WRITE(6,'(A)') '== NUMBER OF OBSERVATIONS =============================================='
  WRITE(6,'(6A12)') 'U','V','T','Q','PS','RH'
  WRITE(6,'(6I12)') iu,iv,it,iq,ips,irh
  WRITE(6,'(A)') '========================================================================'

  RETURN

END SUBROUTINE monit_mean
!-----------------------------------------------------------------------
! Basic modules for observation input
!-----------------------------------------------------------------------
SUBROUTINE get_nobs_mpi(cfile,nn)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(OUT) :: nn
  REAL(r_sngl) :: wk(6)
  INTEGER :: ios
  INTEGER :: iu,iv,it,iq,irh,ips
  INTEGER :: iunit
  LOGICAL :: ex

  nn = 0
  iu = 0
  iv = 0
  it = 0
  iq = 0
  irh = 0
  ips = 0
  iunit=91
  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is accessing a file ',cfile
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
    DO
      READ(iunit,IOSTAT=ios) wk
      IF(ios /= 0) EXIT
      SELECT CASE(NINT(wk(1)))
      CASE(id_u_obs)
        iu = iu + 1
      CASE(id_v_obs)
        iv = iv + 1
      CASE(id_t_obs)
        it = it + 1
      CASE(id_q_obs)
        iq = iq + 1
      CASE(id_rh_obs)
        irh = irh + 1
      CASE(id_ps_obs)
        ips = ips + 1
      END SELECT
      nn = nn + 1
    END DO
    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
    WRITE(6,'(A12,I10)') '          U:',iu
    WRITE(6,'(A12,I10)') '          V:',iv
    WRITE(6,'(A12,I10)') '          T:',it
    WRITE(6,'(A12,I10)') '          Q:',iq
    WRITE(6,'(A12,I10)') '         RH:',irh
    WRITE(6,'(A12,I10)') '         Ps:',ips
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF

  RETURN
END SUBROUTINE get_nobs_mpi

SUBROUTINE read_obs_mpi(cfile,nn,elem,rlon,rlat,rlev,odat,oerr)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(OUT) :: elem(nn) ! element number
  REAL(r_size),INTENT(OUT) :: rlon(nn)
  REAL(r_size),INTENT(OUT) :: rlat(nn)
  REAL(r_size),INTENT(OUT) :: rlev(nn)
  REAL(r_size),INTENT(OUT) :: odat(nn)
  REAL(r_size),INTENT(OUT) :: oerr(nn)
  REAL(r_sngl) :: wk(6)
  INTEGER :: n,iunit

  iunit=91
  WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',cfile
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,nn
    READ(iunit) wk
    SELECT CASE(NINT(wk(1)))
    CASE(id_u_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_v_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_t_obs)
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
    END SELECT
    elem(n) = REAL(wk(1),r_size)
    rlon(n) = REAL(wk(2),r_size)
    rlat(n) = REAL(wk(3),r_size)
    rlev(n) = REAL(wk(4),r_size)
    odat(n) = REAL(wk(5),r_size)
    oerr(n) = REAL(wk(6),r_size)
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs_mpi

END MODULE common_obs_speedy
