PROGRAM obsmake
!=======================================================================
!
! [PURPOSE:] For generating observations
!
! [HISTORY:]
!   01/26/2009 Takemasa Miyoshi  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_obs
  USE common_roms

  IMPLICIT NONE

  LOGICAL,PARAMETER :: verbose=.FALSE.

  INTEGER,PARAMETER :: obkind=5
  INTEGER,PARAMETER :: kmin=10
  INTEGER,PARAMETER :: kmax=40
  INTEGER,PARAMETER :: kint=3
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: error(:)
  REAL(r_sngl) :: wk(6)
  REAL(r_size) :: oberr(obkind)
  INTEGER :: obelm(obkind)
  LOGICAL :: ex
  REAL(r_size) :: a
  INTEGER :: i,j,k,n,nn,ios
  INTEGER :: nobkind,nstation
  CHARACTER(100) :: cdummy

  CHARACTER(7) :: truef='true.nc'
  CHARACTER(7) :: obsf='obs.dat'
  !
  ! Initial settings
  !
  CALL set_common_roms
  !
  ! Read true (nature run)
  !
  ALLOCATE(v3d(nlon,nlat,nlev,nv3d))
  ALLOCATE(v2d(nlon,nlat,nv2d))
  CALL read_grd(truef,v3d,v2d)
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
  ! Locations of station (station.tbl)
  !
  OPEN(10,FILE='station.tbl')
  READ(10,'(A)') cdummy
  READ(10,'(A)') cdummy
  nstation = 0
  DO
    READ(10,'(2I3)',IOSTAT=ios) i,j
    IF(ios /= 0) EXIT
    nstation = nstation + 1
  END DO
  CLOSE(10)
  IF(verbose) PRINT '(A,I)','nstation = ',nstation
  !
  ! Count number of obs
  !
  nobs = 0
  DO n=1,nobkind
    IF(obelm(n) == id_z_obs) THEN
      nobs = nobs + nstation
    ELSE
      nobs = nobs + nstation*((kmax-kmin)/kint+1)
    END IF
  END DO
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
  OPEN(90,FILE=obsf,FORM='unformatted',ACCESS='sequential')
  OPEN(10,FILE='station.tbl')
  READ(10,'(A)') cdummy
  READ(10,'(A)') cdummy
  DO
    READ(10,'(2I4)',IOSTAT=ios) i,j
    IF(ios /= 0) EXIT
    DO n=1,nobkind
      IF(obelm(n) == id_z_obs) THEN
        wk(1)=REAL(obelm(n),r_sngl)
        wk(2)=REAL(i,r_sngl)
        wk(3)=REAL(j,r_sngl)
        wk(4)=REAL(nlev,r_sngl)
        nn = nn+1
        wk(5)=REAL(v2d(i,j,iv2d_z)+error(nn)*oberr(n),r_sngl)
        wk(6)=REAL(oberr(n),r_sngl)
        WRITE(90) wk
      ELSE
        wk(1)=REAL(obelm(n),r_sngl)
        wk(2)=REAL(i,r_sngl)
        wk(3)=REAL(j,r_sngl)
        wk(6)=REAL(oberr(n),r_sngl)
        DO k=kmin,kmax,kint
          wk(4)=REAL(k,r_sngl)
          nn = nn+1
          SELECT CASE(obelm(n))
          CASE(id_u_obs)
            wk(5)=REAL(v3d(i,j,k,iv3d_u)+error(nn)*oberr(n),r_sngl)
          CASE(id_v_obs)
            wk(5)=REAL(v3d(i,j,k,iv3d_v)+error(nn)*oberr(n),r_sngl)
          CASE(id_t_obs)
            wk(5)=REAL(v3d(i,j,k,iv3d_t)+error(nn)*oberr(n),r_sngl)
          CASE(id_s_obs)
            wk(5)=REAL(v3d(i,j,k,iv3d_s)+error(nn)*oberr(n),r_sngl)
          END SELECT
          WRITE(90) wk
        END DO
      END IF
    END DO
  END DO
  CLOSE(90)
  !
  ! Count number of obs
  !
  IF(verbose) THEN
    OPEN(10,FILE=obsf,FORM='unformatted',ACCESS='sequential')
    CALL get_nobs(10)
    CLOSE(10)
  END IF

  DEALLOCATE(v3d,v2d,error)

  STOP
END PROGRAM obsmake
