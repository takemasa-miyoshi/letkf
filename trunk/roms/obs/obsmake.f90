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
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: error(:)
  REAL(r_sngl) :: wk(6)
  REAL(r_size) :: err
  INTEGER :: i,j,k,elm,nn,ios
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
  ! Read station (station.tbl)
  !
  OPEN(10,FILE='station.tbl')
  READ(10,'(A)') cdummy
  READ(10,'(A)') cdummy
  nobs = 0
  DO
    READ(10,'(2I4,I3,I6,ES9.2)',IOSTAT=ios) i,j,k,elm,err
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
  OPEN(90,FILE=obsf,FORM='unformatted',ACCESS='sequential')
  OPEN(10,FILE='station.tbl')
  READ(10,'(A)') cdummy
  READ(10,'(A)') cdummy
  DO
    READ(10,'(2I4,I3,I6,ES9.2)',IOSTAT=ios) i,j,k,elm,err
    IF(ios /= 0) EXIT
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
    WRITE(90) wk
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
