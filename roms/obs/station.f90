PROGRAM station
  USE common
  USE common_obs
  USE common_roms

  IMPLICIT NONE

  INTEGER,PARAMETER :: msw_obsnet=2 !1:regular, 2:custom
  INTEGER,PARAMETER :: obkind=5
  INTEGER :: i,j,k,n,ios
  INTEGER :: nobkind
  LOGICAL :: ex
  REAL(r_size) :: a
  REAL(r_size) :: oberr(obkind)
  INTEGER :: obelm(obkind)
  CHARACTER(100) :: cdummy
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
  ! stdout
  !
  PRINT '(A)','   I   J  K  ELEM     OBER'
  PRINT '(A)','--------------------------'
  SELECT CASE(msw_obsnet)
  CASE(1)
    CALL regular
  CASE(2)
    CALL custom
  END SELECT

  STOP
CONTAINS
!
! regular observing network
!
SUBROUTINE regular
  INTEGER,PARAMETER :: imin=3
  INTEGER,PARAMETER :: imax=nlon-3-2
  INTEGER,PARAMETER :: iint=5
  INTEGER,PARAMETER :: jmin=3
  INTEGER,PARAMETER :: jmax=nlat-2
  INTEGER,PARAMETER :: jint=5
  INTEGER,PARAMETER :: kmin=10
  INTEGER,PARAMETER :: kmax=40
  INTEGER,PARAMETER :: kint=3
  DO k=kmin,kmax,kint
    DO j=jmin,jmax,jint
      DO i=imin,imax,iint
        DO n=1,nobkind
          IF(obelm(n) == id_z_obs) THEN
            IF(k == nlev) PRINT '(2I4,I3,I6,ES9.2)',i,j,k,obelm(n),oberr(n)
          ELSE
            PRINT '(2I4,I3,I6,ES9.2)',i,j,k,obelm(n),oberr(n)
          END IF
        END DO
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE regular
!
! custom observing network
!
SUBROUTINE custom
  INTEGER,PARAMETER :: imin=1
  INTEGER,PARAMETER :: imax=nlon-3
  INTEGER,PARAMETER :: iint=2
  INTEGER,PARAMETER :: jmin=1
  INTEGER,PARAMETER :: jmax=nlat
  INTEGER,PARAMETER :: jint=2
  INTEGER,PARAMETER :: nsounding=5
  INTEGER,PARAMETER :: nglider=10
  INTEGER :: nn
  REAL(r_size),ALLOCATABLE :: rnd(:)
  !
  ! SST
  !
  k = nlev
  DO n=1,nobkind
    IF(obelm(n) == id_t_obs) EXIT
  END DO
  DO j=jmin,jmax,jint
    DO i=imin,imax,iint
      PRINT '(2I4,I3,I6,ES9.2)',i,j,k,obelm(n),oberr(n)
    END DO
  END DO
  !
  ! SSH
  !
  k = nlev
  DO n=1,nobkind
    IF(obelm(n) == id_z_obs) EXIT
  END DO
  DO j=jmin,jmax,jint
    DO i=imin,imax,iint
      PRINT '(2I4,I3,I6,ES9.2)',i,j,k,obelm(n),oberr(n)
    END DO
  END DO
  !
  ! surface flow (only near coast)
  !
  k = nlev
  DO n=1,nobkind
    IF(obelm(n) == id_u_obs) EXIT
  END DO
  DO j=jmin,jmax
    DO i=imax-30,imax
      PRINT '(2I4,I3,I6,ES9.2)',i,j,k,obelm(n),oberr(n)
    END DO
  END DO
  DO n=1,nobkind
    IF(obelm(n) == id_v_obs) EXIT
  END DO
  DO j=jmin,jmax
    DO i=imax-30,imax
      PRINT '(2I4,I3,I6,ES9.2)',i,j,k,obelm(n),oberr(n)
    END DO
  END DO
  !
  ! soundings (random locations, deep)
  !
  n = nsounding*2
  ALLOCATE(rnd(n))
  CALL com_rand(n,rnd)
  DO nn=1,nsounding
    i = CEILING(rnd(nn*2-1)*(imax-2)+1)
    j = CEILING(rnd(nn*2)*(jmax-2)+1)
    DO k=1,nlev
      DO n=1,nobkind
        IF(obelm(n) == id_z_obs) CYCLE
        IF(obelm(n) == id_t_obs .AND. k == nlev) CYCLE
        IF(obelm(n) == id_u_obs .AND. i >= imax-30 .AND. k == nlev) CYCLE
        IF(obelm(n) == id_v_obs .AND. i >= imax-30 .AND. k == nlev) CYCLE
        PRINT '(2I4,I3,I6,ES9.2)',i,j,k,obelm(n),oberr(n)
      END DO
    END DO
  END DO
  DEALLOCATE(rnd)
  !
  ! gliders (random line-like in x direction near coast, shallow)
  !
  ALLOCATE(rnd(nglider))
  CALL com_rand(nglider,rnd)
  DO nn=1,nglider
    j = CEILING(rnd(nn)*(jmax-2)+1)
    DO k=20,nlev,2
      DO i=imax-30,imax
        DO n=1,nobkind
          IF(obelm(n) == id_z_obs) CYCLE
          IF(obelm(n) == id_t_obs .AND. k == nlev) CYCLE
          IF(obelm(n) == id_u_obs .AND. k == nlev) CYCLE
          IF(obelm(n) == id_v_obs .AND. k == nlev) CYCLE
          PRINT '(2I4,I3,I6,ES9.2)',i,j,k,obelm(n),oberr(n)
        END DO
      END DO
    END DO
  END DO
  DEALLOCATE(rnd)

  RETURN
END SUBROUTINE custom
END PROGRAM station
