PROGRAM station
  USE common
  USE common_obs
  USE common_roms

  IMPLICIT NONE

  INTEGER,PARAMETER :: obkind=5
  INTEGER,PARAMETER :: imin=3
  INTEGER,PARAMETER :: imax=nlon-3-2
  INTEGER,PARAMETER :: iint=5
  INTEGER,PARAMETER :: jmin=3
  INTEGER,PARAMETER :: jmax=nlat-2
  INTEGER,PARAMETER :: jint=5
  INTEGER,PARAMETER :: kmin=10
  INTEGER,PARAMETER :: kmax=40
  INTEGER,PARAMETER :: kint=3
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
  STOP
END PROGRAM station
