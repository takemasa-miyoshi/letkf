PROGRAM station
  USE common
  USE common_roms

  IMPLICIT NONE

  INTEGER,PARAMETER :: imin=3
  INTEGER,PARAMETER :: imax=nlon-3-2
  INTEGER,PARAMETER :: iint=5
!  INTEGER,PARAMETER :: jmin=3
!  INTEGER,PARAMETER :: jmax=nlat-2
!  INTEGER,PARAMETER :: jint=5
  INTEGER,PARAMETER :: jmin=7
  INTEGER,PARAMETER :: jmax=nlat-2
  INTEGER,PARAMETER :: jint=20
  INTEGER,PARAMETER :: id_hrf   =1
  INTEGER,PARAMETER :: id_prof  =2
  INTEGER,PARAMETER :: id_sst   =3
  INTEGER :: i,j
  !
  ! stdout
  !
  PRINT '(A)','KIND   I   J'
  PRINT '(A)','------------'

  DO j=jmin,jmax,jint
    DO i=imin,imax,iint
      PRINT '(3I4)',id_hrf,i,j
      PRINT '(3I4)',id_prof,i,j
      PRINT '(3I4)',id_sst,i,j
    END DO
  END DO

  STOP
END PROGRAM station
