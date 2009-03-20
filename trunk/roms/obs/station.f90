PROGRAM station
  IMPLICIT NONE
  INTEGER,PARAMETER :: nlon=111
  INTEGER,PARAMETER :: nlat=114
  INTEGER,PARAMETER :: jj=2
  INTEGER,PARAMETER :: interval=5
  INTEGER :: i,j

  PRINT '(A)','   I   J'
  PRINT '(A)','--------'
  DO j=1+jj,nlat-jj,interval
    DO i=1+jj,nlon-jj,interval
      PRINT '(2I4)',i,j
    END DO
  END DO
  STOP
END PROGRAM station
