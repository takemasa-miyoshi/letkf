PROGRAM station
  IMPLICIT NONE
  INTEGER,PARAMETER :: nlon=96
  INTEGER,PARAMETER :: nlat=nlon/2
  INTEGER,PARAMETER :: jj=3
  INTEGER,PARAMETER :: interval=2
  INTEGER :: i,j

  PRINT '(A)','  I  J'
  PRINT '(A)','------'
  DO j=1+jj,nlat-jj,interval
    DO i=1,nlon,interval
      PRINT '(2I3)',i,j
    END DO
  END DO
  STOP
END PROGRAM station
