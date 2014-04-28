PROGRAM choceantime
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER :: ncid,istat,varid,k
  CHARACTER(5) :: ncinfile='in.nc'
  CHARACTER(6) :: ncoutfile='out.nc'
  REAL :: ocean_time

  istat = NF_OPEN(ncinfile,NF_NOWRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    PRINT '(A)','netCDF OPEN ERROR'
    STOP
  END IF
  istat = NF_INQ_VARID(ncid,'ocean_time',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,ocean_time)
  istat = NF_CLOSE(ncid)

  istat = NF_OPEN(ncoutfile,NF_WRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    PRINT '(A)','netCDF OPEN ERROR'
    STOP
  END IF
  istat = NF_INQ_VARID(ncid,'ocean_time',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,ocean_time)
  istat = NF_CLOSE(ncid)

  STOP
END PROGRAM choceantime
