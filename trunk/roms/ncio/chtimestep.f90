PROGRAM chtimestep
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER :: ncid,istat,varid,k
  CHARACTER(6) :: ncfile='grd.nc'
  INTEGER :: time_step(4)

  time_step(1) = 19
  time_step(2) = 1
  time_step(3) = 0
  time_step(4) = 0

  istat = NF_OPEN(ncfile,NF_WRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    PRINT '(A)','netCDF OPEN ERROR'
    STOP
  END IF
  istat = NF_INQ_VARID(ncid,'time_step',varid)
  istat = NF_PUT_VAR_INT(ncid,varid,time_step)
  istat = NF_CLOSE(ncid)

  STOP
END PROGRAM chtimestep
