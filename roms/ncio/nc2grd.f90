PROGRAM nc2grd
  USE common
  USE common_roms

  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER :: ncid,istat,varid,k
  CHARACTER(8) :: ncinfile='input.nc'
  CHARACTER(9) :: ncoutfile='output.nc'
  REAL(r_sngl) :: v2d(nlon,nlat)
  REAL(r_sngl) :: v3d(nlon,nlat,nlev)

  istat = NF_OPEN(ncinfile,NF_NOWRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    PRINT '(A)','netCDF OPEN ERROR'
    STOP
  END IF
  istat = NF_INQ_VARID(ncid,'salt',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v3d)
  istat = NF_CLOSE(ncid)
  DO k=1,nlev
    print *,SUM(v3d(10,:,k))/REAL(nlat),SUM(v3d(100,:,k))/REAL(nlat)
  END DO

  istat = NF_OPEN(ncoutfile,NF_WRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    PRINT '(A)','netCDF OPEN ERROR'
    STOP
  END IF
  istat = NF_INQ_VARID(ncid,'salt',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v3d)
  istat = NF_CLOSE(ncid)

  STOP
END PROGRAM nc2grd
