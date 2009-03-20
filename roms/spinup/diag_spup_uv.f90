PROGRAM diag_spup_uv
  USE common
  USE common_roms

  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER :: ncid,istat,varid
  CHARACTER(8) :: ncinfile='input.nc'
  REAL(r_size) :: ubar(nlon,nlat)
  REAL(r_size) :: vbar(nlon,nlat)
  REAL(r_size) :: u(nlon,nlat,nlev)
  REAL(r_size) :: v(nlon,nlat,nlev)
  REAL(r_sngl) :: buf4(nlon,nlat,nlev)
  INTEGER :: varst2d(3),varct2d(3)
  INTEGER :: varst3d(4),varct3d(4)
  INTEGER :: i,j,k,it
  REAL(r_size) :: diag(2)

  istat = NF_OPEN(ncinfile,NF_NOWRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    PRINT '(A)','netCDF OPEN ERROR'
    STOP
  END IF

  DO it=1,30
    varst2d = (/1,1,it/)
    varst3d = (/1,1,1,it/)
    !!! ubar
    varct2d = (/nlon-1,nlat,1/)
    istat = NF_INQ_VARID(ncid,'ubar',varid)
    istat = NF_GET_VARA_REAL(ncid,varid,varst2d,varct2d,buf4(1:nlon-1,:,1))
    IF(istat /= NF_NOERR) THEN
!      PRINT '(A)','netCDF READ ERROR (ubar)'
      CYCLE
    END IF
    DO j=1,nlat
      DO i=1,nlon-1
        ubar(i,j) = REAL(buf4(i,j,1),r_size)
      END DO
    END DO
    ubar(nlon,:) = 0.0d0
    !!! vbar
    varct2d = (/nlon,nlat-1,1/)
    istat = NF_INQ_VARID(ncid,'vbar',varid)
    istat = NF_GET_VARA_REAL(ncid,varid,varst2d,varct2d,buf4(:,1:nlat-1,1))
    IF(istat /= NF_NOERR) THEN
!      PRINT '(A)','netCDF READ ERROR (vbar)'
      CYCLE
    END IF
    DO j=1,nlat-1
      DO i=1,nlon
        vbar(i,j) = REAL(buf4(i,j,1),r_size)
      END DO
    END DO
    vbar(:,nlat) = 0.0d0
    !!! u
    varct3d = (/nlon-1,nlat,nlev,1/)
    istat = NF_INQ_VARID(ncid,'u',varid)
    istat = NF_GET_VARA_REAL(ncid,varid,varst3d,varct3d,buf4(1:nlon-1,:,:))
    IF(istat /= NF_NOERR) THEN
!      PRINT '(A)','netCDF READ ERROR (u)'
      CYCLE
    END IF
    u(1:nlon-1,:,:) = REAL(buf4(1:nlon-1,:,:),r_size)
    u(nlon,:,:) = 0.0d0
    !!! v
    varct3d = (/nlon,nlat-1,nlev,1/)
    istat = NF_INQ_VARID(ncid,'v',varid)
    istat = NF_GET_VARA_REAL(ncid,varid,varst3d,varct3d,buf4(:,1:nlat-1,:))
    IF(istat /= NF_NOERR) THEN
!      PRINT '(A)','netCDF READ ERROR (v)'
      CYCLE
    END IF
    v(:,1:nlat-1,:) = REAL(buf4(:,1:nlat-1,:),r_size)
    v(:,nlat,:) = 0.0d0

    diag(1) = SUM(u(:,:,40)**2)/REAL((nlon-1)*nlat,r_size) &
       & + SUM(v(:,:,40)**2)/REAL(nlon*(nlat-1),r_size)
    diag(1) = SQRT(diag(1))

    diag(2) = SUM(ubar(:,:)**2)/REAL((nlon-1)*nlat,r_size) &
       & + SUM(vbar(:,:)**2)/REAL(nlon*(nlat-1),r_size)
    diag(2) = SQRT(diag(2))


    PRINT *,diag(1),diag(2)

  END DO
  istat = NF_CLOSE(ncid)

  STOP
END PROGRAM diag_spup_uv
