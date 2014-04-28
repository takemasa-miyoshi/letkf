PROGRAM breeding
  USE common
  USE common_roms

  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
!  REAL(r_size) :: rescale_norm = 5.0D-1
!  REAL(r_size) :: rescale_norm = 1.0D-2
!  REAL(r_size) :: rescale_norm = 1.0D-3
  REAL(r_size) :: rescale_norm = 5.0D-3
  INTEGER :: ncid,istat,varid,k
  CHARACTER(9) :: naturefile='nature.nc'
  CHARACTER(6) :: ptbfile='ptb.nc'
  CHARACTER(5) :: bvfile='bv.nc'
  REAL(r_size) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size) :: v2d(nlon,nlat,nv2d)
  REAL(r_size) :: v3dp(nlon,nlat,nlev,nv3d)
  REAL(r_size) :: v2dp(nlon,nlat,nv2d)
  INTEGER :: i,j,it
  REAL(r_size) :: anorm,rescale

  CALL read_grd(naturefile,v3d,v2d)
  CALL read_grd(ptbfile,v3dp,v2dp)

  v3dp = v3dp - v3d
  v2dp = v2dp - v2d

  CALL write_grd(bvfile,v3dp,v2dp)

  CALL norm(v3dp,v2dp,anorm)

  rescale = MIN(rescale_norm/anorm,1.0d0)
  PRINT *,anorm,rescale

  v3dp = v3d + v3dp*rescale
  v2dp = v2d + v2dp*rescale

  CALL write_grd(ptbfile,v3dp,v2dp)

  STOP
CONTAINS
SUBROUTINE norm(v3d,v2d,anorm)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: anorm

!  anorm = SUM(v3d(:,:,nlev,iv3d_u)**2)/REAL((nlon-1)*nlat,r_size) &
!      & + SUM(v3d(:,:,nlev,iv3d_v)**2)/REAL(nlon*(nlat-1),r_size)
!  anorm = SUM(v2d(:,:,iv2d_ubar)**2)/REAL((nlon-1)*nlat,r_size) &
!      & + SUM(v2d(:,:,iv2d_vbar)**2)/REAL(nlon*(nlat-1),r_size)
!  anorm = SUM(v2d(:,:,iv2d_z)**2)/REAL(nlon*nlat,r_size)
  anorm = SUM(v3d(:,:,5,iv3d_t)**2)/REAL(nlon*nlat,r_size)
  anorm = SQRT(anorm)

  RETURN
END SUBROUTINE norm
END PROGRAM breeding
