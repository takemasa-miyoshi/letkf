PROGRAM ave_diff
  USE common
  USE common_roms

  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER,PARAMETER :: itstart=41
  INTEGER,PARAMETER :: itend=120
  CHARACTER(14) :: ncinfile1='input1/0000.nc'
  CHARACTER(14) :: ncinfile2='input2/0000.nc'
!                             12345678901234
  CHARACTER(6) :: ncoutfile='ave.nc'
  REAL(r_size) :: v3d1(nlon,nlat,nlev,nv3d)
  REAL(r_size) :: v2d1(nlon,nlat,nv2d)
  REAL(r_size) :: v3d2(nlon,nlat,nlev,nv3d)
  REAL(r_size) :: v2d2(nlon,nlat,nv2d)
  REAL(r_size) :: av3d(nlon,nlat,nlev,nv3d)
  REAL(r_size) :: av2d(nlon,nlat,nv2d)
  INTEGER :: n,it
  REAL(r_size) :: diag(2)

  av3d = 0.0d0
  av2d = 0.0d0
  n=0
  DO it=itstart,itend
    WRITE(ncinfile1(8:11),'(I4.4)') it
    WRITE(ncinfile2(8:11),'(I4.4)') it
    CALL read_grd(ncinfile1,v3d1,v2d1)
    CALL read_grd(ncinfile2,v3d2,v2d2)
    av3d = av3d + (v3d1-v3d2)**2
    av2d = av2d + (v2d1-v2d2)**2
    n = n+1
  END DO
  av3d = av3d / REAL(n,r_size)
  av2d = av2d / REAL(n,r_size)
  CALL write_grd(ncoutfile,av3d,av2d)

  STOP
END PROGRAM ave_diff
