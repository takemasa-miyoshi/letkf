MODULE common_roms
!=======================================================================
!
! [PURPOSE:] Common Information for ROMS
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   02/02/2009 Takemasa Miyoshi  modified for ROMS
!
!=======================================================================
!$USE OMP_LIB
  USE common
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  INTEGER,PARAMETER :: nlon=114
  INTEGER,PARAMETER :: nlat=114
  INTEGER,PARAMETER :: nlev=40
  INTEGER,PARAMETER :: nv3d=4 ! u,v,t,s
  INTEGER,PARAMETER :: nv2d=4 ! z,ubar,vbar,hbl
  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_t=3
  INTEGER,PARAMETER :: iv3d_s=4
  INTEGER,PARAMETER :: iv2d_z=1
  INTEGER,PARAMETER :: iv2d_ubar=2
  INTEGER,PARAMETER :: iv2d_vbar=3
  INTEGER,PARAMETER :: iv2d_hbl=4
  INTEGER,PARAMETER :: nij0=nlon*nlat
  INTEGER,PARAMETER :: nlevall=nlev*nv3d+nv2d
  INTEGER,PARAMETER :: ngpv=nij0*nlevall
  REAL(r_size),SAVE :: lon(nlon,nlat)
  REAL(r_size),SAVE :: lat(nlon,nlat)
  REAL(r_size),SAVE :: fcori(nlon,nlat)
  REAL(r_size),SAVE :: phi0(nlon,nlat)
  CHARACTER(4),SAVE :: element(nv3d+nv2d)

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_roms
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER :: ncid,istat,varid

  WRITE(6,'(A)') 'Hello from set_common_roms'
  !
  ! Elements
  !
  element(iv3d_u) = 'U   '
  element(iv3d_v) = 'V   '
  element(iv3d_t) = 'T   '
  element(iv3d_s) = 'SALT'
  element(nv3d+iv2d_z) = 'ZETA'
  element(nv3d+iv2d_ubar) = 'UBAR'
  element(nv3d+iv2d_vbar) = 'VBAR'
  element(nv3d+iv2d_hbl) = 'HBL '
  !
  ! Lon, Lat, f, orography
  !
  WRITE(6,'(A)') '  >> accessing to file: grd.nc'
  istat = NF_OPEN('grd.nc',NF_NOWRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF OPEN ERROR'
    STOP
  END IF
  istat = NF_INQ_VARID(ncid,'lon_rho',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,lon)
  istat = NF_INQ_VARID(ncid,'lat_rho',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,lat)
  istat = NF_INQ_VARID(ncid,'f',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,fcori)
  istat = NF_INQ_VARID(ncid,'h',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,phi0)
  istat = NF_CLOSE(ncid)

  RETURN
END SUBROUTINE set_common_roms
!-----------------------------------------------------------------------
! File I/O (netCDF)
!-----------------------------------------------------------------------
!-- Read a grid file ---------------------------------------------------
SUBROUTINE read_grd(filename,v3d,v2d)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl) :: buf4(nlon,nlat,nlev)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid

  istat = NF_OPEN(filename,NF_NOWRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF OPEN ERROR'
    STOP
  END IF
  !!! z
  istat = NF_INQ_VARID(ncid,'zeta',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,buf4(:,:,1))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (zeta)'
    STOP
  END IF
  DO j=1,nlat
    DO i=1,nlon
      v2d(i,j,iv2d_z) = REAL(buf4(i,j,1),r_size)
    END DO
  END DO
  !!! ubar
  istat = NF_INQ_VARID(ncid,'ubar',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,buf4(1:nlon-1,:,1))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (ubar)'
    STOP
  END IF
  DO j=1,nlat
    DO i=1,nlon-1
      v2d(i,j,iv2d_ubar) = REAL(buf4(i,j,1),r_size)
    END DO
  END DO
  DO j=1,nlat
    v2d(nlon,j,iv2d_ubar) = 0.0d0
  END DO
  !!! vbar
  istat = NF_INQ_VARID(ncid,'vbar',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,buf4(:,1:nlat-1,1))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (vbar)'
    STOP
  END IF
  DO j=1,nlat-1
    DO i=1,nlon
      v2d(i,j,iv2d_vbar) = REAL(buf4(i,j,1),r_size)
    END DO
  END DO
  DO i=1,nlon
    v2d(i,nlat,iv2d_vbar) = 0.0d0
  END DO
  !!! u
  istat = NF_INQ_VARID(ncid,'u',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,buf4(1:nlon-1,:,:))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (u)'
    STOP
  END IF
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon-1
        v3d(i,j,k,iv3d_u) = REAL(buf4(i,j,k),r_size)
      END DO
    END DO
  END DO
  DO k=1,nlev
    DO j=1,nlat
      v3d(nlon,j,k,iv3d_u) = 0.0d0
    END DO
  END DO
  !!! v
  istat = NF_INQ_VARID(ncid,'v',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,buf4(:,1:nlat-1,:))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (v)'
    STOP
  END IF
  DO k=1,nlev
    DO j=1,nlat-1
      DO i=1,nlon
        v3d(i,j,k,iv3d_v) = REAL(buf4(i,j,k),r_size)
      END DO
    END DO
  END DO
  DO k=1,nlev
    DO i=1,nlon
      v3d(i,nlat,k,iv3d_v) = 0.0d0
    END DO
  END DO
  !!! t
  istat = NF_INQ_VARID(ncid,'temp',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,buf4)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (temp)'
    STOP
  END IF
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        v3d(i,j,k,iv3d_t) = REAL(buf4(i,j,k),r_size)
      END DO
    END DO
  END DO
  !!! s
  istat = NF_INQ_VARID(ncid,'salt',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,buf4)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (salt)'
    STOP
  END IF
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        v3d(i,j,k,iv3d_s) = REAL(buf4(i,j,k),r_size)
      END DO
    END DO
  END DO
  !!! hbl
  istat = NF_INQ_VARID(ncid,'hbl',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,buf4(:,:,1))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (hbl)'
    STOP
  END IF
  DO j=1,nlat
    DO i=1,nlon
      v2d(i,j,iv2d_hbl) = REAL(buf4(i,j,1),r_size)
    END DO
  END DO

  istat = NF_CLOSE(ncid)

  RETURN
END SUBROUTINE read_grd

SUBROUTINE read_grd4(filename,v3d,v2d)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER :: ncid,istat,varid

  istat = NF_OPEN(filename,NF_NOWRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF OPEN ERROR'
    STOP
  END IF
  !!! z
  istat = NF_INQ_VARID(ncid,'zeta',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v2d(:,:,iv2d_z))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (zeta)'
    STOP
  END IF
  !!! ubar
  istat = NF_INQ_VARID(ncid,'ubar',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v2d(1:nlon-1,:,iv2d_ubar))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (ubar)'
    STOP
  END IF
  v2d(nlon,:,iv2d_ubar) = 0.0
  !!! vbar
  istat = NF_INQ_VARID(ncid,'vbar',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v2d(:,1:nlat-1,iv2d_vbar))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (vbar)'
    STOP
  END IF
  v2d(:,nlat,iv2d_vbar) = 0.0
  !!! u
  istat = NF_INQ_VARID(ncid,'u',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v3d(1:nlon-1,:,:,iv3d_u))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (u)'
    STOP
  END IF
  v3d(nlon,:,:,iv3d_u) = 0.0
  !!! v
  istat = NF_INQ_VARID(ncid,'v',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v3d(:,1:nlat-1,:,iv3d_v))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (v)'
    STOP
  END IF
  v3d(:,nlat,:,iv3d_v) = 0.0
  !!! t
  istat = NF_INQ_VARID(ncid,'temp',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_t))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (temp)'
    STOP
  END IF
  !!! s
  istat = NF_INQ_VARID(ncid,'salt',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_s))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (salt)'
    STOP
  END IF
  !!! hbl
  istat = NF_INQ_VARID(ncid,'hbl',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v2d(:,:,iv2d_hbl))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (hbl)'
    STOP
  END IF

  istat = NF_CLOSE(ncid)

  RETURN
END SUBROUTINE read_grd4
!-- Write a grid file -------------------------------------------------
SUBROUTINE write_grd(filename,v3d,v2d)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl) :: buf4(nlon,nlat,nlev)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid

  istat = NF_OPEN(filename,NF_WRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF OPEN ERROR'
    STOP
  END IF
  !!! z
  istat = NF_INQ_VARID(ncid,'zeta',varid)
  DO j=1,nlat
    DO i=1,nlon
      buf4(i,j,1) = REAL(v2d(i,j,iv2d_z),r_sngl)
    END DO
  END DO
  istat = NF_PUT_VAR_REAL(ncid,varid,buf4(:,:,1))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (zeta)'
    STOP
  END IF
  !!! ubar
  istat = NF_INQ_VARID(ncid,'ubar',varid)
  DO j=1,nlat
    DO i=1,nlon-1
      buf4(i,j,1) = REAL(v2d(i,j,iv2d_ubar),r_sngl)
    END DO
  END DO
  istat = NF_PUT_VAR_REAL(ncid,varid,buf4(1:nlon-1,:,1))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (ubar)'
    STOP
  END IF
  !!! vbar
  istat = NF_INQ_VARID(ncid,'vbar',varid)
  DO j=1,nlat-1
    DO i=1,nlon
      buf4(i,j,1) = REAL(v2d(i,j,iv2d_vbar),r_sngl)
    END DO
  END DO
  istat = NF_PUT_VAR_REAL(ncid,varid,buf4(:,1:nlat-1,1))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (vbar)'
    STOP
  END IF
  !!! u
  istat = NF_INQ_VARID(ncid,'u',varid)
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon-1
        buf4(i,j,k) = REAL(v3d(i,j,k,iv3d_u),r_sngl)
      END DO
    END DO
  END DO
  istat = NF_PUT_VAR_REAL(ncid,varid,buf4(1:nlon-1,:,:))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (u)'
    STOP
  END IF
  !!! v
  istat = NF_INQ_VARID(ncid,'v',varid)
  DO k=1,nlev
    DO j=1,nlat-1
      DO i=1,nlon
        buf4(i,j,k) = REAL(v3d(i,j,k,iv3d_v),r_sngl)
      END DO
    END DO
  END DO
  istat = NF_PUT_VAR_REAL(ncid,varid,buf4(:,1:nlat-1,:))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (v)'
    STOP
  END IF
  !!! t
  istat = NF_INQ_VARID(ncid,'temp',varid)
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        buf4(i,j,k) = REAL(v3d(i,j,k,iv3d_t),r_sngl)
      END DO
    END DO
  END DO
  istat = NF_PUT_VAR_REAL(ncid,varid,buf4)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (temp)'
    STOP
  END IF
  !!! s
  istat = NF_INQ_VARID(ncid,'salt',varid)
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        buf4(i,j,k) = REAL(v3d(i,j,k,iv3d_s),r_sngl)
      END DO
    END DO
  END DO
  istat = NF_PUT_VAR_REAL(ncid,varid,buf4)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (salt)'
    STOP
  END IF
  !!! hbl
  istat = NF_INQ_VARID(ncid,'hbl',varid)
  DO j=1,nlat
    DO i=1,nlon
      buf4(i,j,1) = REAL(v2d(i,j,iv2d_hbl),r_sngl)
    END DO
  END DO
  istat = NF_PUT_VAR_REAL(ncid,varid,buf4(:,:,1))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (hbl)'
    STOP
  END IF

  istat = NF_CLOSE(ncid)

  RETURN
END SUBROUTINE write_grd

SUBROUTINE write_grd4(filename,v3d,v2d)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: ncid,istat,varid

  istat = NF_OPEN(filename,NF_WRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF OPEN ERROR'
    STOP
  END IF
  !!! z
  istat = NF_INQ_VARID(ncid,'zeta',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v2d(:,:,iv2d_z))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (zeta)'
    STOP
  END IF
  !!! ubar
  istat = NF_INQ_VARID(ncid,'ubar',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v2d(1:nlon-1,:,iv2d_ubar))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (ubar)'
    STOP
  END IF
  !!! vbar
  istat = NF_INQ_VARID(ncid,'vbar',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v2d(:,1:nlat-1,iv2d_vbar))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (vbar)'
    STOP
  END IF
  !!! u
  istat = NF_INQ_VARID(ncid,'u',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v3d(1:nlon-1,:,:,iv3d_u))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (u)'
    STOP
  END IF
  !!! v
  istat = NF_INQ_VARID(ncid,'v',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v3d(:,1:nlat-1,:,iv3d_v))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (v)'
    STOP
  END IF
  !!! t
  istat = NF_INQ_VARID(ncid,'temp',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_t))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (temp)'
    STOP
  END IF
  !!! s
  istat = NF_INQ_VARID(ncid,'salt',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_s))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (salt)'
    STOP
  END IF
  !!! hbl
  istat = NF_INQ_VARID(ncid,'hbl',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v2d(:,:,iv2d_hbl))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (hbl)'
    STOP
  END IF

  istat = NF_CLOSE(ncid)

  RETURN
END SUBROUTINE write_grd4
!-----------------------------------------------------------------------
! Depth
!-----------------------------------------------------------------------
SUBROUTINE calc_depth(zeta,bottom,depth)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: zeta
  REAL(r_size),INTENT(IN) :: bottom
  REAL(r_size),INTENT(OUT) :: depth(nlev)
  REAL(r_size) :: cs(nlev)
  REAL(r_size) :: surfd
  INTEGER :: k

  cs=(/ -0.921962804129276, -0.783682417800398, -0.666141723246888, -0.566230075492245, -0.481303371858989, -0.409114077354825, -0.347751744988103, -0.295592456917874, -0.251255848428496, -0.213568577393686, -0.181533272470966, -0.154302138258862, -0.131154518892926, -0.111477826314426, -0.0947513284896102, -0.0805323685443561, -0.0684446501117523, -0.0581682788712608, -0.0494312967346868, -0.0420024846361897, -0.0356852434564318, -0.0303123911431393, -0.0257417383369477, -0.0218523254140864, -0.0185412213614628, -0.0157207997682656, -0.0133164198448193, -0.0112644510982221, -0.00951058938897752, -0.00800841980405428, -0.0067181883136592, -0.00560574970434628, -0.00464166394613383, -0.00380041707951966, -0.00305974600362987, -0.00240004929641906, -0.00180386847462664, -0.00125542596534104, -0.000740207561862293, -0.000244578313804319 /)

  DO k=1,nlev
    depth(k) = zeta+(zeta+bottom)* &
    & (10.0d0*(REAL(k,r_size)-REAL(nlev,r_size)-0.5d0)/ &
    & REAL(nlev,r_size)+cs(k)*bottom)/(bottom+10.0d0)
  END DO
  surfd = depth(nlev)
  DO k=1,nlev
    depth(k) = depth(k) - surfd
  END DO

  RETURN
END SUBROUTINE calc_depth
!-----------------------------------------------------------------------
! Monitor
!-----------------------------------------------------------------------
SUBROUTINE monit_grd(v3d,v2d)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: k,n

  DO k=1,nlev
    WRITE(6,'(I2,A)') k,'th level'
    DO n=1,nv3d
      WRITE(6,'(A,2ES10.2)') element(n),MAXVAL(v3d(:,:,k,n)),MINVAL(v3d(:,:,k,n))
    END DO
  END DO

  DO n=1,nv2d
    WRITE(6,'(A,2ES10.2)') element(nv3d+n),MAXVAL(v2d(:,:,n)),MINVAL(v2d(:,:,n))
  END DO

  RETURN
END SUBROUTINE monit_grd
!-----------------------------------------------------------------------
! Ensemble manipulations
!-----------------------------------------------------------------------
SUBROUTINE ensmean_grd(member,nij,v3d,v2d,v3dm,v2dm)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: member
  INTEGER,INTENT(IN) :: nij
  REAL(r_size),INTENT(IN) :: v3d(nij,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij,member,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dm(nij,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dm(nij,nv2d)
  INTEGER :: i,k,m,n

  DO n=1,nv3d
!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO k=1,nlev
      DO i=1,nij
        v3dm(i,k,n) = v3d(i,k,1,n)
        DO m=2,member
          v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,m,n)
        END DO
        v3dm(i,k,n) = v3dm(i,k,n) / REAL(member,r_size)
      END DO
    END DO
!$OMP END PARALLEL DO
  END DO

  DO n=1,nv2d
    DO i=1,nij
      v2dm(i,n) = v2d(i,1,n)
      DO m=2,member
        v2dm(i,n) = v2dm(i,n) + v2d(i,m,n)
      END DO
      v2dm(i,n) = v2dm(i,n) / REAL(member,r_size)
    END DO
  END DO

  RETURN
END SUBROUTINE ensmean_grd

END MODULE common_roms
