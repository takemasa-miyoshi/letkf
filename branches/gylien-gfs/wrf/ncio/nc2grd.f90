PROGRAM nc2grd
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER(4),PARAMETER :: imt_grd_cnst = 10
  INTEGER(4),PARAMETER :: imt_grd      = 11
  INTEGER(4) :: nlon,nlat,nlev
  INTEGER(4) :: ncid,varid,dimids(3),shape(3)
  INTEGER(4) :: reclength,irec
  INTEGER(4) :: i,j,k
  INTEGER(4),ALLOCATABLE :: start(:),count(:)
  REAL(4),PARAMETER :: rd = 287.05
  REAL(4),PARAMETER :: cp = 7.0 / 2.0 * rd
  REAL(4),PARAMETER :: p0 = 1.0e+5
  REAL(4),PARAMETER :: t0 = 300.0
  REAL(4),ALLOCATABLE :: rlon(:,:)
  REAL(4),ALLOCATABLE :: rlat(:,:)
  REAL(4),ALLOCATABLE :: rlon_u(:,:)
  REAL(4),ALLOCATABLE :: rlat_u(:,:)
  REAL(4),ALLOCATABLE :: rlon_v(:,:)
  REAL(4),ALLOCATABLE :: rlat_v(:,:)
  REAL(4),ALLOCATABLE :: fcori(:,:)
  REAL(4),ALLOCATABLE :: hgt(:,:)
  REAL(4),ALLOCATABLE :: ps(:,:)
  REAL(4),ALLOCATABLE :: t2(:,:)
  REAL(4),ALLOCATABLE :: q2(:,:)
  REAL(4),ALLOCATABLE :: u(:,:,:)
  REAL(4),ALLOCATABLE :: v(:,:,:)
  REAL(4),ALLOCATABLE :: w(:,:,:)
  REAL(4),ALLOCATABLE :: prs(:,:,:)
  REAL(4),ALLOCATABLE :: pb(:,:,:)
  REAL(4),ALLOCATABLE :: theta(:,:,:)
  REAL(4),ALLOCATABLE :: t(:,:,:)
  REAL(4),ALLOCATABLE :: ph(:,:,:)
  REAL(4),ALLOCATABLE :: phb(:,:,:)
  REAL(4),ALLOCATABLE :: qv(:,:,:)
  REAL(4),ALLOCATABLE :: qc(:,:,:)
  REAL(4),ALLOCATABLE :: qr(:,:,:)
  REAL(4),ALLOCATABLE :: rainc(:,:)
  REAL(4),ALLOCATABLE :: rainnc(:,:)
!                                123456789012345
  CHARACTER(8) :: ncin        = 'input.nc'
  CHARACTER(9) :: grdout_cnst = 'const.grd'
  CHARACTER(9) :: grdout      = 'guess.grd'
!
! initialization
!
  CALL check_io(NF_OPEN(ncin,NF_NOWRITE,ncid))
  CALL check_io(NF_INQ_VARID(ncid,'U',varid))
  CALL check_io(NF_INQ_VARDIMID(ncid,varid,dimids))
  DO i=1,3
    CALL check_io(NF_INQ_DIMLEN(ncid,dimids(i),shape(i)))
  END DO
  nlon = shape(1)
  nlat = shape(2) + 1
  nlev = shape(3) + 1
  WRITE(6,'(A)') '*** grid information ***'
  WRITE(6,'(3(2X,A,I5))') 'nlon =',nlon,'nlat =',nlat,'nlev =',nlev
  reclength = 4 * nlon * nlat
!
! ALLOCATE valiables
!
  ALLOCATE(rlon(nlon,nlat))
  ALLOCATE(rlat(nlon,nlat))
  ALLOCATE(rlon_u(nlon,nlat))
  ALLOCATE(rlat_u(nlon,nlat))
  ALLOCATE(rlon_v(nlon,nlat))
  ALLOCATE(rlat_v(nlon,nlat))
  ALLOCATE(fcori(nlon,nlat))
  ALLOCATE(hgt(nlon,nlat))
  ALLOCATE(ps(nlon,nlat))
  ALLOCATE(t2(nlon,nlat))
  ALLOCATE(q2(nlon,nlat))
  ALLOCATE(u(nlon,nlat,nlev))
  ALLOCATE(v(nlon,nlat,nlev))
  ALLOCATE(w(nlon,nlat,nlev))
  ALLOCATE(prs(nlon,nlat,nlev))
  ALLOCATE(pb(nlon,nlat,nlev))
  ALLOCATE(t(nlon,nlat,nlev))
  ALLOCATE(theta(nlon,nlat,nlev))
  ALLOCATE(ph(nlon,nlat,nlev))
  ALLOCATE(phb(nlon,nlat,nlev))
  ALLOCATE(qv(nlon,nlat,nlev))
  ALLOCATE(qc(nlon,nlat,nlev))
  ALLOCATE(qr(nlon,nlat,nlev))
  ALLOCATE(rainc(nlon,nlat))
  ALLOCATE(rainnc(nlon,nlat))
!
! read netcdf file (constants)
!
  ALLOCATE(start(3),count(3))
  start = (/ 1,1,1 /)  
  WRITE(6,'(A)') '*** START : READ NETCDF (CNST) ***'
  !!! XLONG
  count = (/ nlon-1,nlat-1,1 /)
  rlon(1:nlon,1:nlat) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'XLONG',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,rlon(1:nlon-1,1:nlat-1)))
  !!! XLAT
  count = (/ nlon-1,nlat-1,1 /)
  rlat(1:nlon,1:nlat) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'XLAT',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,rlat(1:nlon-1,1:nlat-1)))
  !!! XLONG_U
  count = (/ nlon,nlat-1,1 /)
  rlon_u(1:nlon,1:nlat) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'XLONG_U',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,rlon_u(1:nlon,1:nlat-1)))
  !!! XLAT_U
  count = (/ nlon,nlat-1,1 /)
  rlat_u(1:nlon,1:nlat) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'XLAT_U',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,rlat_u(1:nlon,1:nlat-1)))
  !!! XLONG_V
  count = (/ nlon-1,nlat,1 /)
  rlon_v(1:nlon,1:nlat) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'XLONG_V',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,rlon_v(1:nlon-1,1:nlat)))
  !!! XLAT_V
  count = (/ nlon-1,nlat,1 /)
  rlat_v(1:nlon,1:nlat) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'XLAT_V',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,rlat_v(1:nlon-1,1:nlat)))
  !!! F
  count = (/ nlon-1,nlat-1,1 /)
  fcori(1:nlon,1:nlat) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'F',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,fcori(1:nlon-1,1:nlat-1)))
  !!! HGT
  count = (/ nlon-1,nlat-1,1 /)
  hgt(1:nlon,1:nlat) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'HGT',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,hgt(1:nlon-1,1:nlat-1)))
  WRITE(6,'(a)') '***  END  : READ NETCDF (CNST) ***'
  DEALLOCATE(start,count)
!
! read netcdf file (3-D valiables)  
!
  ALLOCATE(start(4),count(4))
  start = (/ 1,1,1,1 /)
  WRITE(6,'(a)') '*** START : READ NETCDF (3D-VAR) ***'
  !!! U
  count = (/ nlon,nlat-1,nlev-1,1 /)
  u(1:nlon,1:nlat,1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'U',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,u(1:nlon,1:nlat-1,1:nlev-1)))
  CALL monit_3d('U',nlon,nlat,nlev,u)
  !!! V
  count = (/ nlon-1,nlat,nlev-1,1 /)
  v(1:nlon,1:nlat,1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'V',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,v(1:nlon-1,1:nlat,1:nlev-1)))
  CALL monit_3d('V',nlon,nlat,nlev,v)
  !!! W
  count = (/ nlon-1,nlat-1,nlev,1 /)
  w(1:nlon,1:nlat,1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'W',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,w(1:nlon-1,1:nlat-1,1:nlev)))
  CALL monit_3d('W',nlon,nlat,nlev,w)
  !!! P
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  prs(1:nlon,1:nlat,1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'P',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,prs(1:nlon-1,1:nlat-1,1:nlev-1)))
  pb(1:nlon,1:nlat,1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'PB',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,pb(1:nlon-1,1:nlat-1,1:nlev-1)))
  DO k=1,nlev-1
    DO j=1,nlat-1
      DO i=1,nlon-1
        prs(i,j,k) = prs(i,j,k) + pb(i,j,k)
      END DO
    END DO
  END DO
  CALL monit_3d('P',nlon,nlat,nlev,prs)
  !!! T
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  t(1:nlon,1:nlat,1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'T',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,theta(1:nlon-1,1:nlat-1,1:nlev-1)))
  theta(1:nlon-1,1:nlat-1,1:nlev-1) = theta(1:nlon-1,1:nlat-1,1:nlev-1) + t0
  CALL monit_3d('THETA',nlon,nlat,nlev,theta)
  DO k=1,nlev-1
    DO j=1,nlat-1
      DO i=1,nlon-1
        t(i,j,k) = theta(i,j,k) * (prs(i,j,k) / p0) ** (rd / cp)
      END DO
    END DO
  END DO
  CALL monit_3d('T',nlon,nlat,nlev,t)
  !!! PH
  count = (/ nlon-1,nlat-1,nlev,1 /)
  ph(1:nlon,1:nlat,1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'PH',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,ph(1:nlon-1,1:nlat-1,1:nlev)))
  phb(1:nlon,1:nlat,1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'PHB',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,phb(1:nlon-1,1:nlat-1,1:nlev)))
  DO k=1,nlev
    DO j=1,nlat-1
      DO i=1,nlon-1
        ph(i,j,k) = ph(i,j,k) + phb(i,j,k)
      END DO
    END DO
  END DO
  CALL monit_3d('PH',nlon,nlat,nlev,ph)
  !!! QVAPOR
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  qv(1:nlon,1:nlat,1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'QVAPOR',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,qv(1:nlon-1,1:nlat-1,1:nlev)))
  CALL monit_3d('QVAPOR',nlon,nlat,nlev,qv)
  !!! QCLOUD
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  qc(1:nlon,1:nlat,1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'QCLOUD',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,qc(1:nlon-1,1:nlat-1,1:nlev)))
  CALL monit_3d('QCLOUD',nlon,nlat,nlev,qc)
  !!! QRAIN
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  qr(1:nlon,1:nlat,1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'QRAIN',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,qr(1:nlon-1,1:nlat-1,1:nlev)))
  CALL monit_3d('QRAIN',nlon,nlat,nlev,qr)
  WRITE(6,'(a)') '***  END  : READ NETCDF (3D-VAR) ***'
  DEALLOCATE(start,count)
!
! read netcdf file (2-D valiables)  
!
  ALLOCATE(start(3),count(3))
  start = (/ 1,1,1 /)
  WRITE(6,'(a)') '*** START : READ NETCDF (2D-VAR) ***'
  !!! PSFC
  count = (/ nlon-1,nlat-1,1 /)
  ps(1:nlon,1:nlat) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'PSFC',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,ps(1:nlon-1,1:nlat-1)))
  !!! T2
  count = (/ nlon-1,nlat-1,1 /)
  t2(1:nlon,1:nlat) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'T2',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,t2(1:nlon-1,1:nlat-1)))
  !!! Q2
  count = (/ nlon-1,nlat-1,1 /)
  q2(1:nlon,1:nlat) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'Q2',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,q2(1:nlon-1,1:nlat-1)))
  !!! RAINC
  count = (/ nlon-1,nlat-1,1 /)
  rainc(1:nlon,1:nlat) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'RAINC',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,rainc(1:nlon-1,1:nlat-1)))
  !!! RAINNC
  count = (/ nlon-1,nlat-1,1 /)
  rainnc(1:nlon,1:nlat) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'RAINNC',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,rainnc(1:nlon-1,1:nlat-1)))
  !!! RAIN = RAINC + RAINNC
  rainc(1:nlon-1,1:nlat-1) = rainc(1:nlon-1,1:nlat-1) + rainnc(1:nlon-1,1:nlat-1)
  WRITE(6,'(a)') '***  END  : READ NETCDF (2D-VAR) ***'
  DEALLOCATE(start,count)
  CALL check_io(NF_CLOSE(ncid))
!
! output grads cntl file
!
  OPEN(imt_grd_cnst,FILE=grdout_cnst,FORM='unformatted',ACCESS='direct',RECL=reclength)
  irec = 1
  !!! Lon,Lat
  WRITE(imt_grd_cnst,REC=irec) ((rlon(i,j),i=1,nlon),j=1,nlat)
  irec = irec + 1
  WRITE(imt_grd_cnst,REC=irec) ((rlat(i,j),i=1,nlon),j=1,nlat)
  irec = irec + 1
  WRITE(imt_grd_cnst,REC=irec) ((rlon_u(i,j),i=1,nlon),j=1,nlat)
  irec = irec + 1
  WRITE(imt_grd_cnst,REC=irec) ((rlat_u(i,j),i=1,nlon),j=1,nlat)
  irec = irec + 1
  WRITE(imt_grd_cnst,REC=irec) ((rlon_v(i,j),i=1,nlon),j=1,nlat)
  irec = irec + 1
  WRITE(imt_grd_cnst,REC=irec) ((rlat_v(i,j),i=1,nlon),j=1,nlat)
  irec = irec + 1
  !!! Coriolis sine latitude term
  WRITE(imt_grd_cnst,REC=irec) ((fcori(i,j),i=1,nlon),j=1,nlat)
  irec = irec + 1
  !!! Terrain height
  WRITE(imt_grd_cnst,REC=irec) ((hgt(i,j),i=1,nlon),j=1,nlat)
  irec = irec + 1
  CLOSE(imt_grd_cnst)
!
! output grads var file
!
  OPEN(imt_grd,FILE=grdout,FORM='unformatted',ACCESS='direct',RECL=reclength)
  irec = 1
  !!! U
  DO k=1,nlev
    WRITE(imt_grd,REC=irec) ((u(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! V
  DO k=1,nlev
    WRITE(imt_grd,REC=irec) ((v(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! W
  DO k=1,nlev
    WRITE(imt_grd,REC=irec) ((w(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! T
  DO k=1,nlev
    WRITE(imt_grd,REC=irec) ((t(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! P
  DO k=1,nlev
    WRITE(imt_grd,REC=irec) ((prs(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! PH
  DO k=1,nlev
    WRITE(imt_grd,REC=irec) ((ph(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! QVAPOR
  DO k=1,nlev
    WRITE(imt_grd,REC=irec) ((qv(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! QCLOUD
  DO k=1,nlev
    WRITE(imt_grd,REC=irec) ((qc(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! QRAIN
  DO k=1,nlev
    WRITE(imt_grd,REC=irec) ((qr(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! PSFC
  WRITE(imt_grd,REC=irec) ((ps(i,j),i=1,nlon),j=1,nlat)
  irec = irec + 1
  !!! T2
  WRITE(imt_grd,REC=irec) ((t2(i,j),i=1,nlon),j=1,nlat)
  irec = irec + 1
  !!! Q2
  WRITE(imt_grd,REC=irec) ((q2(i,j),i=1,nlon),j=1,nlat)
  irec = irec + 1
  !!! RAIN
  WRITE(imt_grd,REC=irec) ((rainc(i,j),i=1,nlon),j=1,nlat)
  irec = irec + 1
  CLOSE(imt_grd)

  STOP
END PROGRAM nc2grd
!-----------------------------------------------------------------------
! Check the status of netcdf io
!-----------------------------------------------------------------------
SUBROUTINE check_io(status)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'  
  INTEGER(4),INTENT(IN) :: status
  
  IF(status /= nf_noerr) THEN
    WRITE(6,*) TRIM(nf_strerror(status))
    STOP 10
  ENDIF
  
  RETURN
END SUBROUTINE check_io
!-----------------------------------------------------------------------
! Monit 3-D variables
!-----------------------------------------------------------------------
SUBROUTINE monit_3d(elem,nlon,nlat,nlev,var)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: elem
  INTEGER(4),INTENT(IN) :: nlon,nlat,nlev
  REAL(4),INTENT(IN) :: var(nlon,nlat,nlev)
  INTEGER(4) :: ied,jed,ked,k
  
  IF(elem == 'U') THEN
    ied = nlon
    jed = nlat - 1
    ked = nlev - 1
  ELSE IF(elem == 'V') THEN
    ied = nlon - 1
    jed = nlat
    ked = nlev - 1
  ELSE IF(elem == 'W') THEN
    ied = nlon - 1
    jed = nlat - 1
    ked = nlev
  ELSE
    ied = nlon - 1
    jed = nlat - 1
    ked = nlev - 1
  END IF

  WRITE(6,'(3a)') ' === minmax monitor of ',elem,' ==='
  IF(elem == 'P') THEN
    DO k=1,ked
      WRITE(6,'(i5,2f12.5)') k,minval(var(1:ied,1:jed,k)) / 100.0,maxval(var(1:ied,1:jed,k))/ 100.0
    END DO
  ELSE
    DO k=1,ked
      WRITE(6,'(i5,2f12.5)') k,minval(var(1:ied,1:jed,k)),maxval(var(1:ied,1:jed,k))
    END DO
  END IF
    
  RETURN
END SUBROUTINE monit_3d
