PROGRAM init_merge
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER(4),PARAMETER :: imt_grd = 50
  INTEGER(4) :: nlon,nlat,nlev
  INTEGER(4) :: ncid,varid,dimids(3),shape(3)
  INTEGER(4) :: reclength,irec
  INTEGER(4) :: i,j,k
  INTEGER(4),ALLOCATABLE :: start(:),count(:)
  REAL(4),PARAMETER :: rd = 287.05
  REAL(4),PARAMETER :: cp = 7.0 / 2.0 * rd
  REAL(4),PARAMETER :: p0 = 1.0e+5
  REAL(4),PARAMETER :: t0 = 300.0
  REAL(4),ALLOCATABLE :: ps(:,:)
  REAL(4),ALLOCATABLE :: u(:,:,:)
  REAL(4),ALLOCATABLE :: v(:,:,:)
  REAL(4),ALLOCATABLE :: w(:,:,:)
  REAL(4),ALLOCATABLE :: prs(:,:,:)
  REAL(4),ALLOCATABLE :: pb(:,:,:)
  REAL(4),ALLOCATABLE :: theta(:,:,:)
  REAL(4),ALLOCATABLE :: t(:,:,:)
  REAL(4),ALLOCATABLE :: ph(:,:,:)
  REAL(4),ALLOCATABLE :: phb(:,:,:)
  REAL(4),ALLOCATABLE :: ph_gs(:,:,:)
  REAL(4),ALLOCATABLE :: qv(:,:,:)
  REAL(4),ALLOCATABLE :: qc(:,:,:)
  REAL(4),ALLOCATABLE :: qr(:,:,:)
  REAL(4),ALLOCATABLE :: t2(:,:)
  REAL(4),ALLOCATABLE :: q2(:,:)
!  REAL(4),ALLOCATABLE :: rainc(:,:)
!  REAL(4),ALLOCATABLE :: rainnc(:,:)
!                          1234567890
  CHARACTER(8) :: ncin  = 'input.nc'
  CHARACTER(8) :: grdin = 'anal.grd'
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
  ALLOCATE(ps(nlon,nlat))
  ALLOCATE(u(nlon,nlat,nlev))
  ALLOCATE(v(nlon,nlat,nlev))
  ALLOCATE(w(nlon,nlat,nlev))
  ALLOCATE(prs(nlon,nlat,nlev))
  ALLOCATE(pb(nlon,nlat,nlev))
  ALLOCATE(t(nlon,nlat,nlev))
  ALLOCATE(theta(nlon,nlat,nlev))
  ALLOCATE(ph(nlon,nlat,nlev))
  ALLOCATE(phb(nlon,nlat,nlev))
  ALLOCATE(ph_gs(nlon,nlat,nlev))
  ALLOCATE(qv(nlon,nlat,nlev))
  ALLOCATE(qc(nlon,nlat,nlev))
  ALLOCATE(qr(nlon,nlat,nlev))
  ALLOCATE(t2(nlon,nlat))
  ALLOCATE(q2(nlon,nlat))
!  ALLOCATE(rainc(nlon,nlat))
!  ALLOCATE(rainnc(nlon,nlat))

!
! READ netcdf file (3-D valiables)
!
  ALLOCATE(start(4),count(4))
  start = (/ 1,1,1,1 /)
  WRITE(6,'(a)') '*** START : READ NETCDF (3D-VAR) ***'
  !!! PB
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  pb(1:nlon,1:nlat,1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'PB',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,pb(1:nlon-1,1:nlat-1,1:nlev-1)))
  !!! PHB
  count = (/ nlon-1,nlat-1,nlev,1 /)
  phb(1:nlon,1:nlat,1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'PHB',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,phb(1:nlon-1,1:nlat-1,1:nlev)))
  !!! PH
  count = (/ nlon-1,nlat-1,nlev,1 /)
  ph_gs(1:nlon,1:nlat,1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'PH',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,ph_gs(1:nlon-1,1:nlat-1,1:nlev)))
  WRITE(6,'(a)') '***  END  : READ NETCDF (3D-VAR) ***'
  DEALLOCATE(start,count)
  CALL check_io(NF_CLOSE(ncid))
!
! READ grads file
!
  OPEN(imt_grd,FILE=grdin,FORM='unformatted',ACCESS='direct',RECL=reclength)
  irec = 1
  !!! U
  DO k=1,nlev
    READ(imt_grd,rec=irec) ((u(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1  
  END DO
  !!! V
  DO k=1,nlev
    READ(imt_grd,rec=irec) ((v(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! W
  DO k=1,nlev
    READ(imt_grd,rec=irec) ((w(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! T
  DO k=1,nlev
    READ(imt_grd,rec=irec) ((t(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! P
  DO k=1,nlev
    READ(imt_grd,rec=irec) ((prs(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! PH
  DO k=1,nlev
    READ(imt_grd,rec=irec) ((ph(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! QVAPOR
  DO k=1,nlev
    READ(imt_grd,rec=irec) ((qv(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! QCLOUD
  DO k=1,nlev
    READ(imt_grd,rec=irec) ((qc(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! QRAIN
  DO k=1,nlev
    READ(imt_grd,rec=irec) ((qr(i,j,k),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO
  !!! PSFC
  READ(imt_grd,rec=irec) ((ps(i,j),i=1,nlon),j=1,nlat)
  !!! T2
  READ(imt_grd,rec=irec) ((t2(i,j),i=1,nlon),j=1,nlat)
  !!! Q2
  READ(imt_grd,rec=irec) ((q2(i,j),i=1,nlon),j=1,nlat)
!  !!! RAIN
!  READ(imt_grd,rec=irec) ((rainc(i,j),i=1,nlon),j=1,nlat)
  CLOSE(imt_grd)
!
! check PH
!
  WRITE(6,*) '@@@ check PH @@@'
  DO k=1,nlev
    DO j=1,nlat-1
      DO i=1,nlon-1
        ph_gs(i,j,k) = ph_gs(i,j,k) + phb(i,j,k)
      END DO
    END DO
  END DO
  DO k=1,nlev
    WRITE(6,*) k,ph_gs(10,70,k),ph(10,70,k),phb(10,70,k)
  END DO
!
! convert P,T,PH
!
 DO k=1,nlev-1
    DO j=1,nlat-1
      DO i=1,nlon-1
        theta(i,j,k) = t(i,j,k) *(p0 /  prs(i,j,k)) ** (rd / cp) - t0
        prs(i,j,k) = prs(i,j,k) - pb(i,j,k)
        ph(i,j,k) = ph(i,j,k) - phb(i,j,k)
      END DO
    END DO
  END DO
  DO j=1,nlat-1
    DO i=1,nlon-1
      ph(i,j,nlev) = ph(i,j,nlev) - phb(i,j,nlev)
    END DO
  END DO
!
! WRITE netcdf file (2-D valiables)
!
  CALL check_io(NF_OPEN(ncin,NF_WRITE,ncid))
  ALLOCATE(start(3),count(3))
  start = (/ 1,1,1 /)
  WRITE(6,'(a)') '*** START : WRITE NETCDF (2D-VAR) ***'
  !!! PSFC
  count = (/ nlon-1,nlat-1,1 /)
  CALL check_io(NF_INQ_VARID(ncid,'PSFC',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,ps(1:nlon-1,1:nlat-1)))
  !!! T2
  count = (/ nlon-1,nlat-1,1 /)
  CALL check_io(NF_INQ_VARID(ncid,'T2',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,t2(1:nlon-1,1:nlat-1)))
  !!! Q2
  count = (/ nlon-1,nlat-1,1 /)
  CALL check_io(NF_INQ_VARID(ncid,'Q2',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,q2(1:nlon-1,1:nlat-1)))
  WRITE(6,'(a)') '***  END  : WRITE NETCDF (2D-VAR) ***'
  DEALLOCATE(start,count)
!
! WRITE netcdf file (3-D valiables)
!
  ALLOCATE(start(4),count(4))
  start = (/ 1,1,1,1 /)
  WRITE(6,'(a)') '*** START : WRITE NETCDF (3D-VAR) ***'
  !!! U
  count = (/ nlon,nlat-1,nlev-1,1 /)
  CALL check_io(NF_INQ_VARID(ncid,'U',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,u(1:nlon,1:nlat-1,1:nlev-1)))
  CALL monit_3d('U',nlon,nlat,nlev,u)
  !!! V
  count = (/ nlon-1,nlat,nlev-1,1 /)
  CALL check_io(NF_INQ_VARID(ncid,'V',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,v(1:nlon-1,1:nlat,1:nlev-1)))
  CALL monit_3d('V',nlon,nlat,nlev,v)
  !!! W
  count = (/ nlon-1,nlat-1,nlev,1 /)
  CALL check_io(NF_INQ_VARID(ncid,'W',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,w(1:nlon-1,1:nlat-1,1:nlev)))
  CALL monit_3d('W',nlon,nlat,nlev,w)
  !!! P
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  CALL check_io(NF_INQ_VARID(ncid,'P',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,prs(1:nlon-1,1:nlat-1,1:nlev-1)))  
  CALL monit_3d('P',nlon,nlat,nlev,prs)
  !!! T
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  CALL check_io(NF_INQ_VARID(ncid,'T',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,theta(1:nlon-1,1:nlat-1,1:nlev-1))) 
  CALL monit_3d('T',nlon,nlat,nlev,theta)
  !!! PH
  count = (/ nlon-1,nlat-1,nlev,1 /)
  CALL check_io(NF_INQ_VARID(ncid,'PH',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,ph(1:nlon-1,1:nlat-1,1:nlev)))
  CALL monit_3d('PH',nlon,nlat,nlev,ph)
  !!! QVAPOR
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  CALL check_io(NF_INQ_VARID(ncid,'QVAPOR',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,qv(1:nlon-1,1:nlat-1,1:nlev)))
  CALL monit_3d('QVAPOR',nlon,nlat,nlev,qv)
  !!! QCLOUD
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  CALL check_io(NF_INQ_VARID(ncid,'QCLOUD',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,qc(1:nlon-1,1:nlat-1,1:nlev)))
  CALL monit_3d('QCLOUD',nlon,nlat,nlev,qc)
  !!! QRAIN
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  CALL check_io(NF_INQ_VARID(ncid,'QRAIN',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,qr(1:nlon-1,1:nlat-1,1:nlev)))
  CALL monit_3d('QRAIN',nlon,nlat,nlev,qr)
  WRITE(6,'(a)') '***  END  : WRITE NETCDF (3D-VAR) ***'
  DEALLOCATE(start,count)
  CALL check_io(NF_CLOSE(ncid))

  STOP
END PROGRAM init_merge
!-----------------------------------------------------------------------
! Check the status of netcdf io
!-----------------------------------------------------------------------
SUBROUTINE check_io(status)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'  
  INTEGER(4),INTENT(IN) :: status
  
  IF(status /= nf_noerr) THEN
    WRITE(6,*) trim(nf_strerror(status))
    STOP 10
  ENDIF
  
  RETURN
END SUBROUTINE check_io
!-----------------------------------------------------------------------
! Monit 3-D variables
!-----------------------------------------------------------------------
SUBROUTINE monit_3d(elem,nlon,nlat,nlev,var)
  IMPLICIT NONE
  CHARACTER(*),INTENT(in) :: elem
  INTEGER(4),INTENT(in) :: nlon,nlat,nlev
  REAL(4),INTENT(in) :: var(nlon,nlat,nlev)
  INTEGER(4) :: ied,jed,ked,k
  
  IF(elem == 'U') THEN
    ied = nlon
    jed = nlat - 1
    ked = nlev - 1
  ELSE IF(elem == 'V') THEN
    ied = nlon - 1
    jed = nlat
    ked = nlev - 1
  ELSE IF(elem == 'W' .or. elem == 'PH') THEN
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
