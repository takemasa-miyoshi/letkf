MODULE common_wrf
!=======================================================================
!
! [PURPOSE:] Common Information for WRF
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   07/14/2010 Takemasa Miyoshi  modified for WRF
!
!=======================================================================
!$USE OMP_LIB
  USE common
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  INTEGER,PARAMETER :: dx=60.0d3 ! grid spacing [m]
  INTEGER,PARAMETER :: nlon=137
  INTEGER,PARAMETER :: nlat=109
  INTEGER,PARAMETER :: nlev=40
  INTEGER,PARAMETER :: nv3d=8 ! u,v,w,t,p,qv,qc,qr
  INTEGER,PARAMETER :: nv2d=1 ! rain
  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_w=3
  INTEGER,PARAMETER :: iv3d_t=4
  INTEGER,PARAMETER :: iv3d_p=5
  INTEGER,PARAMETER :: iv3d_qv=6
  INTEGER,PARAMETER :: iv3d_qc=7
  INTEGER,PARAMETER :: iv3d_qr=8
  INTEGER,PARAMETER :: iv2d_rain=1
  INTEGER,PARAMETER :: nij0=nlon*nlat
  INTEGER,PARAMETER :: nlevall=nlev*nv3d+nv2d
  INTEGER,PARAMETER :: ngpv=nij0*nlevall
  REAL(r_size),SAVE :: lon(nlon,nlat)
  REAL(r_size),SAVE :: lat(nlon,nlat)
  REAL(r_size),SAVE :: lonu(nlon,nlat)
  REAL(r_size),SAVE :: latu(nlon,nlat)
  REAL(r_size),SAVE :: lonv(nlon,nlat)
  REAL(r_size),SAVE :: latv(nlon,nlat)
  REAL(r_size),SAVE :: fcori(nlon,nlat)
  REAL(r_size),SAVE :: phi0(nlon,nlat)
  CHARACTER(4),SAVE :: element(nv3d+nv2d)

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_wrf
  IMPLICIT NONE
  REAL(r_sngl) :: buf4(nlon,nlat,8)
  INTEGER :: i,j,iunit,iolen,n,irec

  WRITE(6,'(A)') 'Hello from set_common_wrf'
  !
  ! Elements
  !
  element(iv3d_u)  = 'U   '
  element(iv3d_v)  = 'V   '
  element(iv3d_w)  = 'W   '
  element(iv3d_t)  = 'T   '
  element(iv3d_p)  = 'P   '
  element(iv3d_qv) = 'QV  '
  element(iv3d_qc) = 'QC  '
  element(iv3d_qr) = 'QR  '
  element(nv3d+iv2d_rain) = 'RAIN'
  !
  ! Lon, Lat, F, phi0
  !
  iunit=12
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE='const.grd',FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)
  irec=1
  DO n=1,8
    READ(iunit,REC=irec) ((buf4(i,j,n),i=1,nlon),j=1,nlat)
    irec=irec+1
  END DO
  CLOSE(iunit)
  lon   = REAL(buf4(:,:,1),r_size)
  lat   = REAL(buf4(:,:,2),r_size)
  lonu  = REAL(buf4(:,:,3),r_size)
  latu  = REAL(buf4(:,:,4),r_size)
  lonv  = REAL(buf4(:,:,5),r_size)
  latv  = REAL(buf4(:,:,6),r_size)
  fcori = REAL(buf4(:,:,7),r_size)
  phi0  = REAL(buf4(:,:,8),r_size)

  RETURN
END SUBROUTINE set_common_wrf
!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------
!-- Read a grid file ---------------------------------------------------
SUBROUTINE read_grd(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl) :: buf4(nlon,nlat)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,irec

  iunit=11
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3d
    DO k=1,nlev
      READ(iunit,REC=irec) buf4
      irec = irec + 1
      v3d(:,:,k,n) = REAL(buf4,r_size)
    END DO
  END DO

  DO n=1,nv2d
    READ(iunit,REC=irec) buf4
    irec = irec + 1
    v2d(:,:,n) = REAL(buf4,r_size)
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE read_grd

SUBROUTINE read_grd4(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec

  iunit=11
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3d
    DO k=1,nlev
      READ(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2d
    READ(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE read_grd4
!-- Write a grid file -------------------------------------------------
SUBROUTINE write_grd(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl) :: buf4(nlon,nlat)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,irec

  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3d
    DO k=1,nlev
      buf4 = REAL(v3d(:,:,k,n),r_sngl)
      WRITE(iunit,REC=irec) buf4
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2d
    buf4 = REAL(v2d(:,:,n),r_sngl)
    WRITE(iunit,REC=irec) buf4
    irec = irec + 1
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_grd

SUBROUTINE write_grd4(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec

  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3d
    DO k=1,nlev
      WRITE(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2d
    WRITE(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_grd4
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

END MODULE common_wrf
