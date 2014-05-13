MODULE common_gfs_pres
!=======================================================================
!
! [PURPOSE:] Common Information for GFS pressure-level I/O
!
! [HISTORY:]
!   12/17/2012 Guo-Yuan Lien     created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_gfs, ONLY: nlon,nlat,nlev,nij0
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  INTEGER,PARAMETER :: nlevp = 12
  REAL(r_sngl),PARAMETER :: levp(nlevp) = &
        (/100000., 92500., 85000., 70000., 50000., 30000., &
           20000., 10000.,  5000.,  2000.,  1000.,   500./)

  INTEGER,PARAMETER :: nv3dp=8  ! 3D variables: u,v,t,q,rh,qc,gph,ozone
  INTEGER,PARAMETER :: nv2dp=10 ! 2D variables: ps,slp,orog,slmsk,tsea,u10m,v10m,t2m,q2m,tprcp
  INTEGER,PARAMETER :: iv3dp_u=1
  INTEGER,PARAMETER :: iv3dp_v=2
  INTEGER,PARAMETER :: iv3dp_t=3
  INTEGER,PARAMETER :: iv3dp_q=4
  INTEGER,PARAMETER :: iv3dp_rh=5
  INTEGER,PARAMETER :: iv3dp_qc=6
  INTEGER,PARAMETER :: iv3dp_gph=7
  INTEGER,PARAMETER :: iv3dp_ozone=8
  INTEGER,PARAMETER :: iv2dp_ps=1
  INTEGER,PARAMETER :: iv2dp_slp=2
  INTEGER,PARAMETER :: iv2dp_orog=3
  INTEGER,PARAMETER :: iv2dp_slmsk=4
  INTEGER,PARAMETER :: iv2dp_tsea=5
  INTEGER,PARAMETER :: iv2dp_u10m=6
  INTEGER,PARAMETER :: iv2dp_v10m=7
  INTEGER,PARAMETER :: iv2dp_t2m=8
  INTEGER,PARAMETER :: iv2dp_q2m=9
  INTEGER,PARAMETER :: iv2dp_tprcp=10
  INTEGER,PARAMETER :: nlevallp=nlevp*nv3dp+nv2dp
  INTEGER,PARAMETER :: ngpvp=nij0*nlevallp

CONTAINS
!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------
!-- Read a grid file ---------------------------------------------------
SUBROUTINE read_grdp(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlevp,nv3dp)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2dp)
  REAL(r_sngl) :: buf4(nlon,nlat)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,irec

  iunit=11
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3dp
    DO k=1,nlevp
      READ(iunit,REC=irec) buf4
      irec = irec + 1
      v3d(:,:,k,n) = REAL(buf4,r_size)
    END DO
  END DO

  DO n=1,nv2dp
    READ(iunit,REC=irec) buf4
    irec = irec + 1
    v2d(:,:,n) = REAL(buf4,r_size)
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE read_grdp

SUBROUTINE read_grd4p(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlevp,nv3dp)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2dp)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec

  iunit=11
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3dp
    DO k=1,nlevp
      READ(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2dp
    READ(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE read_grd4p
!-- Write a grid file -------------------------------------------------
SUBROUTINE write_grdp(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlevp,nv3dp)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2dp)
  REAL(r_sngl) :: buf4(nlon,nlat)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,irec

  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3dp
    DO k=1,nlevp
      buf4 = REAL(v3d(:,:,k,n),r_sngl)
      WRITE(iunit,REC=irec) buf4
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2dp
    buf4 = REAL(v2d(:,:,n),r_sngl)
    WRITE(iunit,REC=irec) buf4
    irec = irec + 1
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_grdp

SUBROUTINE write_grd4p(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlevp,nv3dp)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2dp)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec

  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3dp
    DO k=1,nlevp
      WRITE(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2dp
    WRITE(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_grd4p

END MODULE common_gfs_pres
