MODULE common_gfs
!=======================================================================
!
! [PURPOSE:] Common Information for GFS
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   10/03/2012 Guo-Yuan Lien     modified for GFS model
!   07/01/2013 Daisuke Hotta     added the weight for global averages
!   12/19/2013 Guo-Yuan Lien     modified the weight for global averages
!
!=======================================================================
!$USE OMP_LIB
  USE common
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  INTEGER,PARAMETER :: nlon=192
  INTEGER,PARAMETER :: nlat=94
  INTEGER,PARAMETER :: nlev=64
  INTEGER,PARAMETER :: idrt=4  ! 4: Gaussian grid
  INTEGER,PARAMETER :: nv3dx=6 ! 3D extended variables: [u,v,t,q,qc],p
  INTEGER,PARAMETER :: nv3d=5  ! 3D state variables
  INTEGER,PARAMETER :: nv2dx=6 ! 2D extended variables: [ps],orog,f10m,t2m,q2m,tprcp
  INTEGER,PARAMETER :: nv2d=1  ! 2D state variables
  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_t=3
  INTEGER,PARAMETER :: iv3d_q=4
  INTEGER,PARAMETER :: iv3d_qc=5
  INTEGER,PARAMETER :: iv3d_p=6
  INTEGER,PARAMETER :: iv2d_ps=1
  INTEGER,PARAMETER :: iv2d_orog=2
  INTEGER,PARAMETER :: iv2d_f10m=3
  INTEGER,PARAMETER :: iv2d_t2m=4
  INTEGER,PARAMETER :: iv2d_q2m=5
  INTEGER,PARAMETER :: iv2d_tprcp=6
  INTEGER,PARAMETER :: nij0=nlon*nlat
  INTEGER,PARAMETER :: nlevallx=nlev*nv3dx+nv2dx
  INTEGER,PARAMETER :: nlevall=nlev*nv3d+nv2d
  INTEGER,PARAMETER :: ngpvx=nij0*nlevallx
  INTEGER,PARAMETER :: ngpv=nij0*nlevall
  INTEGER,PARAMETER :: gfs_jcap=62
  INTEGER,PARAMETER :: gfs_ntrac=3
  INTEGER,PARAMETER :: gfs_nvcoord=3
  INTEGER,PARAMETER :: gfs_idsl=2
  INTEGER,PARAMETER :: gfs_idvc=3
  REAL(r_size),PARAMETER :: gfs_vcoord(nlev+1,gfs_nvcoord)=RESHAPE( &
  (/0.000000,0.000000,0.5750000,5.741000,21.51600,55.71200,116.8990,214.0150, &
    356.2230,552.7200,812.4890,1143.988,1554.789,2051.150,2637.553,3316.217, &
    4086.614,4945.029,5884.206,6893.117,7956.908,9057.051,10171.71,11276.35, &
    12344.49,13348.67,14261.43,15056.34,15708.89,16197.32,16503.14,16611.60, &
    16511.74,16197.97,15683.49,14993.07,14154.32,13197.07,12152.94,11054.85, &
    9936.614,8832.537,7777.150,6804.874,5937.050,5167.146,4485.493,3883.052, &
    3351.460,2883.038,2470.788,2108.366,1790.051,1510.711,1265.752,1051.080, &
    863.0580,698.4570,554.4240,428.4340,318.2660,221.9580,137.7900,64.24700, &
    0.000000,1.000000,0.9946712,0.9886266,0.9817423,0.9738676,0.9648276,0.9544341, &
    0.9424911,0.9287973,0.9131510,0.8953550,0.8752236,0.8525907,0.8273188,0.7993097, &
    0.7685147,0.7349452,0.6986829,0.6598870,0.6187996,0.5757467,0.5311348,0.4854433, &
    0.4392108,0.3930182,0.3474685,0.3031641,0.2606854,0.2205702,0.1832962,0.1492688, &
    0.1188122,9.2166908E-02,6.9474578E-02,5.0646842E-02,3.5441618E-02,2.3555880E-02,1.4637120E-02,8.2940198E-03, &
    4.1067102E-03,1.6359100E-03,4.3106001E-04,3.6969999E-05,0.000000,0.000000,0.000000,0.000000, &
    0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
    0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
    0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
    0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
    0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
    0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
    0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
    0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
    0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
    0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
    0.000000,0.000000,0.000000/),(/nlev+1,gfs_nvcoord/) )
  REAL(r_size),SAVE :: lon(nlon)
  REAL(r_size),SAVE :: lat(nlat)
  REAL(r_size),SAVE :: dx(nlat)
  REAL(r_size),SAVE :: dy(nlat)
  REAL(r_size),SAVE :: dy2(nlat)
  REAL(r_size),SAVE :: fcori(nlat)
  REAL(r_size),SAVE :: wg(nlon,nlat)
  CHARACTER(4),SAVE :: element(nv3d+nv2d)

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_gfs
  IMPLICIT NONE
  REAL(r_sngl) :: slat(nlat), wlat(nlat)
  REAL(r_size) :: totalwg, wgtmp, latm1, latm2
  INTEGER :: i,j

!  WRITE(6,'(A)') 'Hello from set_common_gfs'
  !
  ! Elements
  !
  element(iv3d_u)  = 'U   '
  element(iv3d_v)  = 'V   '
  element(iv3d_t)  = 'T   '
  element(iv3d_q)  = 'Q   '
  element(iv3d_qc) = 'QC  '
  element(nv3d+iv2d_ps) = 'PS  '
  !
  ! Lon, Lat
  !
!$OMP PARALLEL DO PRIVATE(i)
  DO i=1,nlon
    lon(i) = 360.d0/nlon*(i-1)
  END DO
!$OMP END PARALLEL DO
  CALL SPLAT(idrt,nlat,slat,wlat)
  do j=1,nlat
    lat(j) = 180.d0/pi*asin(slat(nlat-j+1))
  end do
  !
  ! dx and dy
  !
!$OMP PARALLEL
!$OMP WORKSHARE
  dx(:) = 2.0d0 * pi * re * cos(lat(:) * pi / 180.0d0) / REAL(nlon,r_size)
!$OMP END WORKSHARE

!$OMP DO
  DO i=1,nlat-1
    dy(i) = 2.0d0 * pi * re * (lat(i+1) - lat(i)) / 360.0d0
  END DO
!$OMP END DO
!$OMP END PARALLEL
  dy(nlat) = 2.0d0 * pi * re * (90.0d0 - lat(nlat)) / 180.0d0

!$OMP PARALLEL DO
  DO i=2,nlat
    dy2(i) = (dy(i-1) + dy(i)) * 0.5d0
  END DO
!$OMP END PARALLEL DO
  dy2(1) = (dy(nlat) + dy(1)) * 0.5d0
  !
  ! Corioris parameter
  !
!$OMP PARALLEL WORKSHARE
  fcori(:) = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)
!$OMP END PARALLEL WORKSHARE
  !
  ! Weight for global average
  !
  totalwg = 0.0_r_size
  DO j=1,nlat
    if (j == 1) then
      latm1 = -0.5d0*pi !-90 degree
    else
      latm1 = 0.5d0*(lat(j-1) + lat(j))*pi/180.0d0
    end if
    if (j == nlat) then
      latm2 = 0.5d0*pi !90 degree
    else
      latm2 = 0.5d0*(lat(j) + lat(j+1))*pi/180.0d0
    end if
    wgtmp = abs(sin(latm2) - sin(latm1))
    wg(:,j) = wgtmp
    totalwg = totalwg + wgtmp * nlon
  END DO
  totalwg = 1.0_r_size / totalwg
  wg(:,:) = sqrt(wg(:,:) * totalwg)
  RETURN
END SUBROUTINE set_common_gfs
!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------
!-- Read a grid file ---------------------------------------------------
SUBROUTINE read_grd(filename,v3d,v2d,opt)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER,INTENT(IN) :: opt ! 0: read state variables from a file with only state variables
                            ! 1: read state variables from a file with extended variables
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

  if (opt == 1) irec = nv3dx*nlev+1
  DO n=1,nv2d
    READ(iunit,REC=irec) buf4
    irec = irec + 1
    v2d(:,:,n) = REAL(buf4,r_size)
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE read_grd

SUBROUTINE read_grdx(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3dx)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2dx)
  REAL(r_sngl) :: buf4(nlon,nlat)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,irec

  iunit=11
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3dx
    DO k=1,nlev
      READ(iunit,REC=irec) buf4
      irec = irec + 1
      v3d(:,:,k,n) = REAL(buf4,r_size)
    END DO
  END DO

  DO n=1,nv2dx
    READ(iunit,REC=irec) buf4
    irec = irec + 1
    v2d(:,:,n) = REAL(buf4,r_size)
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE read_grdx

SUBROUTINE read_grd4(filename,v3d,v2d,opt)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER,INTENT(IN) :: opt ! 0: read state variables from a file with only state variables
                            ! 1: read state variables from a file with extended variables
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

  if (opt == 1) irec = nv3dx*nlev+1
  DO n=1,nv2d
    READ(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE read_grd4

SUBROUTINE read_grd4x(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3dx)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2dx)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec

  iunit=11
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3dx
    DO k=1,nlev
      READ(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2dx
    READ(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE read_grd4x
!-- Write a grid file -------------------------------------------------
SUBROUTINE write_grd(filename,v3d,v2d,opt)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER,INTENT(IN) :: opt ! 0: read state variables from a file with only state variables
                            ! 1: read state variables from a file with extended variables
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

  if (opt == 1) irec = nv3dx*nlev+1
  DO n=1,nv2d
    buf4 = REAL(v2d(:,:,n),r_sngl)
    WRITE(iunit,REC=irec) buf4
    irec = irec + 1
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_grd

SUBROUTINE write_grdx(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3dx)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2dx)
  REAL(r_sngl) :: buf4(nlon,nlat)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,irec

  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3dx
    DO k=1,nlev
      buf4 = REAL(v3d(:,:,k,n),r_sngl)
      WRITE(iunit,REC=irec) buf4
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2dx
    buf4 = REAL(v2d(:,:,n),r_sngl)
    WRITE(iunit,REC=irec) buf4
    irec = irec + 1
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_grdx

SUBROUTINE write_grd4(filename,v3d,v2d,opt)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER,INTENT(IN) :: opt ! 0: read state variables from a file with only state variables
                            ! 1: read state variables from a file with extended variables
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

  if (opt == 1) irec = nv3dx*nlev+1
  DO n=1,nv2d
    WRITE(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_grd4

SUBROUTINE write_grd4x(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3dx)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2dx)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec

  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3dx
    DO k=1,nlev
      WRITE(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2dx
    WRITE(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_grd4x
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

END MODULE common_gfs
