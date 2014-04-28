MODULE common_afes
!=======================================================================
!
! [PURPOSE:] Common Information for AFES
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   09/25/2009 Takemasa Miyoshi  modified
!
!=======================================================================
!$USE OMP_LIB
  USE common
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  INTEGER,PARAMETER :: nlon=360
  INTEGER,PARAMETER :: nlat=nlon/2
  INTEGER,PARAMETER :: nlev=48
  INTEGER,PARAMETER :: nv3d=5 ! u,v,t,q,cw
  INTEGER,PARAMETER :: nv2d=1 ! ps
  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_t=3
  INTEGER,PARAMETER :: iv3d_q=4
  INTEGER,PARAMETER :: iv3d_cw=5
  INTEGER,PARAMETER :: iv2d_ps=1
  INTEGER,PARAMETER :: nij0=nlon*nlat
  INTEGER,PARAMETER :: nlevall=nlev*nv3d+nv2d
  INTEGER,PARAMETER :: ngpv=nij0*nlevall
  REAL(r_size),SAVE :: lon(nlon)
  REAL(r_size),SAVE :: lat(nlat)
  REAL(r_size),SAVE :: sig(nlev)
  REAL(r_size),SAVE :: dx(nlat)
  REAL(r_size),SAVE :: dy(nlat)
  REAL(r_size),SAVE :: dy2(nlat)
  REAL(r_size),SAVE :: fcori(nlat)
  REAL(r_size),SAVE :: phi0(nlon,nlat)
  CHARACTER(4),SAVE :: element(nv3d+nv2d)

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_afes
  IMPLICIT NONE
  INTEGER :: i,j

  WRITE(6,'(A)') 'Hello from set_common_afes'
  !
  ! Elements
  !
  element(iv3d_u) = 'U   '
  element(iv3d_v) = 'V   '
  element(iv3d_t) = 'T   '
  element(iv3d_q) = 'Q   '
  element(iv3d_cw)= 'CW  '
  element(nv3d+iv2d_ps) = 'PS  '
  !
  ! Lon
  !
!$OMP PARALLEL DO PRIVATE(i)
  DO i=1,nlon
    lon(i) = 360.d0/nlon*(i-1)
  END DO
!$OMP END PARALLEL DO
  !
  ! Lat
  !
  CALL setlat(lat,nlat)
  !
  ! Sigma levels
  !
  READ(41) sig
  !
  ! dx and dy
  !
!$OMP PARALLEL
!$OMP WORKSHARE
  dx(:) = 2.0d0 * pi * re * COS(lat(:) * pi / 180.0d0) / REAL(nlon,r_size)
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
  fcori(:) = 2.0d0 * r_omega * SIN(lat(:)*pi/180.0d0)
!$OMP END PARALLEL WORKSHARE
  !
  ! Surface geoptential (Read Orography file)
  !
  READ(21) phi0

  RETURN
END SUBROUTINE set_common_afes
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
! p_full
!-----------------------------------------------------------------------
SUBROUTINE calc_pfull(ix,jy,ps,p_full)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ix,jy
  REAL(r_size),INTENT(IN) :: ps(ix,jy)
  REAL(r_size),INTENT(OUT) :: p_full(ix,jy,nlev)
  INTEGER :: i,j,k

!$OMP PARALLEL DO PRIVATE(i,j,k)
  DO k=1,nlev
    DO j=1,jy
      DO i=1,ix
        p_full(i,j,k) = ps(i,j) * sig(k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  RETURN
END SUBROUTINE calc_pfull
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
!-----------------------------------------------------------------------
! set latitude
!-----------------------------------------------------------------------
  subroutine setlat(lat,nlat)
    implicit none
    real(kind=8), dimension(:), intent(out) :: lat
    integer(kind=4), intent(in) :: nlat
    real(kind=8), dimension(size(lat)) :: ww
    real(kind=8) :: pi

    call f_scgaus(lat, ww, nlat)
    pi=atan(1.d0)*4
    lat(:)=-asin(lat(:))/pi*180

  end subroutine setlat


!***********************************************************************
! following programs are from AFES code: leg_macro/fourier_newton_module.f90
!
! calculates Gaussian points and weights.
!
! Source: Based on Swartztrauber (2003)
! Author: T. Enomoto
! Usage:
!   Calculate Gaussian points and weights
!     subroutine f_scgaus(x, w, jMax)
!       x(jMax)  cos(Gaussian colatitudes) = sin(Gaussian latitudes)
!       w(jMax)  Gaussian weights
! History:
!   TE 27 Apr 2003  Fixed a bug in setting SH colatitudes
!   TE 28 Apr 2003  Ported for AFES
!   TE 24 Apr 2003  Implemented Fourier-Legendre formulation.
!***********************************************************************
  subroutine f_scgaus(x, w, jMax)
!
    implicit none
!
! returns sin(Gaussian latitudes) between 1 and -1
! and Gaussian weights.
! NB. Gaussian colatitudes are used during calculation
!
    real(kind=8), dimension(:), intent(out) :: x, w
    integer(kind=4), intent(in) :: jMax
!
    integer(kind=4) :: jMid
    integer(kind=4) :: l, j
    real(kind=8) :: guess, dpn, pi
    real(kind=8), dimension(0:jMax/2) :: an
!-----------------------------------------------------------------------
    pi = acos(-1.d0)
    jMid = jMax/2

    call legendre_init(jMid, an)

    guess = pi/2 - 0.5d0*pi/(jMax + 1)
    call newton(jMid, an, guess, x(jMid))
    guess = 3.d0*x(jMid) - pi
    call newton(jMid, an, guess, x(jMid-1))
    do l = jMid-2, 1, -1
      guess = 2*x(l+1) - x(l+2)
      call newton(jMid, an, guess, x(l))
    end do
    do j = 1, jMid
      call legendre_dP(jMid, an, x(j), dpn)
      w(j) = (2.d0*jMax + 1.d0)/(dpn)**2
    end do

    w(jMax:jMid+1:-1) = w(1:jMid)
    x(1:jMid) = cos(x(1:jMid))
    x(jMax:jMid+1:-1) = -x(1:jMid)

  end subroutine f_scgaus
!***********************************************************************
  subroutine legendre_init(jMid,an)
    implicit none
    integer(kind=4), intent(in) :: jMid
    real(kind=8), dimension(0:jMid), intent(out) :: an
    real(kind=8), parameter :: a0sq = 1.5d0
    integer(kind=4) :: k, l, n
!-----------------------------------------------------------------------
    n = jMid*2
    an(jMid) = a0sq
    do k=2, n
      an(jMid) = (1.d0-1.d0/(4.d0*k*k))*an(jMid)
    end do
    an(jMid) = sqrt(an(jMid))
    do k=1, jMid
      l = 2*k
      an(jMid-k) = (l-1.d0)*(2.d0*n-l+2.d0)/(l*(2.d0*n-l+1.d0)) * an(jMid-k+1)
    end do
    an(0) = 0.5d0*an(0)
  end subroutine legendre_init
!***********************************************************************
  subroutine legendre_P(jMid, an, theta, pn)
    implicit none
    integer(kind=4), intent(in) :: jMid
    real(kind=8), dimension(0:jMid), intent(in) :: an
    real(kind=8), intent(in) :: theta
    real(kind=8), intent(out) :: pn
    integer(kind=4) :: k, l
!-----------------------------------------------------------------------
    pn = 0.d0
    do l=0, jMid
      k=l*2                                   ! k = l*2 + 1 if n odd
      pn = pn + an(l) * cos(k*theta)
    end do
  end subroutine legendre_P
!***********************************************************************
  subroutine legendre_dP(jMid, an, theta, dpn)
    implicit none
    integer(kind=4), intent(in) :: jMid
    real(kind=8), dimension(0:jMid), intent(in) :: an
    real(kind=8), intent(in) :: theta
    real(kind=8), intent(out) :: dpn
    integer(kind=4) :: k, l
!-----------------------------------------------------------------------
    dpn = 0.d0
    do l=1, jMid
      k=l*2                                   ! k = l*2 + 1 if n odd
      dpn = dpn - k * an(l) * sin(k*theta)
    end do
  end subroutine legendre_dP
!***********************************************************************
  subroutine newton(jMid, an, x0, x)
    implicit none
    integer(kind=4), intent(in) :: jMid
    real(kind=8), dimension(0:jMid), intent(in) :: an
    real(kind=8), intent(in) :: x0
    real(kind=8), intent(out) :: x
    integer(kind=4), parameter :: newton_max = 500
    real(kind=8), parameter :: xacc_min = 1.d-15
    real(kind=8) :: xacc
    real(kind=8) :: y, dy
    integer(kind=4) :: i
!-----------------------------------------------------------------------
    xacc = xacc_min
    x = x0
    do i = 1, newton_max
      call legendre_P(jMid, an, x, y)
      call legendre_dP(jMid, an, x, dy)
      y = y/dy
      if (abs(y) < xacc) return
      x = x - y
    end do
    write(*,*) '### Error in newton : Too many refinement.'
    stop
  end subroutine newton
!***********************************************************************

END MODULE common_afes
