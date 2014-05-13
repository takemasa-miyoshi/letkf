!===============================================================================
module ssio_tools
!-------------------------------------------------------------------------------

  use common, only: r_sngl
  use common_gfs
  use common_gfs_pres
  use physcons
  implicit none

  integer, parameter :: lusigf = 11, lusfcf = 12
  integer, parameter :: lusiga = 21, lusfca = 22
  integer, parameter :: lugrdf = 31
  integer, parameter :: lugrda = 41

  real(r_sngl), parameter :: missing = -9.99e33 ! missing value in grads files

contains

!===============================================================================
subroutine jrev (km, vdati, vdato)
!-------------------------------------------------------------------------------

  implicit none

  integer, intent(in)       :: km
  real(r_sngl), intent(in)  :: vdati(nlon,nlat,km) ! implicit reshape
  real(r_sngl), intent(out) :: vdato(nlon,nlat,km)

  vdato(:,:,:) = vdati(:,nlat:1:-1,:)

!-------------------------------------------------------------------------------
end subroutine jrev
!===============================================================================

!===============================================================================
subroutine sfc_extrap (km, p, t, ps, qs, lprate, tsext, tvsext)
!-------------------------------------------------------------------------------

  implicit none
  integer, intent(in)       :: km
  real(r_sngl), intent(in)  :: p(km), t(km)
  real(r_sngl), intent(in)  :: ps, qs
  real(r_sngl), intent(in)  :: lprate
  real(r_sngl), intent(out) :: tsext, tvsext
  real(r_sngl) :: ptarg, dp, dpmin
  integer :: k, ktarg

!-------------------------------------------------------------------------------

  ptarg = ps - 15000.
  dpmin = 1.e6
  do k = km, 1, -1
    dp = abs(p(k) - ptarg)
    if (dp > dpmin) exit
    dpmin = dp
  end do
  ktarg = k + 1

  tsext = t(ktarg) * (ps / p(ktarg)) ** (con_rog * lprate)
  tvsext = tsext * (1. + con_fvirt * qs)

!-------------------------------------------------------------------------------
end subroutine sfc_extrap
!===============================================================================

!===============================================================================
subroutine itpl_p (km, p, var, kmp, levp, undf, varp)
!-------------------------------------------------------------------------------

  implicit none
  integer, intent(in)       :: km
  real(r_sngl), intent(in)  :: p(km), var(km)
  integer, intent(in)       :: kmp
  real(r_sngl), intent(in)  :: levp(kmp)
  real(r_sngl), intent(in)  :: undf
  real(r_sngl), intent(out) :: varp(kmp)
  real(r_sngl) :: logp(km)
  real(r_sngl) :: rk, ak
  integer :: k, kp

!-------------------------------------------------------------------------------

  logp = log(p)
  do kp = 1, kmp
    rk = log(levp(kp))
    if (logp(1) < rk) then
      varp(kp) = undf
      cycle
    end if
    do k = 2, km
      if (logp(k) < rk) exit
    end do
    ak = (rk - logp(k-1)) / (logp(k) - logp(k-1))
    varp(kp) = var(k-1) * (1-ak) + var(k) * ak
  end do

!-------------------------------------------------------------------------------
end subroutine itpl_p
!===============================================================================

!===============================================================================
subroutine getgph (km, dp, tv, ps, orog, gph, pstag)
!-------------------------------------------------------------------------------

  implicit none

  integer, intent(in)       :: km
  real(r_sngl), intent(in)  :: dp(km), tv(km)
  real(r_sngl), intent(in)  :: ps, orog
  real(r_sngl), intent(out) :: gph(km+1)
  real(r_sngl), intent(out) :: pstag(km+1)
  integer :: k

!-------------------------------------------------------------------------------

  gph(1) = orog
  pstag(1) = ps
  do k = 1, km-1
    pstag(k+1) = pstag(k) - dp(k)
    gph(k+1) = gph(k) + con_rog * tv(k) * log(pstag(k) / pstag(k+1))
  enddo
  pstag(km+1) = 0.
  gph(km+1) = 9.99e33

!-------------------------------------------------------------------------------
end subroutine getgph
!===============================================================================

!===============================================================================
subroutine getrh (km, p, sh, t, rh)
!-------------------------------------------------------------------------------

  implicit none

  integer, intent(in)  :: km
  real(r_sngl), intent(in)  :: p(km), sh(km), t(km)
  real(r_sngl), intent(out) :: rh(km)
  real(r_sngl) :: tr, es, shs
  integer :: k

!-------------------------------------------------------------------------------

  do k = 1, km
    tr = con_ttp / t(k)
    es = con_psat * (tr**con_xpona) * exp(con_xponb*(1.-tr))
    es = min(es, p(k))
    shs = con_eps * es / (p(k) + con_epsm1 * es)
    rh(k) = 1.e2 * min(max(sh(k)/shs, 0.), 1.)
!    rh(k) = 1.e2 * sh(k) / shs
  enddo

!-------------------------------------------------------------------------------
end subroutine getrh
!===============================================================================

!===============================================================================
subroutine write_ctl (fid, dset, options, tstart, tint, tnum)
!-------------------------------------------------------------------------------

  implicit none

  integer, intent(in)          :: fid
  character(len=*), intent(in) :: dset
  character(len=*), intent(in) :: options
  character(len=*), intent(in) :: tstart
  character(len=*), intent(in) :: tint
  integer, intent(in)          :: tnum

  integer :: j

!-------------------------------------------------------------------------------

  write (fid, '(2A)') 'dset ^', dset
  write (fid, '(2A)') 'options ', options
  write (fid, '(A,ES12.5)') 'undef ', missing
  write (fid, '(A,I6,A,2F12.6)') 'xdef', nlon, ' linear', 0., 360./real(nlon)
  write (fid, '(A,I6,A)') 'ydef', nlat, ' levels'
  do j = 1, nlat
    write (fid, '(F12.6)') lat(j)
  end do
  write (fid, '(A,I6,A)') 'zdef', nlev, ' linear 1 1'
  write (fid, '(A,I9,4A)') 'tdef', tnum, ' linear ', tstart, ' ', tint
  write (fid, '(A,I6)') 'vars', nv3d+nv2d
  write (fid, '(A,2I4,A)') 'u    ', nlev, 99, '  zonal wind (m/s)'
  write (fid, '(A,2I4,A)') 'v    ', nlev, 99, '  meridional wind (m/s)'
  write (fid, '(A,2I4,A)') 't    ', nlev, 99, '  temperature (K)'
  write (fid, '(A,2I4,A)') 'q    ', nlev, 99, '  specific humidity (kg/kg)'
  write (fid, '(A,2I4,A)') 'qc   ', nlev, 99, '  cloud condensate (kg/kg)'
  write (fid, '(A,2I4,A)') 'ps   ', 0,    99, '  surface pressure (Pa)'
  write (fid, '(A)') 'endvars'

!-------------------------------------------------------------------------------
end subroutine write_ctl
!===============================================================================

!===============================================================================
subroutine write_ctlx (fid, dset, options, tstart, tint, tnum)
!-------------------------------------------------------------------------------

  implicit none

  integer, intent(in)          :: fid
  character(len=*), intent(in) :: dset
  character(len=*), intent(in) :: options
  character(len=*), intent(in) :: tstart
  character(len=*), intent(in) :: tint
  integer, intent(in)          :: tnum

  integer :: j

!-------------------------------------------------------------------------------

  write (fid, '(2A)') 'dset ^', dset
  write (fid, '(2A)') 'options ', options
  write (fid, '(A,ES12.5)') 'undef ', missing
  write (fid, '(A,I6,A,2F12.6)') 'xdef', nlon, ' linear', 0., 360./real(nlon)
  write (fid, '(A,I6,A)') 'ydef', nlat, ' levels'
  do j = 1, nlat
    write (fid, '(F12.6)') lat(j)
  end do
  write (fid, '(A,I6,A)') 'zdef', nlev, ' linear 1 1'
  write (fid, '(A,I9,4A)') 'tdef', tnum, ' linear ', tstart, ' ', tint
  write (fid, '(A,I6)') 'vars', nv3dx+nv2dx
  write (fid, '(A,2I4,A)') 'u    ', nlev, 99, '  zonal wind (m/s)'
  write (fid, '(A,2I4,A)') 'v    ', nlev, 99, '  meridional wind (m/s)'
  write (fid, '(A,2I4,A)') 't    ', nlev, 99, '  temperature (K)'
  write (fid, '(A,2I4,A)') 'q    ', nlev, 99, '  specific humidity (kg/kg)'
  write (fid, '(A,2I4,A)') 'qc   ', nlev, 99, '  cloud condensate (kg/kg)'
  write (fid, '(A,2I4,A)') 'p    ', nlev, 99, '  pressure (Pa)'
  write (fid, '(A,2I4,A)') 'ps   ', 0,    99, '  surface pressure (Pa)'
  write (fid, '(A,2I4,A)') 'orog ', 0,    99, '  orography height (m)'
  write (fid, '(A,2I4,A)') 'f10m ', 0,    99, '  10-meter wind speed over the lowest level wind speed'
  write (fid, '(A,2I4,A)') 't2m  ', 0,    99, '  2-meter temperature (K)'
  write (fid, '(A,2I4,A)') 'q2m  ', 0,    99, '  2-meter specific humidity (kg/kg)'
  write (fid, '(A,2I4,A)') 'tprcp', 0,    99, '  total precipitation rate (mm/h)'
  write (fid, '(A)') 'endvars'

!-------------------------------------------------------------------------------
end subroutine write_ctlx
!===============================================================================

!===============================================================================
subroutine write_ctlp (fid, dset, options, tstart, tint, tnum)
!-------------------------------------------------------------------------------

  implicit none

  integer, intent(in)          :: fid
  character(len=*), intent(in) :: dset
  character(len=*), intent(in) :: options
  character(len=*), intent(in) :: tstart
  character(len=*), intent(in) :: tint
  integer, intent(in)          :: tnum

  integer :: j, k

!-------------------------------------------------------------------------------

  write (fid, '(2A)') 'dset ^', dset
  write (fid, '(2A)') 'options ', options
  write (fid, '(A,ES12.5)') 'undef ', missing
  write (fid, '(A,I6,A,2F12.6)') 'xdef', nlon, ' linear', 0., 360./real(nlon)
  write (fid, '(A,I6,A)') 'ydef', nlat, ' levels'
  do j = 1, nlat
    write (fid, '(F12.6)') lat(j)
  end do
  write (fid, '(A,I6,A)') 'zdef', nlevp, ' levels'
  do k = 1, nlevp
    write (fid, '(F12.6)') levp(k)/100.
  end do
  write (fid, '(A,I9,4A)') 'tdef', tnum, ' linear ', tstart, ' ', tint
  write (fid, '(A,I6)') 'vars', nv3dp+nv2dp
  write (fid, '(A,2I4,A)') 'u    ', nlevp, 99, '  zonal wind (m/s)'
  write (fid, '(A,2I4,A)') 'v    ', nlevp, 99, '  meridional wind (m/s)'
  write (fid, '(A,2I4,A)') 't    ', nlevp, 99, '  temperature (K)'
  write (fid, '(A,2I4,A)') 'q    ', nlevp, 99, '  specific humidity (kg/kg)'
  write (fid, '(A,2I4,A)') 'rh   ', nlevp, 99, '  relative humidity (%)'
  write (fid, '(A,2I4,A)') 'qc   ', nlevp, 99, '  cloud condensate (kg/kg)'
  write (fid, '(A,2I4,A)') 'gph  ', nlevp, 99, '  geopotential height (m)'
  write (fid, '(A,2I4,A)') 'ozone', nlevp, 99, '  ozone (kg/kg)'
  write (fid, '(A,2I4,A)') 'ps   ', 0,     99, '  surface pressure (Pa)'
  write (fid, '(A,2I4,A)') 'slp  ', 0,     99, '  sea level pressure (hPa)'
  write (fid, '(A,2I4,A)') 'orog ', 0,     99, '  orography height (m)'
  write (fid, '(A,2I4,A)') 'slmsk', 0,     99, '  sea-land-ice mask (0-sea, 1-land, 2-ice)'
  write (fid, '(A,2I4,A)') 'tsea ', 0,     99, '  sea surface temperature (K)'
  write (fid, '(A,2I4,A)') 'u10m ', 0,     99, '  10-meter u-wind speed (m/s)'
  write (fid, '(A,2I4,A)') 'v10m ', 0,     99, '  10-meter v-wind speed (m/s)'
  write (fid, '(A,2I4,A)') 't2m  ', 0,     99, '  2-meter temperature (K)'
  write (fid, '(A,2I4,A)') 'q2m  ', 0,     99, '  2-meter specific humidity (kg/kg)'
  write (fid, '(A,2I4,A)') 'tprcp', 0,     99, '  total precipitation rate (mm/h)'
  write (fid, '(A)') 'endvars'

!-------------------------------------------------------------------------------
end subroutine write_ctlp
!===============================================================================

end module ssio_tools
!===============================================================================
