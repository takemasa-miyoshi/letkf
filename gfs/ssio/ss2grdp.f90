!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Convert GFS sig/sfc files to a grid file in pressure coordinate.
!  (must be recomplied for different model resolutions)
!
!  October 2012 - Adapted from ss2gg.f by G.-Y. Lien, University of Maryland
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program ss2grdp

  use sigio_module
  use sfcio_module

  use common, only: r_sngl
  use common_gfs
  use common_gfs_pres
  use ssio_tools
  use physcons

  implicit none

!-------------------------------------------------------------------------------

  logical, parameter :: extrap = .true.     ! if do under-surface extrapolation
  real(r_sngl), parameter :: lprate = 0.005 ! standard lapse rate for under-surface extrapolation (K/m)

!-------------------------------------------------------------------------------

  type(sigio_head) :: headsig
  type(sfcio_head) :: headsfc
  type(sigio_data) :: datasig
  type(sfcio_data) :: datasfc

  real(r_sngl) :: f1(nlon*nlat)
  real(r_sngl) :: f2(nlon*nlat)
  real(r_sngl) :: g1(nlon*nlat,nlev)
  real(r_sngl) :: g2(nlon*nlat,nlev)
  real(r_sngl) :: tt(nlon*nlat,nlev)
  real(r_sngl) :: pm(nlon*nlat,nlev)
  real(r_sngl) :: rh(nlev)
  real(r_sngl) :: gph(nlon*nlat,nlev+1)
  real(r_sngl) :: pstag(nlon*nlat,nlev+1)
  real(r_sngl) :: tsext(nlon*nlat)
  real(r_sngl) :: tvsext(nlon*nlat)

  real(r_sngl) :: p_u(nlon*nlat,nlevp)
  real(r_sngl) :: p_v(nlon*nlat,nlevp)
  real(r_sngl) :: p_t(nlon*nlat,nlevp)
  real(r_sngl) :: p_q(nlon*nlat,nlevp)
  real(r_sngl) :: p_rh(nlon*nlat,nlevp)
  real(r_sngl) :: p_q2(nlon*nlat,nlevp)
  real(r_sngl) :: p_q3(nlon*nlat,nlevp)
  real(r_sngl) :: p_gph(nlon*nlat,nlevp)
  real(r_sngl) :: slp(nlon*nlat)
  real(r_sngl) :: u10m(nlon*nlat)
  real(r_sngl) :: v10m(nlon*nlat)
  
  real(r_sngl) :: v3d(nlon,nlat,nlevp,nv3dp)
  real(r_sngl) :: v2d(nlon,nlat,nv2dp)

  integer :: i, j, ij, k, n, iret, iolen

  character(len=7) :: fsigf, fsfcf, fgrdf
  character(len=120) :: errmsg

!-------------------------------------------------------------------------------
! Open and read files
!-------------------------------------------------------------------------------

  fsigf = 'fort.00'
  fsfcf = 'fort.00'
  fgrdf = 'fort.00'
  write (fsigf(6:7), '(I2)') lusigf
  write (fsfcf(6:7), '(I2)') lusfcf
  write (fgrdf(6:7), '(I2)') lugrdf

  call sigio_srohdc(lusigf, fsigf, headsig, datasig, iret)
  if (iret /= 0) then
    write (errmsg, *) '[Error] Open and read sigma file: ', fsigf
    write (0, '(A)') trim(errmsg)
    stop
  endif

  call sfcio_srohdc(lusfcf, fsfcf, headsfc, datasfc, iret)
  if (iret /= 0) then
    write (errmsg, *) '[Error] Open and read surface file: ', fsfcf
    write (0, '(A)') trim(errmsg)
    stop
  endif
  
  if (headsig%jcap /= gfs_jcap .OR. &
      headsig%lonb /= nlon .OR. &
      headsfc%lonb /= nlon .OR. &
      headsig%latb /= nlat .OR. &
      headsfc%latb /= nlat .OR. &
      headsig%levs /= nlev .OR. &
      headsig%ntrac /= gfs_ntrac) then
    stop '[Error] Resolutions of input files mismatch the pre-compiled settings.'
  end if

!-------------------------------------------------------------------------------
! Convert variables from the sigma file into gridded data, and output
!-------------------------------------------------------------------------------

!--- ps: surface pressure (Pa)
  call SPTEZ(0, gfs_jcap, idrt, nlon, nlat, datasig%ps, f1, 1)
  f1 = exp(f1) * 1.e3
  call jrev(1, f1/100., v2d(:,:,iv2dp_ps))

!--- orog: orography height (m)
  call SPTEZ(0, gfs_jcap, idrt, nlon, nlat, datasig%hs, f2, 1)

!--- p: total pressure at model levels (Pa)
  call SPTEZM(0, gfs_jcap, idrt, nlon, nlat, nlev, datasig%t, g1, 1)
  call sigio_modpr(nlon*nlat, nlon*nlat, nlev, headsig%nvcoord, headsig%idvc, headsig%idsl, &
                   headsig%vcoord, iret, f1, g1, pd=g2, pm=pm)

  do ij = 1, nlon*nlat
    call getgph(nlev, g2(ij,:), g1(ij,:), f1(ij), f2(ij), gph(ij,:), pstag(ij,:))
    call itpl_p(nlev+1, pstag(ij,:), gph(ij,:), nlevp, levp, missing, p_gph(ij,:))
  end do

!--- t: temperature, q: specific humidity (kg/kg), rh: relative humidity (%)
  call SPTEZM(0, gfs_jcap, idrt, nlon, nlat, nlev, datasig%q(:,:,1), g2, 1)
  tt = g1 / (1. + con_fvirt * g2)

  do ij = 1, nlon*nlat
    call getrh(nlev, pm(ij,:), g2(ij,:), tt(ij,:), rh)
    call itpl_p(nlev, pm(ij,:), tt(ij,:), nlevp, levp, missing, p_t(ij,:))
    if (extrap) then
      call itpl_p(nlev, pm(ij,:), g2(ij,:), nlevp, levp, g2(ij,1), p_q(ij,:))
      call itpl_p(nlev, pm(ij,:), rh, nlevp, levp, rh(1), p_rh(ij,:))
    else
      call itpl_p(nlev, pm(ij,:), g2(ij,:), nlevp, levp, missing, p_q(ij,:))
      call itpl_p(nlev, pm(ij,:), rh, nlevp, levp, missing, p_rh(ij,:))
    end if
    call sfc_extrap(nlev, pm(ij,:), tt(ij,:), f1(ij), g2(ij,1), lprate, tsext(ij), tvsext(ij))
  end do

!--- u, v: horizontal winds (m/s)
  call SPTEZMV(0, gfs_jcap, idrt, nlon, nlat, nlev, datasig%d, datasig%z, g1, g2, 1)

  do ij = 1, nlon*nlat
    if (extrap) then
      call itpl_p(nlev, pm(ij,:), g1(ij,:), nlevp, levp, g1(ij,1), p_u(ij,:))
      call itpl_p(nlev, pm(ij,:), g2(ij,:), nlevp, levp, g2(ij,1), p_v(ij,:))
    else
      call itpl_p(nlev, pm(ij,:), g1(ij,:), nlevp, levp, missing, p_u(ij,:))
      call itpl_p(nlev, pm(ij,:), g2(ij,:), nlevp, levp, missing, p_v(ij,:))
    end if
  end do

  do j = 1, nlat
    do i = 1, nlon
      ij = (j-1)*nlon + i
      u10m(ij) = datasfc%f10m(i,j) * g1(ij,1)
      v10m(ij) = datasfc%f10m(i,j) * g2(ij,1)
    end do
  end do

!--- q2: tracer 2 - ozone (kg/kg)
  call SPTEZM(0, gfs_jcap, idrt, nlon, nlat, nlev, datasig%q(:,:,2), g1, 1)

  do ij = 1, nlon*nlat
    if (extrap) then
      call itpl_p(nlev, pm(ij,:), g1(ij,:), nlevp, levp, g1(ij,1), p_q2(ij,:))
    else
      call itpl_p(nlev, pm(ij,:), g1(ij,:), nlevp, levp, missing, p_q2(ij,:))
    end if
  end do

!--- q3: tracer 3 - cloud condensate (kg/kg)
  call SPTEZM(0, gfs_jcap, idrt, nlon, nlat, nlev, datasig%q(:,:,3), g1, 1)

  do ij = 1, nlon*nlat
    if (extrap) then
      call itpl_p(nlev, pm(ij,:), g1(ij,:), nlevp, levp, g1(ij,1), p_q3(ij,:))
    else
      call itpl_p(nlev, pm(ij,:), g1(ij,:), nlevp, levp, missing, p_q3(ij,:))
    end if
  end do

!-------------------------------------------------------------------------------
! Compute SLP, extrapolated geopotential height and temperature
!-------------------------------------------------------------------------------

  do ij = 1, nlon*nlat
    slp(ij) = f1(ij) * (1 + lprate * f2(ij) / tvsext(ij)) ** (1./con_rog/lprate)
  end do
  call jrev(1, slp/100., v2d(:,:,iv2dp_slp))

  if (extrap) then    
    do ij = 1, nlon*nlat
      do k = 1, nlevp
        if (p_gph(ij,k) == missing) then
          p_gph(ij,k) = f2(ij) + tvsext(ij)/lprate * (1. - (levp(k)/f1(ij)) ** (con_rog*lprate))
        end if
        if (p_t(ij,k) == missing) then
          p_t(ij,k) = tsext(ij) + lprate * (f2(ij) - p_gph(ij,k))
        end if
      end do
    end do
  end if

!-------------------------------------------------------------------------------

  call jrev(nlevp, p_u  , v3d(:,:,:,iv3dp_u    ))
  call jrev(nlevp, p_v  , v3d(:,:,:,iv3dp_v    ))
  call jrev(nlevp, p_t  , v3d(:,:,:,iv3dp_t    ))
  call jrev(nlevp, p_q  , v3d(:,:,:,iv3dp_q    ))
  call jrev(nlevp, p_rh , v3d(:,:,:,iv3dp_rh   ))
  call jrev(nlevp, p_q3 , v3d(:,:,:,iv3dp_qc   ))
  call jrev(nlevp, p_gph, v3d(:,:,:,iv3dp_gph  ))
  call jrev(nlevp, p_q2 , v3d(:,:,:,iv3dp_ozone))

!-------------------------------------------------------------------------------
! Convert variables from the surface file into gridded data, and output
!-------------------------------------------------------------------------------

!--- orog: orography height (m)
  call jrev(1, datasfc%orog, v2d(:,:,iv2dp_orog))

!--- slmsk: sea-land-ice mask (0-sea, 1-land, 2-ice)
  call jrev(1, datasfc%slmsk, v2d(:,:,iv2dp_slmsk))

!--- sst: SST (K)
  call jrev(1, datasfc%tsea, v2d(:,:,iv2dp_tsea))

!--- u10m, v10m: 10-meter wind speeds
  call jrev(1, u10m, v2d(:,:,iv2dp_u10m))
  call jrev(1, v10m, v2d(:,:,iv2dp_v10m))

!--- t2m: 2-meter temperature (K)
  call jrev(1, datasfc%t2m, v2d(:,:,iv2dp_t2m))

!--- q2m: 2-meter specific humidity (kg/kg)
  call jrev(1, datasfc%q2m, v2d(:,:,iv2dp_q2m))

!--- tprcp: total precipitation rate (mm/s)
  call jrev(1, datasfc%tprcp*3600., v2d(:,:,iv2dp_tprcp))
  
!-------------------------------------------------------------------------------
! Write grid file
!-------------------------------------------------------------------------------

  call write_grd4p(fgrdf, v3d, v2d)

!-------------------------------------------------------------------------------

end program ss2grdp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
