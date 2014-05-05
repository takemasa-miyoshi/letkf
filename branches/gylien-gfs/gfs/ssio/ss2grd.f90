!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Convert GFS sig/sfc files to a grid file.
!  (must be recomplied for different model resolutions)
!
!  October 2012 - Adapted from ss2gg.f by G.-Y. Lien, University of Maryland
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program ss2grd

  use sigio_module
  use sfcio_module

  use common, only: r_sngl
  use common_gfs
  use ssio_tools
  use physcons

  implicit none

!-------------------------------------------------------------------------------

  type(sigio_head) :: headsig
  type(sfcio_head) :: headsfc
  type(sigio_data) :: datasig
  type(sfcio_data) :: datasfc

  real(r_sngl) :: f1(nlon*nlat)
  real(r_sngl) :: f2(nlon*nlat)
  real(r_sngl) :: f3(nlon*nlat)
  real(r_sngl) :: g1(nlon*nlat,nlev)
  real(r_sngl) :: g2(nlon*nlat,nlev)
  real(r_sngl) :: pm(nlon*nlat,nlev)

  real(r_sngl) :: v3d(nlon,nlat,nlev,nv3dx)
  real(r_sngl) :: v2d(nlon,nlat,nv2dx)

  integer :: j, k, n, iret, iolen

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
  call jrev(1, f1, v2d(:,:,iv2d_ps))

!--- p: total pressure at model levels (Pa)
  call SPTEZM(0, gfs_jcap, idrt, nlon, nlat, nlev, datasig%t, g1, 1)
  call sigio_modpr(nlon*nlat, nlon*nlat, nlev, headsig%nvcoord, headsig%idvc, headsig%idsl, &
                   headsig%vcoord, iret, f1, g1, pm=pm)
  call jrev(nlev, pm, v3d(:,:,:,iv3d_p))

!--- t: temperature, q: specific humidity (kg/kg)
  call SPTEZM(0, gfs_jcap, idrt, nlon, nlat, nlev, datasig%q(:,:,1), g2, 1)
  g1 = g1 / (1. + con_fvirt * g2)
  call jrev(nlev, g1, v3d(:,:,:,iv3d_t))
  call jrev(nlev, g2, v3d(:,:,:,iv3d_q))

!--- u, v: horizontal winds (m/s)
  call SPTEZMV(0, gfs_jcap, idrt, nlon, nlat, nlev, datasig%d, datasig%z, g1, g2, 1)
  call jrev(nlev, g1, v3d(:,:,:,iv3d_u))
  call jrev(nlev, g2, v3d(:,:,:,iv3d_v))

!--- q3: tracer 3 - cloud condensate (kg/kg)
  call SPTEZM(0, gfs_jcap, idrt, nlon, nlat, nlev, datasig%q(:,:,3), g1, 1)
  call jrev(nlev, g1, v3d(:,:,:,iv3d_qc))

!-------------------------------------------------------------------------------
! Convert variables from the surface file into gridded data, and output
!-------------------------------------------------------------------------------

!--- orog: orography heightmake c (m)
  call jrev(1, datasfc%orog, v2d(:,:,iv2d_orog))

!--- f10m: 10-meter wind speed over lowest model wind speed
  call jrev(1, datasfc%f10m, v2d(:,:,iv2d_f10m))

!--- t2m: 2-meter temperature (K)
  call jrev(1, datasfc%t2m, v2d(:,:,iv2d_t2m))

!--- q2m: 2-meter specific humidity (kg/kg)
  call jrev(1, datasfc%q2m, v2d(:,:,iv2d_q2m))

!--- tprcp: total precipitation rate (mm/s)
  call jrev(1, datasfc%tprcp*3600., v2d(:,:,iv2d_tprcp))

!-------------------------------------------------------------------------------
! Write grid file
!-------------------------------------------------------------------------------

  call write_grd4x(fgrdf, v3d, v2d)

!-------------------------------------------------------------------------------

end program ss2grd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
