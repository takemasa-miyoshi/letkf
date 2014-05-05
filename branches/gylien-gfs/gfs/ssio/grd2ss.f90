!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Convert a grid file to GFS sig/sfc formats, cycle forecasts into analyses
!  (must be recomplied for different model resolutions)
!
!  fsigf/fsfcf: (GFS sig/sfc) [Input]        Input all other forecast fields (from GFS 6-hour forecast)
!  fsiga/fsfca: (GFS sig/sfc) [Input/Output] Input ozone and tsea (from other reference analysis);
!                                            Input head info;
!                                            Output our analysis
!  fgrda:       (LETKF grd)   [Input]        Input state variables (u,v,t,q,ps) (from LETKF analysis)
!
!  October 2012 - Created by G.-Y. Lien, University of Maryland
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program grd2ss

  use sigio_module
  use sfcio_module

  use common, only: r_sngl
  use common_gfs
  use ssio_tools
  use physcons

  implicit none

!-------------------------------------------------------------------------------

  type(sigio_head) :: headsigf, headsiga
  type(sfcio_head) :: headsfcf, headsfca
  type(sigio_data) :: datasigf, datasiga
  type(sfcio_data) :: datasfcf, datasfca

  real(r_sngl) :: f1(nlon*nlat)
  real(r_sngl) :: g1(nlon*nlat,nlev)
  real(r_sngl) :: g2(nlon*nlat,nlev)

  real(r_sngl) :: v3d(nlon,nlat,nlev,nv3d)
  real(r_sngl) :: v2d(nlon,nlat,nv2d)

  integer :: j, k, n, iret, iolen

  character(len=7) :: fsigf, fsfcf, fsiga, fsfca, fgrda
  character(len=120) :: errmsg

!-------------------------------------------------------------------------------
! Open and read files
!-------------------------------------------------------------------------------

  fsigf = 'fort.00'
  fsfcf = 'fort.00'
  fsiga = 'fort.00'
  fsfca = 'fort.00'
  fgrda = 'fort.00'
  write (fsigf(6:7), '(I2)') lusigf
  write (fsfcf(6:7), '(I2)') lusfcf
  write (fsiga(6:7), '(I2)') lusiga
  write (fsfca(6:7), '(I2)') lusfca
  write (fgrda(6:7), '(I2)') lugrda

  call sigio_srohdc(lusigf, fsigf, headsigf, datasigf, iret)
  if (iret /= 0) then
    write (errmsg, *) '[Error] Open and read sigma file: ', fsigf
    write (0, '(A)') trim(errmsg)
    stop
  endif

  call sigio_srohdc(lusiga, fsiga, headsiga, datasiga, iret)
  if (iret /= 0) then
    write (errmsg, *) '[Error] Open and read sigma file: ', fsiga
    write (0, '(A)') trim(errmsg)
    stop
  endif

  call sfcio_srohdc(lusfcf, fsfcf, headsfcf, datasfcf, iret)
  if (iret /= 0) then
    write (errmsg, *) '[Error] Open and read surface file: ', fsfcf
    write (0, '(A)') trim(errmsg)
    stop
  endif

  call sfcio_srohdc(lusfca, fsfca, headsfca, datasfca, iret)
  if (iret /= 0) then
    write (errmsg, *) '[Error] Open and read surface file: ', fsfca
    write (0, '(A)') trim(errmsg)
    stop
  endif

  if (headsigf%jcap /= gfs_jcap .OR. &
      headsigf%lonb /= nlon .OR. &
      headsfcf%lonb /= nlon .OR. &
      headsigf%latb /= nlat .OR. &
      headsfcf%latb /= nlat .OR. &
      headsigf%levs /= nlev .OR. &
      headsigf%ntrac /= gfs_ntrac .OR. &
      headsiga%jcap /= gfs_jcap .OR. &
      headsiga%lonb /= nlon .OR. &
      headsfca%lonb /= nlon .OR. &
      headsiga%latb /= nlat .OR. &
      headsfca%latb /= nlat .OR. &
      headsiga%levs /= nlev .OR. &
      headsiga%ntrac /= gfs_ntrac) then
    stop '[Error] Resolutions of input files mismatch the pre-compiled settings.'
  end if

!-------------------------------------------------------------------------------
! Read analysis grid file
!-------------------------------------------------------------------------------

  call read_grd4(fgrda, v3d, v2d, 0)

!-------------------------------------------------------------------------------
! Convert variables from gridded data to sig/sfc data, update the forecast
!-------------------------------------------------------------------------------

!--- ps: surface pressure (Pa)
  call jrev(1, v2d(:,:,iv2d_ps), f1)
  f1 = log(f1 / 1.e3)
  call SPTEZ(0, gfs_jcap, idrt, nlon, nlat, datasigf%ps, f1, -1)

!--- tv: virtual temperature, q: specific humidity (kg/kg)
  call jrev(nlev, v3d(:,:,:,iv3d_t), g1)
  call jrev(nlev, v3d(:,:,:,iv3d_q), g2)
  g1 = g1 * (1. + con_fvirt * g2)
  call SPTEZM(0, gfs_jcap, idrt, nlon, nlat, nlev, datasigf%t, g1, -1)
  call SPTEZM(0, gfs_jcap, idrt, nlon, nlat, nlev, datasigf%q(:,:,1), g2, -1)

!--- u, v: horizontal winds (m/s)
  call jrev(nlev, v3d(:,:,:,iv3d_u), g1)
  call jrev(nlev, v3d(:,:,:,iv3d_v), g2)
  call SPTEZMV(0, gfs_jcap, idrt, nlon, nlat, nlev, datasigf%d, datasigf%z, g1, g2, -1)

!--- q3: tracer 3 - cloud condensate (kg/kg)
  call jrev(nlev, v3d(:,:,:,iv3d_qc), g1)
  call SPTEZM(0, gfs_jcap, idrt, nlon, nlat, nlev, datasigf%q(:,:,3), g1, -1)

!-------------------------------------------------------------------------------
! Copy variables from the reference analyses to the forecast
!-------------------------------------------------------------------------------

!--- Sigma data: ozone
  datasigf%q(:,:,2) = datasiga%q(:,:,2)

!--- Surface data: tsea
  datasfcf%tsea = datasfca%tsea

!-------------------------------------------------------------------------------
! Write sig/sfc files
!-------------------------------------------------------------------------------

  call sigio_swohdc(lusiga, fsiga, headsiga, datasigf, iret) ! write forecast data (datasigf) into analysis
  if (iret /= 0) then
    write (errmsg, *) '[Error] Overwrite sigma file: ', fsiga
    write (0, '(A)') trim(errmsg)
    stop
  endif

  call sfcio_swohdc(lusfca, fsfca, headsfca, datasfcf, iret) ! write forecast data (datasfcf) into analysis
  if (iret /= 0) then
    write (errmsg, *) '[Error] Overwrite surface file: ', fsfca
    write (0, '(A)') trim(errmsg)
    stop
  endif

!-------------------------------------------------------------------------------

end program grd2ss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
