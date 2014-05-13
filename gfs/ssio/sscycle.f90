!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Cycle forecasts into analyses
!  (must be recomplied for different model resolutions)
!
!  October 2012 - Created by G.-Y. Lien, University of Maryland
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program sscycle

  use sigio_module
  use sfcio_module

  use common, only: r_sngl
  use common_gfs, only: nlon, nlat, nlev, gfs_jcap, gfs_ntrac
  use ssio_tools, only: lusigf, lusfcf, lusiga, lusfca

  implicit none

!-------------------------------------------------------------------------------

  type(sigio_head) :: headsigf, headsiga
  type(sfcio_head) :: headsfcf, headsfca
  type(sigio_data) :: datasigf, datasiga
  type(sfcio_data) :: datasfcf, datasfca

  integer :: j, k, n, iret, iolen

  character(len=7) :: fsigf, fsfcf, fsiga, fsfca
  character(len=120) :: errmsg

!-------------------------------------------------------------------------------
! Open and read files
!-------------------------------------------------------------------------------

  fsigf = 'fort.00'
  fsfcf = 'fort.00'
  fsiga = 'fort.00'
  fsfca = 'fort.00'
  write (fsigf(6:7), '(I2)') lusigf
  write (fsfcf(6:7), '(I2)') lusfcf
  write (fsiga(6:7), '(I2)') lusiga
  write (fsfca(6:7), '(I2)') lusfca

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
! Copy variables from the new analyses to the forecasts
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

end program sscycle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
