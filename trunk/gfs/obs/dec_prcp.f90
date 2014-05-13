PROGRAM dec_prcp
!
! NOTE: output goes to fort.9X
!
  USE common
  USE common_gfs
  USE common_obs_gfs

  IMPLICIT NONE

  INTEGER :: iobs
  REAL(r_sngl) :: wk(7)

  integer :: i, j, it
  INTEGER :: iunit, iolen

  integer, parameter :: pplat1 = 22
  integer, parameter :: pplat2 = 73
  integer, parameter :: ppnlat = pplat2-pplat1+1

  real, parameter :: obserr_r = 0.5
  real, parameter :: obserr_min = 0.05

  integer, parameter :: obtype_pr = 21

  real(r_sngl) :: pr(nlon,ppnlat)
!  real(r_sngl) :: tdiff(nlon,ppnlat)


  character (len=80) :: prcpfile

  call set_common_gfs

  prcpfile = 'tmpa.dat'

  wk(1) = id_rain_obs
  wk(4) = -9999.
  wk(7) = REAL(obtype_pr)

  iobs = 0

  INQUIRE(IOLENGTH=iolen) iolen

  open (11, file=trim(prcpfile), form='unformatted', access='direct', recl=nlon*ppnlat*iolen)
  read (11, rec=1) ((pr(i,j),i=1,nlon),j=1,ppnlat)
  close (11)


  DO j = 1, ppnlat
  do i = 1, nlon

    wk(2) = lon(i)
    wk(3) = lat(j+pplat1-1)

    wk(5) = pr(i,j)
    wk(6) = max(wk(5) * obserr_r, obserr_min)

    if (wk(5) >= 0.) then

!      it = ceiling(tdiff(i,j) / 60. - 0.5)
!      if (it < -1) it = -1
!      if (it >  1) it =  1
!      iunit = 90 + it
      iunit = 90


      iobs = iobs + 1
!      write(*,'(I6,1x,7F10.4)') iobs, wk(1:7)
      WRITE(iunit) wk


    end if

  end do
  END DO

  write (*,*) iobs


END PROGRAM dec_prcp
