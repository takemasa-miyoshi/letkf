!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  module common_precip
!
!  module for precipitation assimilation
!
!  created  Jan. 2012, Guo-Yuan Lien, UMD
!  adopted to GFS-LETKF and modified, May 2013, Guo-Yuan Lien, UMD
!  modified, Spetember 2013, Guo-Yuan Lien, UMD
!
!  function dinvnorm(p) modified from Ren-Raw Chen, 
!    Rutgers University in New Brunswick, New Jersey
!    http://home.online.no/~pjacklam/notes/invnorm/
!
!-------------------------------------------------------------------------------
!
!  subroutine read_ppcdf     (cdffile_m, cdffile_o, ppcdf_m, ppcdf_o, ppzero_m, ppzero_o)
!  subroutine read_ppmask    (maskfile, ppmask)
!  function   pptrans_normal (pp, ppcdf, ppzero)
!  function   pptrans_log    (pp)
!  subroutine pptrans_normal_mdzero_def (pp_ens, ppcdf, ppzero,           zero_mem, ym, sigma)
!  function   pptrans_normal_mdzero     (pp,     ppcdf, ppzero, ppzero_m, zero_mem, ym, sigma)
!  function   compact_tail   (pos_cdf)
!  function   dinvnorm       (p)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module common_precip

  use common
  use common_gfs, only: nlon, nlat
  use common_letkf, only: nbv
  implicit none

!-------------------------------------------------------------------------------

  integer, parameter :: ncdf = 200   ! # of cdf bins

  real(r_size), parameter :: ppzero_thres = 0.06d0 ! threshold of no precipitation
  real(r_size), parameter :: mask_thres = 0.35d0   ! threshold of assimilation area wrt. the mask file

  integer, parameter :: opt_pptrans = 3  ! 0: no transformation
                                         ! 1: log transformation
                                         ! 2: Gaussian transformation with median zero rain
                                         ! 3: Gaussian transformation with modified median zero rain
  real(r_size), parameter :: log_trans_tiny = 0.6d0
  real(r_size), parameter :: gausstail_thres = 0.001d0

  integer, parameter :: opt_ppobserr = 2 ! 0: original obserr form  obs data file
                                         ! 1: transformed obserr from obs data file
                                         ! 2: constant obserr
  real(r_size), parameter :: const_ppobserr = 0.5d0
  real(r_size), parameter :: min_ppobserr = 0.1d0

!-------------------------------------------------------------------------------

  integer, parameter :: pp_bg_nlev = 2
  integer, parameter :: pp_bg_levs(pp_bg_nlev-1) = &
                        (/24/)
  integer, parameter :: pp_ob_nlev = 2
  real(r_size), parameter :: pp_ob_levs(pp_ob_nlev-1) = &
                        (/ppzero_thres/)
  logical, parameter :: pp_criterion(pp_bg_nlev,pp_ob_nlev) = reshape((/ &
!           bg1   , bg2
           .false.,.true., &  ! ob1
           .false.,.true.  &  ! ob2
           /), (/pp_bg_nlev,pp_ob_nlev/))

!  integer, parameter :: pp_bg_nlev = 4
!  integer, parameter :: pp_bg_levs(pp_bg_nlev-1) = &
!                        (/20,24,28/)
!  integer, parameter :: pp_ob_nlev = 3
!  real(r_size), parameter :: pp_ob_levs(pp_ob_nlev-1) = &
!                        (/ppzero_thres,1.d0/)
!  logical, parameter :: pp_criterion(pp_bg_nlev,pp_ob_nlev) = reshape((/ &
!!           bg1   , bg2,  , bg3   , bg4
!           .false.,.false.,.false.,.true., &  ! ob1
!           .false.,.false.,.true. ,.true., &  ! ob2
!           .false.,.true. ,.true. ,.true.  &  ! ob3
!           /), (/pp_bg_nlev,pp_ob_nlev/))

!-------------------------------------------------------------------------------

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_ppcdf (cdffile_m, cdffile_o, ppcdf_m, ppcdf_o, ppzero_m, ppzero_o)

  implicit none

  character(len=*), intent(in) :: cdffile_m
  character(len=*), intent(in) :: cdffile_o
  real(r_size), intent(out) :: ppcdf_m(nlon,nlat,0:ncdf)
  real(r_size), intent(out) :: ppcdf_o(nlon,nlat,0:ncdf)
  real(r_size), intent(out) :: ppzero_m(nlon,nlat)
  real(r_size), intent(out) :: ppzero_o(nlon,nlat)

  real(r_sngl) :: ppcdf_ms(nlon,nlat,0:ncdf)
  real(r_sngl) :: ppzero_ms(nlon,nlat)
  real(r_sngl) :: ppcdf_os(nlon,nlat,0:ncdf)
  real(r_sngl) :: ppzero_os(nlon,nlat)
  integer :: i, j, b, iolen
  logical :: ex

  inquire (iolength=iolen) iolen

  inquire (file=trim(cdffile_m), exist=ex)
  if (ex) then
    open (90, file=trim(cdffile_m), status='old', form='unformatted', &
              access='direct', recl=iolen*nlon*nlat)
    do b = 0, ncdf
      read (90, rec=(b+1)) ((ppcdf_ms(i,j,b), i=1,nlon), j=1,nlat)
    end do
    read (90, rec=(2*(ncdf+1)+2)) ((ppzero_ms(i,j), i=1,nlon), j=1,nlat)
    close (90)
    ppcdf_m = real(ppcdf_ms, r_size)
    ppzero_m = real(ppzero_ms, r_size)
  else
    write (6,'(3A)') "CDF file ", cdffile_m, " does not exist -- skipped"
  end if

  inquire (file=trim(cdffile_o), exist=ex)
  if (ex) then
    open (91, file=trim(cdffile_o), status='old', form='unformatted', &
              access='direct', recl=iolen*nlon*nlat)
    do b = 0, ncdf
      read (91, rec=(b+1)) ((ppcdf_os(i,j,b), i=1,nlon), j=1,nlat)
    end do
    read (91, rec=(2*(ncdf+1)+2)) ((ppzero_os(i,j), i=1,nlon), j=1,nlat)
    close (91)
    ppcdf_o = real(ppcdf_os, r_size)
    ppzero_o = real(ppzero_os, r_size)
  else
    write (6,'(3A)') "CDF file ", cdffile_o, " does not exist -- skipped"
  end if

end subroutine read_ppcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_ppmask (maskfile, ppmask)

  implicit none

  character(len=*), intent(in) :: maskfile
  real(r_size), intent(out) :: ppmask(nlon,nlat)

  real(r_sngl) :: ppmask_s(nlon,nlat)
  integer :: i, j, iolen
  logical :: ex

  inquire (iolength=iolen) iolen

  inquire (file=trim(maskfile), exist=ex)
  if (ex) then
    open (92, file=trim(maskfile), status='old', form='unformatted', &
              access='direct', recl=iolen*nlon*nlat)
    read (92, rec=1) ((ppmask_s(i,j), i=1,nlon), j=1,nlat)
    close (92)
    ppmask = real(ppmask_s, r_size)
  else
    write (6,'(3A)') "Mask file ", maskfile, " does not exist -- skipped"
    ppmask = 1.0e10  ! All data are used.
  end if

end subroutine read_ppmask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function pptrans_normal (pp, ppcdf, ppzero)

  implicit none

  real(r_size) :: pptrans_normal

  real(r_size), intent(in) :: pp
  real(r_size), intent(in) :: ppcdf(0:ncdf)
  real(r_size), intent(in) :: ppzero

  real(r_size) :: pos_cdf, rr
  integer :: b

  if (ppcdf(0) < -1.0d0 .or. ppzero < -1.0d0) then
    write (*, *) '[Error] Wrong input CDF.'
    stop
  end if

  if (pp < ppzero_thres) then

    if (opt_pptrans == 2) then
      pos_cdf = ppzero * 0.5d0
    else
      write (*, *) '[Error] Unsupported transformation method.'
      stop
    end if

  else ! [pp >= ppzero_thres]

    if (pp < ppcdf(0)) then
      pos_cdf = 0.0d0
    else
      do b = 1, ncdf+1
        if (b > ncdf) then
          pos_cdf = 1.0d0
          exit
        else if (pp < ppcdf(b)) then
          rr = (pp - ppcdf(b-1)) / (ppcdf(b) - ppcdf(b-1))
          pos_cdf = ((1.0d0-rr) * real(b-1,r_size) + rr * real(b,r_size)) / real(ncdf,r_size)
          exit
        end if
      end do
    end if

  end if

  pptrans_normal = dinvnorm(compact_tail(pos_cdf))

end function pptrans_normal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function pptrans_log (pp)

  implicit none

  real(r_size) :: pptrans_log
  real(r_size), intent(in) :: pp

  if (pp < ppzero_thres) then
    pptrans_log = log(log_trans_tiny)
  else
    pptrans_log = log(pp + log_trans_tiny)
  end if

end function pptrans_log

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pptrans_normal_mdzero_def (pp_ens, ppcdf, ppzero, zero_mem, ym, sigma)

  implicit none

  real(r_size), intent(inout) :: pp_ens(nbv)
  real(r_size), intent(in)    :: ppcdf(0:ncdf)
  real(r_size), intent(in)    :: ppzero
  integer, intent(out)        :: zero_mem
  real(r_size), intent(out)   :: ym
  real(r_size), intent(out)   :: sigma

  real(r_size) :: pos_cdf
  real(r_size) :: ppzero_b, pprain_b
  real(r_size) :: y_trace, y_trace_b
  real(r_size) :: alpha, beta
  logical :: zero(nbv)
  integer :: n

!------
!  real(r_size) :: pp_ens_ori(nbv)
!  pp_ens_ori = pp_ens
!------

  if (ppcdf(0) < -1.0d0 .or. ppzero < -1.0d0) then
    write (*, *) '[Error] Wrong input CDF.'
    stop
  end if
  if (opt_pptrans /= 3) then
    write (*, *) '[Error] Unsupported transformation method.'
    stop
  end if

  beta = 0.0d0
  zero_mem = 0
  zero = .false.
  do n = 1, nbv
    if (pp_ens(n) < ppzero_thres) then
      zero_mem = zero_mem + 1
      zero(n) = .true.
    else
      pp_ens(n) = pptrans_normal(pp_ens(n), ppcdf, ppzero)
      beta = beta + pp_ens(n)
    end if
  end do
  beta = beta / real(nbv, r_size)
  ppzero_b = real(zero_mem, r_size) / real(nbv, r_size)
  pprain_b = 1.0d0 - ppzero_b

  y_trace = dinvnorm(compact_tail(ppzero))
  y_trace_b = dinvnorm(compact_tail(ppzero_b))

  alpha = 0.0d0 - exp(0.0d0 - 0.5d0*y_trace_b*y_trace_b) / sqrt(2.0d0*pi)
  ym = (alpha * y_trace + beta * y_trace_b) / (alpha + pprain_b * y_trace_b)
  sigma = (pprain_b * y_trace - beta) / (alpha + pprain_b * y_trace_b)

  do n = 1, nbv
    if (zero(n)) then
      pos_cdf = ppzero_b * 0.5d0
      pp_ens(n) = ym + sigma * dinvnorm(compact_tail(pos_cdf))
    end if
  end do

!------
!  print *, '----'
!  print *, ppzero, ppzero_b, dinvnorm(pos_cdf)
!  print *, y_trace, ym, sigma
!  do n = 1, nbv
!    print *, pp_ens_ori(n), pp_ens(n), zero(n)
!  end do
!------

end subroutine pptrans_normal_mdzero_def

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function pptrans_normal_mdzero (pp, ppcdf, ppzero, ppzero_m, zero_mem, ym, sigma)

  implicit none

  real(r_size) :: pptrans_normal_mdzero

  real(r_size), intent(in) :: pp
  real(r_size), intent(in) :: ppcdf(0:ncdf)
  real(r_size), intent(in) :: ppzero
  real(r_size), intent(in) :: ppzero_m
  integer, intent(in)      :: zero_mem
  real(r_size), intent(in) :: ym
  real(r_size), intent(in) :: sigma

  real(r_size) :: pos_cdf, rr
  integer :: b

  if (ppcdf(0) < -1.0d0 .or. ppzero < -1.0d0) then
    write (*, *) '[Error] Wrong input CDF.'
    stop
  end if
  if (opt_pptrans /= 3) then
    write (*, *) '[Error] Unsupported transformation method.'
    stop
  end if

  if (pp < ppzero_thres) then

    pos_cdf = ppzero * 0.5d0

  else ! [pp >= ppzero_thres]

    if (pp < ppcdf(0)) then
      pos_cdf = 0.0d0
    else
      do b = 1, ncdf+1
        if (b > ncdf) then
          pos_cdf = 1.0d0
          exit
        else if (pp < ppcdf(b)) then
          rr = (pp - ppcdf(b-1)) / (ppcdf(b) - ppcdf(b-1))
          pos_cdf = ((1.0d0-rr) * real(b-1,r_size) + rr * real(b,r_size)) / real(ncdf,r_size)
          exit
        end if
      end do
    end if

  end if

!------
!  print *, '---###'
!  print *, ppzero, ppzero_m, zero_mem, ym, sigma
!  print *, pos_cdf
!------

  if (pos_cdf < ppzero_m) then
    pos_cdf = (pos_cdf / ppzero_m) * (real(zero_mem, r_size) / real(nbv, r_size))
    pptrans_normal_mdzero = ym + sigma * dinvnorm(compact_tail(pos_cdf))
  else
    pptrans_normal_mdzero = dinvnorm(compact_tail(pos_cdf))
  end if

!------
!  print *, pos_cdf
!  print *, pptrans_normal_mdzero
!------

end function pptrans_normal_mdzero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function compact_tail (pos_cdf)

  implicit none

  real(r_size) :: compact_tail
  real(r_size), intent(in) :: pos_cdf

  compact_tail = pos_cdf
  if (compact_tail < gausstail_thres        ) compact_tail = gausstail_thres
  if (compact_tail > 1.0d0 - gausstail_thres) compact_tail = 1.0d0 - gausstail_thres

end function compact_tail

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ren-raw chen, rutgers business school
! normal inverse
! translate from http://home.online.no/~pjacklam/notes/invnorm
! a routine written by john herrero

real*8 function dinvnorm(p)
      real*8 p,p_low,p_high
      real*8 a1,a2,a3,a4,a5,a6
      real*8 b1,b2,b3,b4,b5
      real*8 c1,c2,c3,c4,c5,c6
      real*8 d1,d2,d3,d4
      real*8 z,q,r
      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1-p_low
      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
201   q=dsqrt(-2*dlog(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
      goto 204
301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
202   q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      goto 204
302   if((p.gt.p_high).and.(p.lt.1)) goto 203
203   q=dsqrt(-2*dlog(1-p))
      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
204   dinvnorm=z
      return
end function dinvnorm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module common_precip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
