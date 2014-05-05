!===============================================================================
module superob_tools
!-------------------------------------------------------------------------------

  use common
  use common_gfs
  use common_obs_gfs
  use obs_tools
  use obs_gfs_ext, only: nlevl
  implicit none

contains

!===============================================================================
! method:
!  0: do nothing
!  1: leave only one observation with minimum errors and closest to the grid center
!  2: average all observations
!  3: average observations with minimum errors
!-------------------------------------------------------------------------------
subroutine superob_select (naobs, i, j, k, ai, aj, ak, alon, alat, alev, adat, aerr, method)
!-------------------------------------------------------------------------------

  implicit none
  integer, intent(inout) :: naobs
  integer, intent(in) :: i, j, k
  real(r_size), intent(inout) :: ai(naobs)
  real(r_size), intent(inout) :: aj(naobs)
  real(r_size), intent(inout) :: ak(naobs)
  real(r_size), intent(inout) :: alon(naobs)
  real(r_size), intent(inout) :: alat(naobs)
  real(r_size), intent(inout) :: alev(naobs)
  real(r_size), intent(inout) :: adat(naobs)
  real(r_size), intent(inout) :: aerr(naobs)
  integer, intent(in) :: method

  integer :: idx(naobs)
  integer :: n, n_use, n_min
  real(r_size) :: dx, dy, dz, dist, dist_min, err_min
  real(r_size) :: aai, aaj, aak, aalon, aalat, aalev, aadat

  if (naobs <= 1) return
  if (method <= 0 .or. method >= 4) return

!-------------------------------------------------------------------------------

  if (method == 2) then
    n_use = naobs
    err_min = 1.e20  
    do n = 1, naobs
      idx(n) = n
      if (aerr(n) < err_min) err_min = aerr(n)
    end do
  else if (method == 1 .or. method == 3) then
    n_use = 0
    err_min = 1.e20  
    do n = 1, naobs
      if (aerr(n) == err_min) then
        n_use = n_use + 1
        idx(n_use) = n
      else if (aerr(n) < err_min) then
        n_use = 1
        idx(n_use) = n
        err_min = aerr(n)
      end if
    end do

    if (method == 1) then
      dist_min = 1.e20
      do n = 1, n_use
        dx = ai(idx(n)) - real(i, r_size)
        dy = aj(idx(n)) - real(j, r_size)
        if (k == 0) then
          dist = sqrt(dx*dx + dy*dy)
        else
          dz = ak(idx(n)) - real(k, r_size)
          dist = sqrt(dx*dx + dy*dy + dz*dz)
        end if
        if (dist < dist_min) then
          n_min = idx(n)
          dist_min = dist
        end if
      end do
      n_use = 1
      idx(1) = n_min
    end if
  end if

!-------------------------------------------------------------------------------
!  print *, '===', i, j, k
!  print *, ai(:)
!  print *, aj(:)
!  print *, ak(:)
!  print *, aerr(:)
!  print *, '---'
!  print *, n_use
!  print *, idx(1:n_use)

  naobs = 1

  if (n_use > 1) then
    aai = 0.d0
    aaj = 0.d0
    aak = 0.d0
    aalon = 0.d0
    aalat = 0.d0
    aalev = 0.d0
    aadat = 0.d0
    do n = 1, n_use
      aai   = aai   + ai  (idx(n))
      aaj   = aaj   + aj  (idx(n))
      aak   = aak   + ak  (idx(n))
      aalon = aalon + alon(idx(n))
      aalat = aalat + alat(idx(n))
      aalev = aalev + alev(idx(n))
      aadat = aadat + adat(idx(n))
    end do
    ai  (1) = aai   / real(n_use, r_size)
    aj  (1) = aaj   / real(n_use, r_size)
    ak  (1) = aak   / real(n_use, r_size)
    alon(1) = aalon / real(n_use, r_size)
    alat(1) = aalat / real(n_use, r_size)
    alev(1) = aalev / real(n_use, r_size)
    adat(1) = aadat / real(n_use, r_size)
    aerr(1) = err_min
  else if (n_use == 1) then
    ai  (1) = ai  (idx(1))
    aj  (1) = aj  (idx(1))
    ak  (1) = ak  (idx(1))
    alon(1) = alon(idx(1))
    alat(1) = alat(idx(1))
    alev(1) = alev(idx(1))
    adat(1) = adat(idx(1))
    aerr(1) = aerr(idx(1))
  end if

!  print *, ai(1), aj(1), ak(1), aerr(1)

!-------------------------------------------------------------------------------
end subroutine superob_select
!===============================================================================

!-------------------------------------------------------------------------------
subroutine print_obsnum_cal (naobs, aelm, atyp)
!-------------------------------------------------------------------------------

  implicit none
  integer, intent(in) :: naobs
  integer, intent(in) :: aelm(naobs)
  integer, intent(in) :: atyp(naobs)

  integer :: obscount(nid_obs,nobtype+1)
  integer :: obstotal
  integer :: itype, n

  obscount = 0
  obstotal = 0

  do n = 1, naobs
    obscount(aelm(n),atyp(n)) = obscount(aelm(n),atyp(n)) + 1
    obstotal = obstotal + 1
  end do

  call print_obsnum(obstotal, obscount)

!-------------------------------------------------------------------------------
end subroutine print_obsnum_cal
!===============================================================================

!-------------------------------------------------------------------------------
subroutine print_obsnum_cal_qc (naobs, aelm, atyp, aqc)
!-------------------------------------------------------------------------------

  implicit none
  integer, intent(in) :: naobs
  integer, intent(in) :: aelm(naobs)
  integer, intent(in) :: atyp(naobs)
  integer, intent(in) :: aqc(naobs)

  integer :: obscount(nid_obs,nobtype+1)
  integer :: obstotal
  integer :: itype, n

  obscount = 0
  obstotal = 0

  do n = 1, naobs
    if (aqc(n) == 0) then
      obscount(aelm(n),atyp(n)) = obscount(aelm(n),atyp(n)) + 1
      obstotal = obstotal + 1
    end if
  end do
 
  call print_obsnum(obstotal, obscount)

!-------------------------------------------------------------------------------
end subroutine print_obsnum_cal_qc
!===============================================================================

!-------------------------------------------------------------------------------
subroutine write_obs_den (naobs, aii, ajj, akk, aelm, atype, aqc, filename, qctype)
!-------------------------------------------------------------------------------

  implicit none
  integer, intent(in) :: naobs
  integer, intent(in) :: aii(naobs)
  integer, intent(in) :: ajj(naobs)
  integer, intent(in) :: akk(naobs)
  integer, intent(in) :: aelm(naobs)
  integer, intent(in) :: atype(naobs)
  integer, intent(in) :: aqc(naobs)
  character (len=*), intent(in) :: filename
  integer, intent(in) :: qctype  ! 1: before superob; 2: after superob

  integer :: nobsden(nlon,nlat,nlev+1,nid_obs+1)
  integer :: nobsden_sum(nlon,nlat,nlev+1,nid_obs+1)
  real(r_sngl) :: nobsden_r(nlon,nlat,nlev+1,nid_obs+1)
  integer :: iolen, reclen
  integer :: n, j, k, itype

  inquire (iolength = iolen) iolen
  reclen = iolen * nlon * nlat * (nlev+1) * (nid_obs+1)
  open (21, file = trim(filename), form = 'unformatted', access = 'direct', recl = reclen)

  nobsden_sum(:,:,:,:) = 0
  do itype = 1, nobtype
    nobsden(:,:,:,:) = 0
    do n = 1, naobs
      if (atype(n) /= itype) cycle
      if (qctype == 1 .and. aqc(n) /= 0 .and. aqc(n) /= 9) cycle
      if (qctype == 2 .and. aqc(n) /= 0) cycle
      k = akk(n)
      if (k > 0) then
        k = k - nlevl + 1
        if (k < 2) k = 2
      else
        k = 1
      end if
      nobsden(aii(n),ajj(n),k,aelm(n)  ) = nobsden(aii(n),ajj(n),k,aelm(n)  ) + 1
      nobsden(aii(n),ajj(n),k,nid_obs+1) = nobsden(aii(n),ajj(n),k,nid_obs+1) + 1
    end do
    nobsden_r = real(nobsden, 4)
    write (21, rec = itype) nobsden_r
    nobsden_sum = nobsden_sum + nobsden
  end do

  nobsden_r = real(nobsden_sum, 4)
  write (21, rec = nobtype+1) nobsden_r
  close (21)
    
!-------------------------------------------------------------------------------

  open (22, file = trim(filename) // '.ctl', form = 'formatted', access = 'sequential')
  write (22, '(2A)') 'dset ^', filename
  write (22, '(A)') 'options byteswapped'
  write (22, '(A)') 'undef -9.99E+33'
  write (22, '(A,I6,A,2F12.6)') 'xdef', nlon, ' linear', 0., 360./real(nlon)
  write (22, '(A,I6,A)') 'ydef', nlat, ' levels'
  do j = 1, nlat 
    write (22, '(F12.6)') lat(j)
  end do
  write (22, '(A,I6,A)') 'zdef', nlev+1, ' linear 1 1'
  write (22, '(A,I6,A)') 'tdef', nid_obs+1, ' linear 00z01jan2000 1dy'
  write (22, '(A,I6,A)') 'edef', nobtype+1, ' names'
  do itype = 1, nobtype
    write (22, '(A)') obtypelist(itype)
  end do
  write (22, '(A)') 'TOTAL'
  write (22, '(A)') 'vars     1'
  write (22, '(A,I6,A)') 'den  ', nlev+1, ' 99 observation number'
  write (22, '(A)') 'endvars'
  close (22)

!-------------------------------------------------------------------------------
end subroutine write_obs_den
!===============================================================================

end module superob_tools
!===============================================================================
