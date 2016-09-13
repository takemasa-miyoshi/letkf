!===============================================================================
! superob: main program
!-------------------------------------------------------------------------------
!
!  Procedure of the algorithm:
!   1. general removal
!   2. vertical (only valid for multi-level observations) superob
!   3. temporal superob
!   4. all grid (horizontal+vertical) superob
!
!  Refer variables obmethod_g, obmethod_v, obmethod_t, and obmethod_h
!  for superob method settings
!
!===============================================================================
program superob
!-------------------------------------------------------------------------------

  use common
  use common_gfs
  use common_obs_gfs
  use superob_tools
  use obs_gfs_ext
  use sigio_module
  implicit none

  integer,parameter :: nslots = 7 ! number of time slots for 4d-letkf
  integer,parameter :: nbslot = 4 ! basetime slot

  character(9)  :: obsfile = 'obsTT.dat'
  character(9)  :: supfile = 'supTT.dat'
  character(10) :: obsdenfile = 'obsden.grd'
  character(13) :: obsdenfile_so = 'obsden_so.grd'
  character(11) :: refpresfile = 'refpres.dat'

  logical, parameter :: print_log = .true.
  logical, parameter :: write_density_grid = .false.

  real(r_sngl) :: refv3d(nlon,nlat,nlev,nv3d)
  real(r_sngl) :: refv2d(nlon,nlat,nv2d)
  real(r_size) :: tmpps(nlon*nlat)
  real(r_size) :: tmptv(nlon*nlat,nlev)
  real(r_size) :: tmpp(nlon*nlat,nlev)
  real(r_size) :: refpres(nlon,nlat,nlevx)

  real(r_size), allocatable :: tmpelm(:)
  real(r_size), allocatable :: tmplon(:)
  real(r_size), allocatable :: tmplat(:)
  real(r_size), allocatable :: tmplev(:)
  real(r_size), allocatable :: tmpdat(:)
  real(r_size), allocatable :: tmperr(:)
  real(r_size), allocatable :: tmptyp(:)
  integer,      allocatable :: tmpelmi(:)
  integer,      allocatable :: tmptypi(:)
  integer,      allocatable :: tmploti(:)
  real(r_size), allocatable :: tmpri(:)
  real(r_size), allocatable :: tmprj(:)
  real(r_size), allocatable :: tmprk(:)
  integer,      allocatable :: tmpi(:)
  integer,      allocatable :: tmpj(:)
  integer,      allocatable :: tmpk(:)
  integer,      allocatable :: tmpqc(:)

  real(r_size), allocatable :: tmp2lon(:)
  real(r_size), allocatable :: tmp2lat(:)
  real(r_size), allocatable :: tmp2lev(:)
  real(r_size), allocatable :: tmp2dat(:)
  real(r_size), allocatable :: tmp2err(:)
  integer,      allocatable :: tmp2elmi(:)
  integer,      allocatable :: tmp2typi(:)
  integer,      allocatable :: tmp2loti(:)
  real(r_size), allocatable :: tmp2ri(:)
  real(r_size), allocatable :: tmp2rj(:)
  real(r_size), allocatable :: tmp2rk(:)
  integer,      allocatable :: tmp2i(:)
  integer,      allocatable :: tmp2j(:)
  integer,      allocatable :: tmp2k(:)
  integer,      allocatable :: tmp2qc(:)

  integer, parameter :: max_platform_reclen = 1000
  real(r_size) :: tmp3lon(max_platform_reclen)
  real(r_size) :: tmp3lat(max_platform_reclen)
  real(r_size) :: tmp3lev(max_platform_reclen)
  real(r_size) :: tmp3dat(max_platform_reclen)
  real(r_size) :: tmp3err(max_platform_reclen)
  real(r_size) :: tmp3ri (max_platform_reclen)
  real(r_size) :: tmp3rj (max_platform_reclen)
  real(r_size) :: tmp3rk (max_platform_reclen)
  real(r_size) :: tmp3k  (max_platform_reclen)

  integer :: nobsgrd(nlon,nlat,0:nlevx,nslots,nid_obs)
  integer :: nobsgrdac(0:nlon,nlat,0:nlevx,nslots,nid_obs)
  integer :: iobsgrd(nlon,nlat,0:nlevx,nslots,nid_obs)

  integer :: obmethod_g(nobtype,nid_obs)
  integer :: obmethod_h(nobtype,nid_obs)
  integer :: obmethod_v(nobtype,nid_obs)
  integer :: obmethod_t(nobtype,nid_obs)

  integer :: nobs, nobs_thistype
  integer :: nobslots(nslots)
  integer :: nobsrdc_bd, nobsrdc_g, nobsrdc_tp, nobsrdc_v, nobsrdc_t, nobsrdc_h

  integer :: tdist, tdist_min
  integer :: slot_priority(nslots)
  logical :: slotused(nslots)

  integer :: i, j, k, islot, itype, ielm, accu
  integer :: n, nn, nnn, m, mm
  integer :: n1, n2, n3, n4, nnn_before, nnn_after
  integer :: is, is2, islot2, m2
  integer :: recstart, recend
  integer :: iret
  logical :: cont
  real(r_size) :: levout, datout, errout

!-------------------------------------------------------------------------------
! -- general
! 0: do nothing
! 1: do not use observations with these type/variables
  obmethod_g(:,:) = 0

! -- vertical (only for multi-level observations / require same observation records are stored continuously)
! 0: do nothing
! 1: leave only one observation with minimum errors and closest to the grid center
! 2: average all observations
! 3: average observations with minimum errors
  obmethod_v(:,:) = 0
  obmethod_v(1,:) = 2
  obmethod_v(5,:) = 2
  obmethod_v(6,:) = 2

! -- temporal
! 0: do nothing
! 1: leave only one observation closest to the base time
! 2: leave only one observation at the base time and do not use any observations off the base time
  obmethod_t(:,:) = 0
  obmethod_t(5,:) = 1
  obmethod_t(6,:) = 1
  obmethod_t(8,:) = 1
  obmethod_t(9,:) = 2
  obmethod_t(10,:) = 2

! -- all grid (horizontal+vertical)
! 0: do nothing
! 1: leave only one observation with minimum errors and closest to the grid center
! 2: average all observations
! 3: average observations with minimum errors
  obmethod_h(:,:) = 3
  obmethod_h(1,:) = 0
  obmethod_h(5,:) = 0
  obmethod_h(6,:) = 0

!===============================================================================

  call set_common_gfs

  nobsrdc_bd = 0
  nobsrdc_g = 0
  nobsrdc_tp = 0
  nobsrdc_h = 0
  nobsrdc_v = 0
  nobsrdc_t = 0

  slotused(:) = .false.
  do islot = 1, nslots
    tdist_min = 99999999
    do is = 1, nslots
      tdist = abs(is - nbslot)
      if ((.not. slotused(is)) .and. tdist < tdist_min) then
        slot_priority(islot) = is
        tdist_min = tdist
      end if
    end do
    slotused(slot_priority(islot)) = .true.
  end do
  
!-------------------------------------------------------------------------------

  open (81, file = 'log.1.type',     form = 'formatted', access = 'sequential')
  open (82, file = 'log.2.boundary', form = 'formatted', access = 'sequential')
  open (83, file = 'log.3.vertical', form = 'formatted', access = 'sequential')
  open (84, file = 'log.4.temporal', form = 'formatted', access = 'sequential')
  open (85, file = 'log.5.allgrid',  form = 'formatted', access = 'sequential')

!-------------------------------------------------------------------------------

  call read_grd4(refpresfile,refv3d,refv2d,0)

  tmpps = reshape(refv2d(:,:,iv2d_ps),(/nlon*nlat/))
  tmptv = reshape(refv3d(:,:,:,iv3d_t) * (1.0d0 + fvirt * refv3d(:,:,:,iv3d_q)),(/nlon*nlat,nlev/))
  call sigio_modprd(nlon*nlat,nlon*nlat,nlev,gfs_nvcoord,gfs_idvc,gfs_idsl, &
                    gfs_vcoord,iret,tmpps,tmptv,pm=tmpp)
  refpres(:,:,1+nlevl:nlev+nlevl) = reshape(tmpp,(/nlon,nlat,nlev/))

  do k = nlevl, 1, -1
    refpres(:,:,k) = 2. * refpres(:,:,k+1) - refpres(:,:,k+2)
  end do

!-------------------------------------------------------------------------------

  do islot=1,nslots
    write (obsfile(4:5), '(I2.2)') islot
    call get_nobs(obsfile, 7, nobslots(islot))
  end do
  nobs = sum(nobslots)

!-------------------------------------------------------------------------------

  allocate ( tmpelm (nobs) )
  allocate ( tmplon (nobs) )
  allocate ( tmplat (nobs) )
  allocate ( tmplev (nobs) )
  allocate ( tmpdat (nobs) )
  allocate ( tmperr (nobs) )
  allocate ( tmptyp (nobs) )
  allocate ( tmpelmi(nobs) )
  allocate ( tmptypi(nobs) )
  allocate ( tmploti(nobs) )
  allocate ( tmpri  (nobs) )
  allocate ( tmprj  (nobs) )
  allocate ( tmprk  (nobs) )
  allocate ( tmpi   (nobs) )
  allocate ( tmpj   (nobs) )
  allocate ( tmpk   (nobs) )
  allocate ( tmpqc  (nobs) )

  allocate ( tmp2lon (nobs) )
  allocate ( tmp2lat (nobs) )
  allocate ( tmp2lev (nobs) )
  allocate ( tmp2dat (nobs) )
  allocate ( tmp2err (nobs) )
  allocate ( tmp2elmi(nobs) )
  allocate ( tmp2typi(nobs) )
  allocate ( tmp2loti(nobs) )
  allocate ( tmp2ri  (nobs) )
  allocate ( tmp2rj  (nobs) )
  allocate ( tmp2rk  (nobs) )
  allocate ( tmp2i   (nobs) )
  allocate ( tmp2j   (nobs) )
  allocate ( tmp2k   (nobs) )
  allocate ( tmp2qc  (nobs) )

!===============================================================================

  write (*, *)
  write (*, '(A)') '================================================'
  write (*, '(A)') ' read observations and general/vertical superob '
  write (*, '(A)') '================================================'

  tmpqc(:) = 0
  nn = 0

!------------------------------
timeslots: do islot = 1, nslots
!------------------------------

    if (nobslots(islot) == 0) cycle

    write (obsfile(4:5), '(I2.2)') islot
    write (*, '(2A)') 'read obs file: ', obsfile
    call read_obs(obsfile, nobslots(islot), &
                  tmpelm(nn+1:nn+nobslots(islot)), tmplon(nn+1:nn+nobslots(islot)), &
                  tmplat(nn+1:nn+nobslots(islot)), tmplev(nn+1:nn+nobslots(islot)), &
                  tmpdat(nn+1:nn+nobslots(islot)), tmperr(nn+1:nn+nobslots(islot)), &
                  tmptyp(nn+1:nn+nobslots(islot)))
    tmploti(nn+1:nn+nobslots(islot)) = islot
    do nnn = nn+1, nn+nobslots(islot)
      tmpelmi(nnn) = uid_obs(nint(tmpelm(nnn)))
      tmptypi(nnn) = nint(tmptyp(nnn))
    end do

    !-----------------------------------------------------------------------

    cont = .false.
    recstart = nn + 1
    do n = 1, nobslots(islot)

      call phys2ijk_ext(refpres, tmpelm(nn+n), tmplon(nn+n), tmplat(nn+n), tmplev(nn+n), &
                        tmpri(nn+n), tmprj(nn+n), tmprk(nn+n))

      tmpi(nn+n) = floor(tmpri(nn+n)+0.5)
      tmpj(nn+n) = floor(tmprj(nn+n)+0.5)
      tmpk(nn+n) = floor(tmprk(nn+n)+0.5)

      if (tmpi(nn+n) < 1) tmpi(nn+n) = tmpi(nn+n) + nlon
      if (tmpi(nn+n) > nlon) tmpi(nn+n) = tmpi(nn+n) - nlon

      ! surface observation types (ADPSFC/SFCSHP/SFCBOG/SPSSMI/QKSWND)
      if (tmptypi(nn+n) == 8 .or. tmptypi(nn+n) == 9 .or. tmptypi(nn+n) == 10 .or. &
          tmptypi(nn+n) == 11 .or. tmptypi(nn+n) == 15) then
        tmpk(nn+n) = 0
      ! surface observation variables (PS/PP/TCobs)
      else if (nint(tmpelm(nn+n)) > 9999) then
        tmpk(nn+n) = 0
      else
        if (tmpk(nn+n) < 1) tmpk(nn+n) = 1
        if (ceiling(tmprk(nn+n)) > nlevx) then
          tmpqc(nn+n) = 3
          nobsrdc_bd = nobsrdc_bd + 1
          if (print_log) then
            write (82, '(A,A6,A,A3,A,I3,3(A,F7.2),A,F7.2)') &
              '[', obtypelist(tmptypi(nn+n)), ':', obelmlist(tmpelmi(nn+n)), &
              '] slot=', islot, ' lon=', tmplon(nn+n), ' lat=', tmplat(nn+n), ' lev=', tmplev(nn+n)/100., &
              ' | obs removed because too high'
          end if
        end if
      end if
      
      if (ceiling(tmprj(nn+n)) < 2 .and. tmpqc(nn+n) == 0) then
        tmpqc(nn+n) = 2
        nobsrdc_bd = nobsrdc_bd + 1
        if (print_log) then
          if (tmpk(nn+n) == 0) then
            write (82, '(A,A6,A,A3,A,I3,3(A,F7.2),A,F7.2)') &
              '[', obtypelist(tmptypi(nn+n)), ':', obelmlist(tmpelmi(nn+n)), &
              '] slot=', islot, ' lon=', tmplon(nn+n), ' lat=', tmplat(nn+n), ' alt=', tmplev(nn+n), &
              ' | obs removed because lat <', lat(1)
          else
            write (82, '(A,A6,A,A3,A,I3,3(A,F7.2),A,F7.2)') &
              '[', obtypelist(tmptypi(nn+n)), ':', obelmlist(tmpelmi(nn+n)), &
              '] slot=', islot, ' lon=', tmplon(nn+n), ' lat=', tmplat(nn+n), ' lev=', tmplev(nn+n)/100., &
              ' | obs removed because lat <', lat(1)
          end if
        end if
      end if
      if (ceiling(tmprj(nn+n)) > nlat .and. tmpqc(nn+n) == 0) then
        tmpqc(nn+n) = 2
        nobsrdc_bd = nobsrdc_bd + 1
        if (print_log) then
          if (tmpk(nn+n) == 0) then
            write (82, '(A,A6,A,A3,A,I3,3(A,F7.2),A,F7.2)') &
              '[', obtypelist(tmptypi(nn+n)), ':', obelmlist(tmpelmi(nn+n)), &
              '] slot=', islot, ' lon=', tmplon(nn+n), ' lat=', tmplat(nn+n), ' alt=', tmplev(nn+n), &
              ' | obs removed because lat >', lat(nlat)
          else
            write (82, '(A,A6,A,A3,A,I3,3(A,F7.2),A,F7.2)') &
              '[', obtypelist(tmptypi(nn+n)), ':', obelmlist(tmpelmi(nn+n)), &
              '] slot=', islot, ' lon=', tmplon(nn+n), ' lat=', tmplat(nn+n), ' lev=', tmplev(nn+n)/100., &
              ' | obs removed because lat >', lat(nlat)
          end if
        end if
      end if

      if (obmethod_g(tmptypi(nn+n),tmpelmi(nn+n)) == 1 .and. tmpqc(nn+n) == 0) then
        tmpqc(nn+n) = 4
        nobsrdc_g = nobsrdc_g + 1
        if (print_log) then
          if (tmpk(nn+n) == 0) then
            write (81, '(A,A6,A,A3,A,I3,3(A,F7.2),A,F7.2)') &
              '[', obtypelist(tmptypi(nn+n)), ':', obelmlist(tmpelmi(nn+n)), &
              '] slot=', islot, ' lon=', tmplon(nn+n), ' lat=', tmplat(nn+n), ' alt=', tmplev(nn+n), &
              ' | obs removed because not to use this type/variables'
          else
            write (81, '(A,A6,A,A3,A,I3,3(A,F7.2),A,F7.2)') &
              '[', obtypelist(tmptypi(nn+n)), ':', obelmlist(tmpelmi(nn+n)), &
              '] slot=', islot, ' lon=', tmplon(nn+n), ' lat=', tmplat(nn+n), ' lev=', tmplev(nn+n)/100., &
              ' | obs removed because not to use this type/variables'
          end if
        end if
      end if

      if ((tmpelmi(nn+n) < 1 .or. tmpelmi(nn+n) > nid_obs .or. &
           tmptypi(nn+n) < 1 .or. tmptypi(nn+n) > nobtype) .and. tmpqc(nn+n) == 0) then
        tmpqc(nn+n) = 5
        nobsrdc_tp = nobsrdc_tp + 1
        if (print_log) then
          if (tmpk(nn+n) == 0) then
            write (81, '(A,I3,A,I6,A,I3,3(A,F7.2),A)') &
              'typeid=', tmptypi(nn+n), ' varid=', nint(tmpelm(nn+n)), &
              ' slot=', islot, ' lon=', tmplon(nn+n), ' lat=', tmplat(nn+n), ' alt=', tmplev(nn+n), &
              ' | obs removed because of unsupported type/variables'
          else
            write (81, '(A,I3,A,I6,A,I3,3(A,F7.2),A)') &
              'typeid=', tmptypi(nn+n), ' varid=', nint(tmpelm(nn+n)), &
              ' slot=', islot, ' lon=', tmplon(nn+n), ' lat=', tmplat(nn+n), ' lev=', tmplev(nn+n)/100., &
              ' | obs removed because of unsupported type/variables'
          end if
        end if
      end if

    !-----------------------------------------------------------------------

      if (n > 1) then
        if (obmethod_v(tmptypi(nn+n),tmpelmi(nn+n)) > 0 .and. &
            tmplon(nn+n) == tmplon(recstart) .and. tmplat(nn+n) == tmplat(recstart) .and. &
            tmpelmi(nn+n) == tmpelmi(recstart) .and. tmptypi(nn+n) == tmptypi(recstart)) then
          cont = .true.
        else

    !---- vertical superob for one multi-level record ----------------------
          if (cont) then
            recend = nn + n - 1

            mm = 0
            do k = 0, nlevx
              m = 0
              do nnn = recstart, recend
                if (tmpk(nnn) == k .and. tmpqc(nnn) == 0) then
                  m = m + 1
                  tmp3lon(mm+m) = tmplon(nnn)
                  tmp3lat(mm+m) = tmplat(nnn)
                  tmp3lev(mm+m) = tmplev(nnn)
                  tmp3dat(mm+m) = tmpdat(nnn)
                  tmp3err(mm+m) = tmperr(nnn)
                  tmp3ri (mm+m) = tmpri (nnn)
                  tmp3rj (mm+m) = tmprj (nnn)
                  tmp3rk (mm+m) = tmprk (nnn)
                  tmp3k  (mm+m) = k
                end if
              end do
              if (m > 1) then
                n1 = mm + 1
                n2 = mm + m
                i = tmpi(recstart)
                j = tmpj(recstart)
                call superob_select(m, i, j, k, &
                                    tmp3ri(n1:n2) , tmp3rj(n1:n2) , tmp3rk(n1:n2) , &
                                    tmp3lon(n1:n2), tmp3lat(n1:n2), tmp3lev(n1:n2), &                                       
                                    tmp3dat(n1:n2), tmp3err(n1:n2), &
                                    obmethod_v(tmptypi(recstart),tmpelmi(recstart)))
              end if
              mm = mm + m
            end do

            tmpqc(recstart:recend) = 9
            if (mm > 0) then
              n3 = recstart
              n4 = recstart + mm - 1
              tmplon(n3:n4) = tmp3lon(1:mm)
              tmplat(n3:n4) = tmp3lat(1:mm)
              tmplev(n3:n4) = tmp3lev(1:mm)
              tmpdat(n3:n4) = tmp3dat(1:mm)
              tmperr(n3:n4) = tmp3err(1:mm)
              tmpri (n3:n4) = tmp3ri (1:mm)
              tmprj (n3:n4) = tmp3rj (1:mm)
              tmprk (n3:n4) = tmp3rk (1:mm)
              tmpk  (n3:n4) = tmp3k  (1:mm)
              tmpqc (n3:n4) = 0
            end if

            if (mm < recend-recstart+1) then
              nobsrdc_v = nobsrdc_v + recend-recstart+1 - mm
              if (print_log) then
                write (83, '(A,A6,A,A3,A,I3,4(A,F7.2),2(A,I4))') &
                  '[', obtypelist(tmptypi(nnn)), ':', obelmlist(tmpelmi(nnn)), &
                  '] slot=', islot, ' lon=', tmplon(recstart), ' lat=', tmplat(recstart), &
                  ' lev=', maxval(tmplev(recstart:recend))/100., '~', minval(tmplev(recstart:recend))/100., &
                  ' | obs # reduced from', recend-recstart+1, ' ->', mm
              end if
            end if

            cont = .false.
          end if ! (cont)
          recstart = nn + n
    !-----------------------------------------------------------------------

        end if
      end if ! [ n > 1 ]

    end do ! [ n = 1, nobslots(islot) ]

    !---- repeat vertical superob for one multi-level record once at the end of the loop (if cont)
    if (cont) then
      recend = nn + n

      mm = 0
      do k = 0, nlevx
        m = 0
        do nnn = recstart, recend
          if (tmpk(nnn) == k .and. tmpqc(nnn) == 0) then
            m = m + 1
            tmp3lon(mm+m) = tmplon(nnn)
            tmp3lat(mm+m) = tmplat(nnn)
            tmp3lev(mm+m) = tmplev(nnn)
            tmp3dat(mm+m) = tmpdat(nnn)
            tmp3err(mm+m) = tmperr(nnn)
            tmp3ri (mm+m) = tmpri (nnn)
            tmp3rj (mm+m) = tmprj (nnn)
            tmp3rk (mm+m) = tmprk (nnn)
            tmp3k  (mm+m) = k
          end if
        end do
        if (m > 1) then
          n1 = mm + 1
          n2 = mm + m
          i = tmpi(recstart)
          j = tmpj(recstart)
          call superob_select(m, i, j, k, &
                              tmp3ri(n1:n2) , tmp3rj(n1:n2) , tmp3rk(n1:n2) , &
                              tmp3lon(n1:n2), tmp3lat(n1:n2), tmp3lev(n1:n2), &                                       
                              tmp3dat(n1:n2), tmp3err(n1:n2), &
                              obmethod_v(tmptypi(recstart),tmpelmi(recstart)))
        end if
        mm = mm + m
      end do

      tmpqc(recstart:recend) = 9
      if (mm > 0) then
        n3 = recstart
        n4 = recstart + mm - 1
        tmplon(n3:n4) = tmp3lon(1:mm)
        tmplat(n3:n4) = tmp3lat(1:mm)
        tmplev(n3:n4) = tmp3lev(1:mm)
        tmpdat(n3:n4) = tmp3dat(1:mm)
        tmperr(n3:n4) = tmp3err(1:mm)
        tmpri (n3:n4) = tmp3ri (1:mm)
        tmprj (n3:n4) = tmp3rj (1:mm)
        tmprk (n3:n4) = tmp3rk (1:mm)
        tmpk  (n3:n4) = tmp3k  (1:mm)
        tmpqc (n3:n4) = 0
      end if
      
      if (mm < recend-recstart+1) then
        nobsrdc_v = nobsrdc_v + recend-recstart+1 - mm
        if (print_log) then
          write (83, '(A,A6,A,A3,A,I3,4(A,F7.2),2(A,I4))') &
            '[', obtypelist(tmptypi(nnn)), ':', obelmlist(tmpelmi(nnn)), &
            '] slot=', islot, ' lon=', tmplon(recstart), ' lat=', tmplat(recstart), &
            ' lev=', maxval(tmplev(recstart:recend))/100., '~', minval(tmplev(recstart:recend))/100., &
            ' | obs # reduced from', recend-recstart+1, ' ->', mm
        end if
      end if
    end if ! (cont)
    !-----------------------------------------------------------------------

    nn = nn + nobslots(islot)

  !---------------
  end do timeslots ! [ islot = 1, nslots ]
  !---------------

! tmpqc:
!  0: pass qc
!  2: dropped because of latitude
!  3: dropped because of pressure
!  4: dropped because not to use these type/variables
!  5: dropped because of unsupported type/variables
!  9: dropped because of vertical superob
!-------------------------------------------------------------------------------

  write (*, *)
  write (*, '(A)') 'statistics of original obs file:'

  call print_obsnum_cal(nobs, tmpelmi, tmptypi)

  if (write_density_grid) then
    write (*, *)
    write (*, '(2A)') 'write observation density file before superob: ', obsdenfile
    call write_obs_den(nobs, tmpi, tmpj, tmpk, tmpelmi, tmptypi, tmpqc, obsdenfile, 1)
  end if ! [ write_density_grid ]

!-------------------------------------------------------------------------------

  write (*, *)
  write (*, '(A,I10)') 'observation number reduced because outside model boundary:', nobsrdc_bd
  write (*, '(A,I10)') 'observation number reduced because not to use these type/variables:', nobsrdc_g
  write (*, '(A,I10)') 'observation number reduced because unsupported type:', nobsrdc_tp
  write (*, '(A,I10)') 'observation number reduced by vertical superob:', nobsrdc_v
  write (*, *)
  write (*, '(A)') 'statistics after general and vertical superob:'

  call print_obsnum_cal_qc(nobs, tmpelmi, tmptypi, tmpqc)

!===============================================================================

  write (*, *)
  write (*, '(A)') '================================================'
  write (*, '(A)') ' temporal superob                               '
  write (*, '(A)') '================================================'

  tmp2qc(:) = 0
  nnn = 0

  !--------------------
  do itype = 1, nobtype
  !--------------------

    write (*, '(2A)') 'processing ', obtypelist(itype)

    nobsgrd(:,:,:,:,:) = 0
    nobsgrdac(:,:,:,:,:) = 0
    iobsgrd(:,:,:,:,:) = 0

    do n = 1, nobs
      if (tmptypi(n) /= itype) cycle
      if (tmpqc(n) == 0) then
        nobsgrd(tmpi(n),tmpj(n),tmpk(n),tmploti(n),tmpelmi(n)) &
          = nobsgrd(tmpi(n),tmpj(n),tmpk(n),tmploti(n),tmpelmi(n)) + 1
      end if
    end do

    accu = 0
    do ielm = 1, nid_obs
      do islot = 1, nslots
        do k = 0, nlevx
          do j = 1, nlat
            nobsgrdac(0,j,k,islot,ielm) = accu
            do i = 1, nlon
              accu = accu + nobsgrd(i,j,k,islot,ielm)
              nobsgrdac(i,j,k,islot,ielm) = accu
            end do
          end do
        end do
      end do
    end do
    nobs_thistype = accu

    !-----------------------------------------------------------------------

    if (nobs_thistype > 0) then

      do n = 1, nobs
        if (tmptypi(n) /= itype) cycle

        if (tmpqc(n) == 0) then
          iobsgrd(tmpi(n),tmpj(n),tmpk(n),tmploti(n),tmpelmi(n)) &
            = iobsgrd(tmpi(n),tmpj(n),tmpk(n),tmploti(n),tmpelmi(n)) + 1
          
          n1 = nnn + nobsgrdac(tmpi(n)-1,tmpj(n),tmpk(n),tmploti(n),tmpelmi(n)) &
                   + iobsgrd(tmpi(n),tmpj(n),tmpk(n),tmploti(n),tmpelmi(n))
          tmp2lon (n1) = tmplon (n)
          tmp2lat (n1) = tmplat (n)
          tmp2lev (n1) = tmplev (n)
          tmp2dat (n1) = tmpdat (n)
          tmp2err (n1) = tmperr (n)
          tmp2elmi(n1) = tmpelmi(n)
          tmp2typi(n1) = tmptypi(n)
          tmp2loti(n1) = tmploti(n)
          tmp2ri  (n1) = tmpri  (n)
          tmp2rj  (n1) = tmprj  (n)
          tmp2rk  (n1) = tmprk  (n)
          tmp2i   (n1) = tmpi   (n)
          tmp2j   (n1) = tmpj   (n)
          tmp2k   (n1) = tmpk   (n)
        end if
      end do

    !-----------------------------------------------------------------------

      nn = 0

      do ielm = 1, nid_obs

        if (obmethod_t(itype,ielm) > 0) then
          do k = 0, nlevx
            do j = 1, nlat
              do i = 1, nlon

    !-----------------------------------------------------------------------
                do is = 1, nslots
                  islot = slot_priority(is)

                  do m = 1, nobsgrd(i,j,k,islot,ielm)
                    n1 = nnn + nobsgrdac(i-1,j,k,islot,ielm) + m

        !---------------------------------------------------------------
                    if (obmethod_t(itype,ielm) == 2 .and. is > 1) then

                      tmp2qc(n1) = 2
                      nobsrdc_t = nobsrdc_t + 1
                      if (print_log) then
                        if (k == 0) then ! surface observations
                          write (84, '(A,A6,A,A3,A,I3,3(A,F7.2),A)') &
                            '[', obtypelist(itype), ':', obelmlist(ielm), &
                            '] slot=', islot, ' lon=', tmp2lon(n1), ' lat=', tmp2lat(n1), ' alt=', tmp2lev(n1), &
                            ' | obs removed because it is off the base time'
                        else
                          write (84, '(A,A6,A,A3,A,I3,3(A,F7.2),A)') &
                            '[', obtypelist(itype), ':', obelmlist(ielm), &
                            '] slot=', islot, ' lon=', tmp2lon(n1), ' lat=', tmp2lat(n1), ' lev=', tmp2lev(n1)/100., &
                            ' | obs removed because it is off the base time'
                        end if
                      end if

        !---------------------------------------------------------------
                    else ! [ obmethod_t(itype,ielm) == 2 .and. is > 1 ]

 !----------------------------------
 superob_time_search: do is2 = 1, is
 !----------------------------------
                        islot2 = slot_priority(is2)
                        mm = nobsgrd(i,j,k,islot2,ielm)
                        if (is2 == is) mm = m - 1

                        do m2 = 1, mm
                          n2 = nnn + nobsgrdac(i-1,j,k,islot2,ielm) + m2

                          if (tmp2lon(n1) == tmp2lon(n2) .and. tmp2lat(n1) == tmp2lat(n2) .and. &
                              tmp2lev(n1) == tmp2lev(n2)) then
                            tmp2qc(n1) = 1
                            nobsrdc_t = nobsrdc_t + 1
                            if (print_log) then
                              if (k == 0) then ! surface observations
                                write (84, '(A,A6,A,A3,A,I3,3(A,F7.2),A,I3)') &
                                  '[', obtypelist(itype), ':', obelmlist(ielm), &
                                  '] slot=', islot, ' lon=', tmp2lon(n1), ' lat=', tmp2lat(n1), ' alt=', tmp2lev(n1), &
                                  ' | obs removed because it also exists in slot', islot2
                              else
                                write (84, '(A,A6,A,A3,A,I3,3(A,F7.2),A,I3)') &
                                  '[', obtypelist(itype), ':', obelmlist(ielm), &
                                  '] slot=', islot, ' lon=', tmp2lon(n1), ' lat=', tmp2lat(n1), ' lev=', tmp2lev(n1)/100., &
                                  ' | obs removed because it also exists in slot', islot2
                              end if
                            end if
                            exit superob_time_search
                          end if
                        end do
                      !-------------------------
                      end do superob_time_search ! [ is2 = 1, is ]
                      !-------------------------

                    end if ! [ obmethod_t(itype,ielm) == 2 .and. is > 1 ]
        !---------------------------------------------------------------

                  end do ! [ m = 1, nobsgrd(i,j,k,islot,ielm) ]
                end do   ! [ is = 1, nslots ]
    !-----------------------------------------------------------------------

              end do ! [ i = 1, nlon ]
            end do   ! [ j = 1, nlat ]
          end do     ! [ k = 0, nlevx ]
        end if

      end do ! [ ielm = 1, nid_obs ]

      nnn = nnn + nobs_thistype

    end if ! [ nobs_thistype > 0 ]

  !-----
  end do ! [ itype = 1, nobtype ]
  !-----

! tmpqc:
!  0: pass qc
!  1: dropped because of temporal superob
!  2: dropped because it is off the base time
!-------------------------------------------------------------------------------

  nobs = nnn
  tmplon (1:nobs) = tmp2lon (1:nobs)
  tmplat (1:nobs) = tmp2lat (1:nobs)
  tmplev (1:nobs) = tmp2lev (1:nobs)
  tmpdat (1:nobs) = tmp2dat (1:nobs)
  tmperr (1:nobs) = tmp2err (1:nobs)
  tmpelmi(1:nobs) = tmp2elmi(1:nobs)
  tmptypi(1:nobs) = tmp2typi(1:nobs)
  tmploti(1:nobs) = tmp2loti(1:nobs)
  tmpri  (1:nobs) = tmp2ri  (1:nobs)
  tmprj  (1:nobs) = tmp2rj  (1:nobs)
  tmprk  (1:nobs) = tmp2rk  (1:nobs)
  tmpi   (1:nobs) = tmp2i   (1:nobs)
  tmpj   (1:nobs) = tmp2j   (1:nobs)
  tmpk   (1:nobs) = tmp2k   (1:nobs)
  tmpqc  (1:nobs) = tmp2qc  (1:nobs)

  write (*, '(A,I10)') 'observation number reduced by temporal superob:', nobsrdc_t
  write (*, *)
  write (*, '(A)') 'statistics after temporal superob:'

  call print_obsnum_cal_qc(nobs, tmpelmi(1:nobs), tmptypi(1:nobs), tmpqc(1:nobs))
  
!===============================================================================

  write (*, *)
  write (*, '(A)') '================================================'
  write (*, '(A)') ' all grid superob (horizontal + vertical)       '
  write (*, '(A)') '================================================'

  nnn_before = 0
  nnn_after = 0

  !--------------------
  do itype = 1, nobtype
  !--------------------

    write (*, '(2A)') 'processing ', obtypelist(itype)

    nobsgrd(:,:,:,:,:) = 0
    nobsgrdac(:,:,:,:,:) = 0
    iobsgrd(:,:,:,:,:) = 0

    do n = 1, nobs
      if (tmptypi(n) /= itype) cycle
      if (tmpqc(n) == 0) then
        nobsgrd(tmpi(n),tmpj(n),tmpk(n),tmploti(n),tmpelmi(n)) &
          = nobsgrd(tmpi(n),tmpj(n),tmpk(n),tmploti(n),tmpelmi(n)) + 1
      end if
    end do

    accu = 0
    do ielm = 1, nid_obs
      do islot = 1, nslots
        do k = 0, nlevx
          do j = 1, nlat
            nobsgrdac(0,j,k,islot,ielm) = accu
            do i = 1, nlon
              accu = accu + nobsgrd(i,j,k,islot,ielm)
              nobsgrdac(i,j,k,islot,ielm) = accu
            end do
          end do
        end do
      end do
    end do
    nobs_thistype = accu

    !-----------------------------------------------------------------------

    if (nobs_thistype > 0) then

      do n = 1, nobs
        if (tmptypi(n) /= itype) cycle

        if (tmpqc(n) == 0) then
          iobsgrd(tmpi(n),tmpj(n),tmpk(n),tmploti(n),tmpelmi(n)) &
            = iobsgrd(tmpi(n),tmpj(n),tmpk(n),tmploti(n),tmpelmi(n)) + 1
          
          n1 = nnn_before + nobsgrdac(tmpi(n)-1,tmpj(n),tmpk(n),tmploti(n),tmpelmi(n)) &
                          + iobsgrd(tmpi(n),tmpj(n),tmpk(n),tmploti(n),tmpelmi(n))
          tmp2lon (n1) = tmplon (n)
          tmp2lat (n1) = tmplat (n)
          tmp2lev (n1) = tmplev (n)
          tmp2dat (n1) = tmpdat (n)
          tmp2err (n1) = tmperr (n)
!          tmp2elmi(n1) = tmpelmi(n)
!          tmp2typi(n1) = tmptypi(n)
!          tmp2loti(n1) = tmploti(n)
          tmp2ri  (n1) = tmpri  (n)
          tmp2rj  (n1) = tmprj  (n)
          tmp2rk  (n1) = tmprk  (n)
          tmp2i   (n1) = tmpi   (n)
          tmp2j   (n1) = tmpj   (n)
          tmp2k   (n1) = tmpk   (n)
        end if
      end do

    !-----------------------------------------------------------------------

      nn = 0

      do ielm = 1, nid_obs
        do islot = 1, nslots

          if (obmethod_h(itype,ielm) > 0) then

            do k = 0, nlevx
              do j = 1, nlat
                do i = 1, nlon
                  m = nobsgrd(i,j,k,islot,ielm)
                  if (m > 1) then
                    n1 = nnn_before + nobsgrdac(i-1,j,k,islot,ielm) + 1
                    n2 = nnn_before + nobsgrdac(i,j,k,islot,ielm)
                    call superob_select(m, i, j, k, &
                                        tmp2ri(n1:n2) , tmp2rj(n1:n2) , tmp2rk(n1:n2) , &
                                        tmp2lon(n1:n2), tmp2lat(n1:n2), tmp2lev(n1:n2), &
                                        tmp2dat(n1:n2), tmp2err(n1:n2), obmethod_h(itype,ielm))
                  end if

                  if (m > 0) then
                    n1 = nnn_before + nobsgrdac(i-1,j,k,islot,ielm) + 1
                    n2 = nnn_before + nobsgrdac(i-1,j,k,islot,ielm) + m
                    n3 = nnn_after + nn + 1
                    n4 = nnn_after + nn + m
                    tmp2lon (n3:n4) = tmp2lon(n1:n2)
                    tmp2lat (n3:n4) = tmp2lat(n1:n2)
                    tmp2lev (n3:n4) = tmp2lev(n1:n2)
                    tmp2dat (n3:n4) = tmp2dat(n1:n2)
                    tmp2err (n3:n4) = tmp2err(n1:n2)
!                    tmp2elmi(n3:n4) = tmp2elmi(n1:n2)
!                    tmp2typi(n3:n4) = tmp2typi(n1:n2)
!                    tmp2loti(n3:n4) = tmp2loti(n1:n2)
                    tmp2elmi(n3:n4) = ielm
                    tmp2typi(n3:n4) = itype
                    tmp2loti(n3:n4) = islot
                    tmp2ri  (n3:n4) = tmp2ri (n1:n2)
                    tmp2rj  (n3:n4) = tmp2rj (n1:n2)
                    tmp2rk  (n3:n4) = tmp2rk (n1:n2)
!                    tmp2i   (n3:n4) = tmp2i  (n1:n2)
!                    tmp2j   (n3:n4) = tmp2j  (n1:n2)
!                    tmp2k   (n3:n4) = tmp2k  (n1:n2)
                    tmp2i   (n3:n4) = i
                    tmp2j   (n3:n4) = j
                    tmp2k   (n3:n4) = k
                  end if

                  if (print_log .and. m < nobsgrd(i,j,k,islot,ielm)) then
                    if (k == 0) then ! surface observations
                      write (85, '(A,A6,A,A3,A,I3,2(A,F7.2),12x,2(A,I4))') &
                        '[', obtypelist(itype), ':', obelmlist(ielm), &
                        '] slot=', islot, ' lon=', lon(i), ' lat=', lat(j), &
                        ' | obs # reduced from', nobsgrd(i,j,k,islot,ielm), ' ->', m
                    else
                      write (85, '(A,A6,A,A3,A,I3,3(A,F7.2),2(A,I4))') &
                        '[', obtypelist(itype), ':', obelmlist(ielm), &
                        '] slot=', islot, ' lon=', lon(i), ' lat=', lat(j), ' lev=', refpres(i,j,k)/100., &
                        ' | obs # reduced from', nobsgrd(i,j,k,islot,ielm), ' ->', m
                    end if
                  end if

                  nn = nn + m
                end do ! [ i = 1, nlon ]
              end do   ! [ j = 1, nlat ]
            end do     ! [ k = 0, nlevx ]

          else ! [ obmethod_h(itype,ielm) > 0 ]

            m = nobsgrdac(nlon,nlat,nlevx,islot,ielm) - nobsgrdac(0,1,0,islot,ielm)
            if (m > 0) then
              n1 = nnn_before + nobsgrdac(i-1,1,0,islot,ielm) + 1
              n2 = nnn_before + nobsgrdac(nlon,nlat,nlevx,islot,ielm)
              n3 = nnn_after + nn + 1
              n4 = nnn_after + nn + m
              tmp2lon (n3:n4) = tmp2lon(n1:n2)
              tmp2lat (n3:n4) = tmp2lat(n1:n2)
              tmp2lev (n3:n4) = tmp2lev(n1:n2)
              tmp2dat (n3:n4) = tmp2dat(n1:n2)
              tmp2err (n3:n4) = tmp2err(n1:n2)
!              tmp2elmi(n3:n4) = tmp2elmi(n1:n2)
!              tmp2typi(n3:n4) = tmp2typi(n1:n2)
!              tmp2loti(n3:n4) = tmp2loti(n1:n2)
              tmp2elmi(n3:n4) = ielm
              tmp2typi(n3:n4) = itype
              tmp2loti(n3:n4) = islot
              tmp2ri  (n3:n4) = tmp2ri (n1:n2)
              tmp2rj  (n3:n4) = tmp2rj (n1:n2)
              tmp2rk  (n3:n4) = tmp2rk (n1:n2)
              tmp2i   (n3:n4) = tmp2i  (n1:n2)
              tmp2j   (n3:n4) = tmp2j  (n1:n2)
              tmp2k   (n3:n4) = tmp2k  (n1:n2)
            end if
            nn = nn + m

          end if ! [ obmethod_h(itype,ielm) > 0 ]

        end do ! [ islot = 1, nslots ]
      end do   ! [ ielm = 1, nid_obs ]

      nnn_before = nnn_before + nobs_thistype
      nnn_after  = nnn_after  + nn
      nobsrdc_h = nobsrdc_h + nobs_thistype - nn

    end if ! [ nobs_thistype > 0 ]

  !-----
  end do ! [ itype = 1, nobtype ]
  !-----

  nobs = nnn_after

  write (*, '(A,I10)') 'observation number reduced by all grid superob:', nobsrdc_h
  write (*, *)
  write (*, '(A)') 'statistics after all grid superob:'

  call print_obsnum_cal(nobs, tmp2elmi(1:nobs), tmp2typi(1:nobs))

!===============================================================================

  write (*, *)
  do islot = 1, nslots
    write (supfile(4:5),'(I2.2)') islot
    open (92, file = supfile, form = 'unformatted', access = 'sequential')
    write (*, '(2A)') 'write obs file after superob: ', supfile

    do n = 1, nobs
      if (tmp2loti(n) == islot) then
        levout = tmp2lev(n)
        datout = tmp2dat(n)
        errout = tmp2err(n)
        select case (elem_uid(tmp2elmi(n)))
        case (id_u_obs)
          levout = levout * 0.01 ! Pa -> hPa
        case (id_v_obs)
          levout = levout * 0.01 ! Pa -> hPa
        case (id_t_obs)
          levout = levout * 0.01 ! Pa -> hPa
        case (id_tv_obs)
          levout = levout * 0.01 ! Pa -> hPa
        case (id_q_obs)
          levout = levout * 0.01 ! Pa -> hPa
        case (id_ps_obs)
          datout = datout * 0.01 ! Pa -> hPa
          errout = errout * 0.01 ! Pa -> hPa
        case (id_rh_obs)
          levout = levout * 0.01 ! Pa -> hPa
          datout = datout * 100.0 ! percent output
          errout = errout * 100.0 ! percent output
        case (id_tcmip_obs)
          datout = datout * 0.01 ! Pa -> hPa
          errout = errout * 0.01 ! Pa -> hPa
        end select

        write (92) real(elem_uid(tmp2elmi(n)), r_sngl), &
                   real(tmp2lon(n), r_sngl), &
                   real(tmp2lat(n), r_sngl), &
                   real(levout, r_sngl), &
                   real(datout, r_sngl), &
                   real(errout, r_sngl), &
                   real(tmp2typi(n), r_sngl)
      end if
    end do

    close (92)
  end do

!-------------------------------------------------------------------------------

  if (write_density_grid) then
    write (*, *)
    write (*, '(2A)') 'write observation density file after superob: ', obsdenfile_so
    call write_obs_den(nobs, tmp2i(1:nobs), tmp2j(1:nobs), tmp2k(1:nobs), &
                       tmp2elmi(1:nobs), tmp2typi(1:nobs), tmp2qc(1:nobs), obsdenfile_so, 2)
  end if ! [ write_density_grid ]

!===============================================================================

  deallocate ( tmp2lon  )
  deallocate ( tmp2lat  )
  deallocate ( tmp2lev  )
  deallocate ( tmp2dat  )
  deallocate ( tmp2err  )
  deallocate ( tmp2elmi )
  deallocate ( tmp2typi )
  deallocate ( tmp2loti )
  deallocate ( tmp2ri   )
  deallocate ( tmp2rj   )
  deallocate ( tmp2rk   )
  deallocate ( tmp2i    )
  deallocate ( tmp2j    )
  deallocate ( tmp2k    )
  deallocate ( tmp2qc   )

  deallocate ( tmpelm  )
  deallocate ( tmplon  )
  deallocate ( tmplat  )
  deallocate ( tmplev  )
  deallocate ( tmpdat  )
  deallocate ( tmperr  )
  deallocate ( tmptyp  )
  deallocate ( tmpelmi )
  deallocate ( tmptypi )
  deallocate ( tmploti )
  deallocate ( tmpri   )
  deallocate ( tmprj   )
  deallocate ( tmprk   )
  deallocate ( tmpi    )
  deallocate ( tmpj    )
  deallocate ( tmpk    )
  deallocate ( tmpqc   )

!-------------------------------------------------------------------------------

  close (81)
  close (82)
  close (83)
  close (84)
  close (85)

!-------------------------------------------------------------------------------
end program
!===============================================================================
