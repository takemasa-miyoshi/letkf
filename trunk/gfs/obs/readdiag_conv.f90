program readdiag_conv
!===============================================================================
!
! [PURPOSE:] Read GSI diag_conv file and write out LETKF obs2 format
!
! [HISTORY:]
!   05/14/2013 yang           created
!   07/25/2013 yang           add data usage flag in the wk
!   08/01/2013 yang           change w(7) as h(x). Edit the document
!   09/10/2013 Guo-Yuan Lien  modified
!
!===============================================================================
!
! [GSI diag_conv file format]
!
!  cdiag(ii)      station id
!  rdiag(1,ii)    observation type
!  rdiag(2,ii)    observation subtype
!  rdiag(3,ii)    observation latitude (degrees)
!  rdiag(4,ii)    observation longitude (degrees)
!  rdiag(5,ii)    station elevation (meters)
!  rdiag(6,ii)    observation pressure (hPa)
!  rdiag(7,ii)    observation height (meters)
!  rdiag(8,ii)    obs time (hours relative to analysis time)
!  rdiag(9,ii)    input prepbufr qc or event mark
!  rdiag(10,ii)   setup qc or event mark (currently qtflg only)
!  rdiag(11,ii)   read_prepbufr data usage flag
!  rdiag(12,ii)   analysis usage flag (1=use, -1=not used)
!  rdiag(13,ii)   nonlinear qc relative weight
!  rdiag(14,ii)   prepbufr inverse obs error
!  rdiag(15,ii)   read_prepbufr inverse obs error
!  rdiag(16,ii)   final inverse observation error
!  rdiag(17,ii)   observation value
!  rdiag(18,ii)   obs-ges used in analysis
!  rdiag(19,ii)   obs-ges w/o bias correction (future slot)
!
!===============================================================================

  use common
  use common_gfs
  use common_obs_gfs
  use obs_tools
  implicit none

  character(9) :: obsinfile='obsin.dat'    ! input GSI diag file
  character(10) :: obsoutfile='obsout.dat' ! output LETKF obs2 file
  integer, parameter :: luin = 21
  integer, parameter :: luout = 22
  integer, parameter :: luerr = 23

  real(r_sngl), allocatable :: rdiag(:,:)  
  character(8), allocatable :: cdiag(:)
  character(3) :: dtype
  integer :: nchar, nreal, ii, mype, idate
  integer :: nobs, nobsqc
  integer :: nobst(nid_obs,nobtype+1)
  integer :: nobstqc(nid_obs,nobtype+1)
  integer :: ios

  nobst = 0
  nobstqc = 0

  open (luin, file=obsinfile, status='old', form='unformatted', &
              access='sequential', iostat=ios)
  if (ios /= 0) then
    write (0, '(2A)') '[Error] Open input file: ', obsinfile
    stop
  endif
  open (luout, file=obsoutfile, status='unknown', form='unformatted', &
               access='sequential', convert='big_endian', iostat=ios)
  if (ios /= 0) then
    write (0, '(2A)') '[Error] Open output file: ', obsoutfile
    stop
  endif

  read (luin) idate
  do 
    read (luin, iostat=ios) dtype, nchar, nreal, ii, mype
    if (ios /= 0) exit 

    allocate ( cdiag(ii) )
    allocate ( rdiag(nreal,ii) )

    read(luin, iostat=ios) cdiag, rdiag
    call output_conv(nreal, ii, rdiag, dtype, luout, luerr, nobst, nobstqc)

    deallocate ( cdiag )
    deallocate ( rdiag )
  end do    ! end reading diag file 

  close (luin)
  close (luout)

  nobs = sum(nobst)
  write (*, *)
  write (*, '(A)') 'Observation count in the input diag_conv file:'
  call print_obsnum(nobs, nobst)

  nobsqc = sum(nobstqc)
  write (*, *)
  write (*, '(A)') 'Observation count passed the quality control:'
  call print_obsnum(nobsqc, nobstqc)

end program readdiag_conv

!===============================================================================

subroutine output_conv (nreal, ii, rdiag, dtype, luout, luerr, nobst, nobstqc)

  use common, only: r_sngl
  use common_obs_gfs
  implicit none

  integer, intent(in)      :: nreal, ii
  real(r_sngl), intent(in) :: rdiag(nreal,ii)
  character(3), intent(in) :: dtype
  integer, intent(in)      :: luout
  integer, intent(in)      :: luerr
  integer, intent(inout)   :: nobst(nid_obs,nobtype+1)
  integer, intent(inout)   :: nobstqc(nid_obs,nobtype+1)

  real(r_sngl), parameter :: min_inv_err = 1.0e-7
  real(r_sngl) :: wk(10)
  integer :: i, iqm, itypebufr, itype, iuid, iuidv, id_obs, iqc

  select case (dtype)
    case (' uv')
      id_obs = id_u_obs
      iuidv = uid_obs(id_v_obs)
    case ('  t')
      id_obs = id_t_obs
    case ('  q')
      id_obs = id_q_obs
    case (' ps')
      id_obs = id_ps_obs
    case default
      write (luerr, '(A)') '* Unsupported observation type'
      write (luerr, '(3A,I8,A)') &
            '*   dtype=', dtype, ', total:', ii, ' obs skiped'
      return
  end select
  iuid = uid_obs(id_obs)
  wk = 0.0

!-------------------------------------------------------------------------------
  do i = 1, ii
!-------------------------------------------------------------------------------
    iqc = 1

    wk(1) = id_obs
    wk(2) = rdiag(4,i)
    wk(3) = rdiag(3,i)
    if (id_obs == id_ps_obs) then
      wk(4) = rdiag(5,i)
    else
      wk(4) = rdiag(6,i)
    end if

    wk(5) = rdiag(17,i)
    if (rdiag(16,i) < min_inv_err) then
      wk(6) = 1.0e10
      iqc = 0
      write (luerr, '(A)') '* Observation error is too large'
      write (luerr, '(3A,F6.2,A,F6.2,A,F7.2,A,E8.3)') &
            '*   dtype=', dtype, ', lon=', wk(2), ', lat=', wk(3), &
            ', lev=', wk(4), ', obserr^-1=', rdiag(16,i)
    else
      wk(6) = 1.0 / rdiag(16,i)
    end if

! refer to http://www.emc.ncep.noaa.gov/mmb/data_processing/prepbufr.doc/table_2.htm
    itypebufr = nint(rdiag(1,i))
    select case (itypebufr)
      case (120,132,220:221,232)
        itype = 1  ! ADPUPA
      case (122,222)
        itype = 1  ! ADPUPA (not used)
        iqc = 0
        write (luerr, '(A)') '* Observation platform is not used'
        write (luerr, '(3A,I3)') '*   dtype=', dtype, ', report type=', itypebufr
      case (133,233)
        itype = 2  ! AIRCAR
      case (130:131,135,230:231,235)
        itype = 3  ! AIRCFT
      case (134,234)
        itype = 3  ! AIRCFT (not used)
        iqc = 0
        write (luerr, '(A)') '* Observation platform is not used'
        write (luerr, '(3A,I3)') '*   dtype=', dtype, ', report type=', itypebufr
      case (241:243,245:246,250:254,257:259)
        itype = 4  ! SATWND
      case (240,244,247:249,255:256)
        itype = 4  ! SATWND (not used)
        iqc = 0
        write (luerr, '(A)') '* Observation platform is not used'
        write (luerr, '(3A,I3)') '*   dtype=', dtype, ', report type=', itypebufr
      case (223,228:229)
        itype = 5  ! PROFLR
      case (227)
        itype = 5  ! PROFLR (not used)
        iqc = 0
        write (luerr, '(A)') '* Observation platform is not used'
        write (luerr, '(3A,I3)') '*   dtype=', dtype, ', report type=', itypebufr
      case (224)
        itype = 6  ! VADWND
      case (181,183,187,281,284,287)
        itype = 8  ! ADPSFC
      case (180,182,280,282)
        itype = 9  ! SFCSHP
      case (191)
        itype = 10 ! SFCBOG (not used)
        iqc = 0
        write (luerr, '(A)') '* Observation platform is not used'
        write (luerr, '(3A,I3)') '*   dtype=', dtype, ', report type=', itypebufr
      case (150,152,283)
        itype = 11 ! SPSSMI
      case (111,210)
        itype = 12 ! SYNDAT (not used)
        iqc = 0
        write (luerr, '(A)') '* Observation platform is not used'
        write (luerr, '(3A,I3)') '*   dtype=', dtype, ', report type=', itypebufr
      case (286)
        itype = 13 ! ERS1DA
      case (151,156:175)
        itype = 14 ! GOESND (not used)
        iqc = 0
        write (luerr, '(A)') '* Observation platform is not used'
        write (luerr, '(3A,I3)') '*   dtype=', dtype, ', report type=', itypebufr
      case (285)
        itype = 15 ! QKSWND
      case (188,288)
        itype = 16 ! MSONET (not used)
        iqc = 0
        write (luerr, '(A)') '* Observation platform is not used'
        write (luerr, '(3A,I3)') '*   dtype=', dtype, ', report type=', itypebufr
      case (153)
        itype = 17 ! GPSIPW
      case (126)
        itype = 18 ! RASSDA
      case (289)
        itype = 19 ! WDSATR
      case (290)
        itype = 20 ! ASCATW
      case default
        itype = nobtype + 1 ! others
    end select
    wk(7) = real(itype, r_sngl)

    wk(8) = rdiag(8,ii)
    wk(9) = wk(5) - rdiag(18,i)
    if (id_obs == id_ps_obs) then
      wk(9) = wk(9) * 100.0
    end if

    iqm = nint(rdiag(9,i))
    if (iqm < 0 .or. iqm > 2) then
      iqc = 0
      write (luerr, '(A)') '* PREPBUFR quality mark is bad'
      write (luerr, '(3A,F6.2,A,F6.2,A,F7.2,A,I3)') &
            '*   dtype=', dtype, ', lon=', wk(2), ', lat=', wk(3), &
            ', lev=', wk(4), ', mark=', iqm
    end if
    wk(10) = real(iqc, r_sngl)

    nobst(iuid,itype) = nobst(iuid,itype) + 1
    if (iqc == 1) then
      nobstqc(iuid,itype) = nobstqc(iuid,itype) + 1
    end if
    write (luout) wk

    ! V observation
    if (dtype == ' uv') then
      wk(1) = id_v_obs
      wk(5) = rdiag(20,i)
      wk(9) = wk(5) - rdiag(21,i)

      nobst(iuidv,itype) = nobst(iuidv,itype) + 1
      if (iqc == 1) then
        nobstqc(iuidv,itype) = nobstqc(iuidv,itype) + 1
      end if
      write (luout) wk
    end if
!-------------------------------------------------------------------------------
  end do ! [ i = 1, ii ]
!-------------------------------------------------------------------------------

end subroutine output_conv

!===============================================================================
