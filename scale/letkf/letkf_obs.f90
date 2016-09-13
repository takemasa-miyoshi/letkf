MODULE letkf_obs
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!   10/04/2012 Guo-Yuan Lien     modified for GFS model
!   08/30/2013 Guo-Yuan Lien     separating obs operator, following changes by Takemasa MIYOSHI
!   01/01/2014 Guo-Yuan Lien     add EFSO, following changes by Daisuke HOTTA
!
!=======================================================================
!$USE OMP_LIB
  USE common
  use common_nml
  USE common_mpi
  USE common_scale
  USE common_obs_scale
  USE common_mpi_scale
  USE common_letkf
!  USE common_scalelib
!  USE common_precip

  IMPLICIT NONE
  PUBLIC

!  real(r_size),parameter :: dist_zero_fac = SQRT(10.0d0/3.0d0) * 2.0d0            ! 3.651483717
!  real(r_size),parameter :: dist_zero_fac_square = dist_zero_fac * dist_zero_fac  ! 13.33333333
  real(r_size),parameter :: dist_zero_fac = 3.651483717
  real(r_size),parameter :: dist_zero_fac_square = 13.33333333

  real(r_size),save :: dlon_zero
  real(r_size),save :: dlat_zero

  type(obs_info),allocatable,save :: obs(:)
  type(obs_da_value),save :: obsda
  type(obs_da_value),allocatable,save :: obsda2(:)  ! sorted
                                                    !!!!!! need to add %err and %dat if they can be determined in letkf_obs.f90
  integer,save :: nobs_ext

  integer,allocatable,save :: nobsgrd(:,:,:)
  integer,allocatable,save :: nobsgrd2(:,:,:)
  integer,save :: nobstotalg
  integer,save :: nobstotal

CONTAINS
!-----------------------------------------------------------------------
! Initialize
!-----------------------------------------------------------------------
SUBROUTINE set_letkf_obs
  use scale_grid, only: &
    DX, &
    DY
  use scale_grid_index, only: &
    IHALO,JHALO
  use scale_process, only: &
!    MPI_COMM_d => LOCAL_COMM_WORLD, &
    PRC_myrank
  use scale_rm_process, only: &
    PRC_NUM_X, &
    PRC_NUM_Y


  IMPLICIT NONE
!  REAL(r_size),PARAMETER :: gross_error=10.0d0 !!!!! move to namelist
  INTEGER :: n,i,j,ierr,im,iof,iidx

  integer :: mem_ref
!  CHARACTER(8) :: cdffile_m='cdfm.grd'         ! GYL, PRECIP assimilation
!  CHARACTER(8) :: cdffile_o='cdfo.grd'         ! GYL
!  CHARACTER(10) :: maskfile='ppmask.grd'       ! GYL
!  REAL(r_size) :: ppcdf_m(nlon,nlat,0:ncdf)    ! GYL
!  REAL(r_size) :: ppcdf_o(nlon,nlat,0:ncdf)    ! GYL
!  REAL(r_size) :: ppzero_m(nlon,nlat)          ! GYL
!  REAL(r_size) :: ppzero_o(nlon,nlat)          ! GYL
!  REAL(r_size) :: ppmask(nlon,nlat)            ! GYL
!  INTEGER :: pp_mem, il, bg_lev, ob_lev        ! GYL
!  INTEGER :: pp_ntotal(pp_bg_nlev,pp_ob_nlev)  ! GYL
!  INTEGER :: zero_mem                          ! GYL
!  REAL(r_size) :: ym, sigma                    ! GYL
!  REAL(r_size) :: ri, rj                       ! GYL
!  REAL(r_size) :: tmpdat_o, obserr_p, obserr_n ! GYL
!  INTEGER :: ii, jj                            ! GYL

  integer :: it,ip
  REAL(r_size),allocatable :: bufr(:,:)
  INTEGER,allocatable :: bufri(:)
  INTEGER,allocatable :: bufri2(:,:,:)
#ifdef H08
  REAL(r_size),allocatable :: bufr2(:) ! H08
  REAL(r_size):: ch_num ! H08
!  REAL(r_size),allocatable :: hx_sprd(:) ! H08
#endif
  integer :: iproc,jproc
  integer,allocatable :: nnext(:,:)

!---
  integer :: ns
  integer, allocatable :: nr(:), nrt(:), obsidx(:)
  type(obs_da_value) :: obsbufs, obsbufr
  integer :: ip2, imin1,imax1,jmin1,jmax1,imin2,imax2,jmin2,jmax2
!---




  real(r_size),allocatable :: tmpelm(:)
  INTEGER :: monit_nobs(nid_obs)
  REAL(r_size) :: bias(nid_obs)
  REAL(r_size) :: rmse(nid_obs)

  character(filelenmax) :: obsdafile
  character(11) :: obsda_suffix = '.000000.dat'

  type(obs_da_value) :: obsda_ext


  WRITE(6,'(A)') 'Hello from set_letkf_obs'


  !!!!!! changes for different observation types.... (do not communicate all observaitons in the same way...)
  dlon_zero = max(SIGMA_OBS, SIGMA_OBS_RADAR, SIGMA_OBS_RADAR_OBSNOREF, SIGMA_OBS_TC) * dist_zero_fac / DX
  dlat_zero = max(SIGMA_OBS, SIGMA_OBS_RADAR, SIGMA_OBS_RADAR_OBSNOREF, SIGMA_OBS_TC) * dist_zero_fac / DY
!  dlon_zero = max(SIGMA_OBS, SIGMA_OBS_RADAR, SIGMA_OBS_RADAR_OBSNOREF, SIGMA_OBS_RAIN) * dist_zero_fac / DX
!  dlat_zero = max(SIGMA_OBS, SIGMA_OBS_RADAR, SIGMA_OBS_RADAR_OBSNOREF, SIGMA_OBS_RAIN) * dist_zero_fac / DY
#ifdef H08
  dlon_zero = max(dlon_zero,SIGMA_OBS_H08 * dist_zero_fac / DX) ! H08
  dlat_zero = max(dlat_zero,SIGMA_OBS_H08 * dist_zero_fac / DY) ! H08
#endif


! Read externally processed observations
!-----------------------------------

  do it = 1, nitmax
    im = proc2mem(1,it,myrank+1)
    if (im >= 1 .and. im <= MEMBER) then
      if (it == 1) then
        WRITE(6,'(A,I10)') 'Internally processed observations: ', obsda%nobs - nobs_ext
        WRITE(6,'(A,I10)') 'Externally processed observations: ', nobs_ext
        WRITE(6,'(A,I10)') 'Total                observations: ', obsda%nobs
      end if

      if (OBSDA_IN .and. nobs_ext > 0) then
        obsda_ext%nobs = nobs_ext
        call obs_da_value_allocate(obsda_ext,0)
        write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is reading externally processed observations for member ', &
              im, ', subdomain id #', proc2mem(2,it,myrank+1)
        call file_member_replace(im, OBSDA_IN_BASENAME, obsdafile)
        write (obsda_suffix(2:7),'(I6.6)') proc2mem(2,it,myrank+1)
        call read_obs_da(trim(obsdafile)//obsda_suffix,obsda_ext,0)

        if (OBSDA_OUT) then
          write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is appending observations for member ', &
                im, ', subdomain id #', proc2mem(2,it,myrank+1)
          call file_member_replace(im, OBSDA_OUT_BASENAME, obsdafile)
!          write (obsda_suffix(2:7),'(I6.6)') proc2mem(2,it,myrank+1)
          call write_obs_da(trim(obsdafile)//obsda_suffix,obsda_ext,0,append=.true.)
        end if

        ! variables without an ensemble dimension
        if (it == 1) then
          obsda%set(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%set
          obsda%idx(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%idx
          obsda%ri(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%ri
          obsda%rj(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%rj
          obsda%qc(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%qc
#ifdef H08
          obsda%lev(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%lev
          obsda%val2(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%val2
#endif
        else
          if (maxval(abs(obsda%set(obsda%nobs-nobs_ext+1:obsda%nobs) - obsda_ext%set)) > 0) then
            write (6,'(A)') 'error: obsda%set are inconsistent among the ensemble'
            stop
          end if
          if (maxval(abs(obsda%idx(obsda%nobs-nobs_ext+1:obsda%nobs) - obsda_ext%idx)) > 0) then
            write (6,'(A)') 'error: obsda%idx are inconsistent among the ensemble'
            stop
          end if
          if (maxval(abs(obsda%ri(obsda%nobs-nobs_ext+1:obsda%nobs) - obsda_ext%ri)) > 1.e-6) then
            write (6,'(A)') 'error: obsda%ri are inconsistent among the ensemble'
            stop
          end if
          if (maxval(abs(obsda%rj(obsda%nobs-nobs_ext+1:obsda%nobs) - obsda_ext%rj)) > 1.e-6) then
            write (6,'(A)') 'error: obsda%rj are inconsistent among the ensemble'
            stop
          end if
          obsda%qc(obsda%nobs-nobs_ext+1:obsda%nobs) = max(obsda%qc(obsda%nobs-nobs_ext+1:obsda%nobs), obsda_ext%qc)
#ifdef H08
          obsda%lev(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda%lev(obsda%nobs-nobs_ext+1:obsda%nobs) + obsda_ext%lev
          obsda%val2(obsda%nobs-nobs_ext+1:obsda%nobs) = obsda%val2(obsda%nobs-nobs_ext+1:obsda%nobs) + obsda_ext%val2
#endif
        end if

        ! variables with an ensemble dimension
        obsda%ensval(im,obsda%nobs-nobs_ext+1:obsda%nobs) = obsda_ext%val
      end if ! [ OBSDA_IN .and. nobs_ext > 0 ]
    end if ! [ im >= 1 .and. im <= MEMBER ]
  end do ! [ it = 1, nitmax ]

! All_reduce the observations
!-----------------------------------
!
! if the number of processors is greater then the ensemble size,
! broadcast the observation indices and real grid numbers
! from myrank_e=MEMBER-1 to the rest of processors that didn't read anything.
!
!####################### now broadcast to all, need to be corrected.
  if (nprocs_e > MEMBER) then

!      ALLOCATE(ranks(nprocs_e-MEMBER+1))
!      do n = MEMBER, nprocs_e
!        ranks(n-MEMBER+1) = n-1
!      end do
!      call MPI_Comm_group(MPI_COMM_e,MPI_G_e,ierr)
!      call MPI_GROUP_INCL(MPI_G_e,nprocs_e-MEMBER+1,ranks,MPI_G_obstmp,ierr)
!      call MPI_COMM_CREATE(MPI_COMM_e,MPI_G_obstmp,MPI_COMM_obstmp,ierr)

!      IF(myrank_e+1 >= MEMBER) THEN
!        CALL MPI_BARRIER(MPI_COMM_obstmp,ierr)
!        call MPI_BCAST(obsda%nobs, 1, MPI_INTEGER, 0, MPI_COMM_obstmp, ierr)
!!    print *, myrank, obsda%nobs

!        if (myrank_e+1 > MEMBER) then
!          CALL obs_da_value_allocate(obsda,MEMBER)
!        end if

!        CALL MPI_BARRIER(MPI_COMM_obstmp,ierr)
!        call MPI_BCAST(obsda%set, obsda%nobs, MPI_INTEGER, 0, MPI_COMM_obstmp, ierr)
!        call MPI_BCAST(obsda%idx, obsda%nobs, MPI_INTEGER, 0, MPI_COMM_obstmp, ierr)
!        call MPI_BCAST(obsda%ri, obsda%nobs, MPI_r_size, 0, MPI_COMM_obstmp, ierr)
!        call MPI_BCAST(obsda%rj, obsda%nobs, MPI_r_size, 0, MPI_COMM_obstmp, ierr)
!!        CALL MPI_BARRIER(MPI_COMM_obstmp,ierr)
!      end if

!      deallocate(ranks)
!      call MPI_Comm_free(MPI_COMM_obstmp,ierr)


      CALL MPI_BARRIER(MPI_COMM_e,ierr)
      call MPI_BCAST(obsda%nobs, 1, MPI_INTEGER, 0, MPI_COMM_e, ierr)
!    print *, myrank, obsda%nobs

      if (myrank_e+1 > MEMBER) then
        CALL obs_da_value_allocate(obsda,MEMBER)
      end if

      CALL MPI_BARRIER(MPI_COMM_e,ierr)
      call MPI_BCAST(obsda%set, obsda%nobs, MPI_INTEGER, 0, MPI_COMM_e, ierr)
      call MPI_BCAST(obsda%idx, obsda%nobs, MPI_INTEGER, 0, MPI_COMM_e, ierr)
      call MPI_BCAST(obsda%ri, obsda%nobs, MPI_r_size, 0, MPI_COMM_e, ierr)
      call MPI_BCAST(obsda%rj, obsda%nobs, MPI_r_size, 0, MPI_COMM_e, ierr)
!        CALL MPI_BARRIER(MPI_COMM_e,ierr)


  end if ! [ nprocs_e > MEMBER ]


! obsda%val not used; averaged from obsda%ensval later

  allocate (bufr(MEMBER,obsda%nobs))
  bufr = 0.0d0
  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  CALL MPI_ALLREDUCE(obsda%ensval,bufr,obsda%nobs*MEMBER,MPI_r_size,MPI_SUM,MPI_COMM_e,ierr)
  obsda%ensval = bufr
  deallocate(bufr)

  allocate (bufri(obsda%nobs))
  bufri = 0
  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  CALL MPI_ALLREDUCE(obsda%qc,bufri,obsda%nobs,MPI_INTEGER,MPI_MAX,MPI_COMM_e,ierr)
  obsda%qc = bufri
  deallocate(bufri)

#ifdef H08
!-- H08
! calculate the ensemble mean of obsda%lev
!
  allocate (bufr2(obsda%nobs))
  bufr2 = 0.0d0
  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  CALL MPI_ALLREDUCE(obsda%lev,bufr2,obsda%nobs,MPI_r_size,MPI_SUM,MPI_COMM_e,ierr)
  obsda%lev = bufr2
  deallocate(bufr2)

  obsda%lev = obsda%lev / REAL(MEMBER,r_size)

! calculate the ensemble mean of obsda%val2 (clear sky BT)
!
  allocate (bufr2(obsda%nobs))
  bufr2 = 0.0d0
  CALL MPI_BARRIER(MPI_COMM_e,ierr)
  CALL MPI_ALLREDUCE(obsda%val2,bufr2,obsda%nobs,MPI_r_size,MPI_SUM,MPI_COMM_e,ierr)
  obsda%val2 = bufr2
  deallocate(bufr2)

  obsda%val2 = obsda%val2 / REAL(MEMBER,r_size)

!-- H08
#endif


!    call MPI_Comm_free(MPI_COMM_e,ierr)


!!                                                                               ! GYL, PRECIP assimilation
!! reading precipitation transformation definition and mask                      ! GYL
!!                                                                               ! GYL
!  if (opt_pptrans >= 2) then                                                    ! GYL
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',cdffile_m          ! GYL
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',cdffile_o          ! GYL
!    call read_ppcdf(cdffile_m, cdffile_o, ppcdf_m, ppcdf_o, ppzero_m, ppzero_o) ! GYL
!  end if                                                                        ! GYL
!  WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',maskfile             ! GYL
!  call read_ppmask(maskfile, ppmask)                                            ! GYL
!  pp_ntotal = 0                                                                 ! GYL

!
! preprocess 'obs'
!-----------------------------------
  do iof = 1, OBS_IN_NUM
    do n = 1, obs(iof)%nobs
      if (obs(iof)%elm(n) == id_radar_ref_obs) then
        if (obs(iof)%dat(n) >= 0.0d0 .and. obs(iof)%dat(n) < 1.0d10) then
          if (obs(iof)%dat(n) < MIN_RADAR_REF) then
            obs(iof)%dat(n) = MIN_RADAR_REF_DBZ + LOW_REF_SHIFT
          else
            obs(iof)%dat(n) = 10.0d0 * log10(obs(iof)%dat(n))
          end if
        else
          obs(iof)%dat(n) = undef
        end if
        if (USE_OBSERR_RADAR_REF) then
          obs(iof)%err(n) = OBSERR_RADAR_REF
        end if
      end if

      if (USE_OBSERR_RADAR_VR .AND. obs(iof)%elm(n) == id_radar_vr_obs) then
        obs(iof)%err(n) = OBSERR_RADAR_VR
      end if
    end do ! [ n = 1, obs(iof)%nobs ]
  end do ! [ iof = 1, OBS_IN_NUM ]


! Compute perturbation and departure
! gross error check
!-----------------------------------

  allocate(tmpelm(obsda%nobs))

#ifdef H08
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i,iof,iidx,mem_ref,ch_num)
#else
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i,iof,iidx,mem_ref)
#endif
  do n = 1, obsda%nobs
    IF(obsda%qc(n) > 0) CYCLE

    iof = obsda%set(n)
    iidx = obsda%idx(n)

    tmpelm(n) = obs(iof)%elm(iidx)


!!!###### PRECIP assimilation ######
!!    if (tmpelm(n) == id_rain_obs) then

!!      CALL phys2ij(tmplon(n),tmplat(n),ri,rj)
!!      ii = CEILING(ri-0.5) ! nearest point
!!      jj = CEILING(rj-0.5) ! nearest point
!!      if (ii < 1)    ii = ii + nlon
!!      if (ii > nlon) ii = ii - nlon
!!      if (jj < 1)    jj = 1
!!      if (jj > nlat) jj = nlat

!!      if (ppmask(ii,jj) < mask_thres) then
!!        tmpqc(n) = 0
!!!        write (6,'(A)') '* Precipitation not used because of the mask file'
!!!        write (6,'(A,F6.2,A,F6.2,A)') &
!!!              '*  (lon,lat)=(',tmplon(n),',',tmplat(n),')'
!!!        cycle
!!      end if

!!      pp_mem = 0
!!      do i = 1, MEMBER
!!        if (tmphdxf(n,i) >= ppzero_thres) then
!!          pp_mem = pp_mem + 1
!!        end if
!!      end do

!!      bg_lev = pp_bg_nlev
!!      do il = 1, pp_bg_nlev-1
!!        if (pp_mem < pp_bg_levs(il)) then
!!          bg_lev = il
!!          exit
!!        end if
!!      end do
!!      ob_lev = pp_ob_nlev
!!      do il = 1, pp_ob_nlev-1
!!        if (tmpdat(n) < pp_ob_levs(il)) then
!!          ob_lev = il
!!          exit
!!        end if
!!      end do
!!      pp_ntotal(bg_lev,ob_lev) = pp_ntotal(bg_lev,ob_lev) + 1
!!      tmpqc(n) = 100 + pp_mem  !! For precip, qc = 100 + number of members with precip

!!      if (.not. pp_criterion(bg_lev,ob_lev)) then
!!        tmpqc(n) = 0
!!!        write (6,'(A)') '* Precipitation does not fit assimilation criterion'
!!!        write (6,'(A,F6.2,A,F6.2,A,I3,A,I2,A,F7.3,A,I2,A)') &
!!!              '*  (lon,lat)=(',tmplon(n),',',tmplat(n),'), pp_mem=', &
!!!              pp_mem, '(', bg_lev, '), pp_obs=', tmpdat(n), '(', ob_lev, ')'
!!        cycle
!!      end if

!!      if (opt_pptrans >= 1) then
!!        tmpdat_o = tmpdat(n)
!!        tmplev(n) = tmpdat_o  !! For precip, lev = original observed value if transformation is used
!!        if (opt_pptrans == 1) then ! log transformation
!!          do i = 1, MEMBER
!!            tmphdxf(n,i) = pptrans_log(tmphdxf(n,i))
!!          end do
!!          tmpdat(n) = pptrans_log(tmpdat_o)
!!        else if (opt_pptrans == 2) then ! Gaussian transformation with median zero rain
!!          do i = 1, MEMBER
!!            tmphdxf(n,i) = pptrans_normal(tmphdxf(n,i), ppcdf_m(ii,jj,:), ppzero_m(ii,jj))
!!          end do
!!          tmpdat(n) = pptrans_normal(tmpdat_o, ppcdf_o(ii,jj,:), ppzero_o(ii,jj))
!!        else if (opt_pptrans == 3) then ! Gaussian transformation with modified median zero rain
!!          call pptrans_normal_mdzero_def(tmphdxf(n,:), ppcdf_m(ii,jj,:), ppzero_m(ii,jj), zero_mem, ym, sigma)
!!          tmpdat(n) = pptrans_normal_mdzero(tmpdat_o, ppcdf_o(ii,jj,:), ppzero_o(ii,jj), ppzero_m(ii,jj), zero_mem, ym, sigma)
!!        end if

!!        if (opt_ppobserr == 1) then ! transformed obserr from obs data file
!!          if (opt_pptrans == 1) then ! log transformation
!!            tmperr(n) = tmperr(n) / (tmpdat_o + log_trans_tiny)
!!            if (tmperr(n) < min_ppobserr) tmperr(n) = min_ppobserr
!!          else if (opt_pptrans == 2) then ! Gaussian transformation with median zero rain
!!            obserr_p = pptrans_normal(tmpdat_o+tmperr(n), ppcdf_o(ii,jj,:), ppzero_o(ii,jj)) - tmpdat(n)
!!            if (obserr_p < min_ppobserr) obserr_p = min_ppobserr
!!            obserr_n = tmpdat(n) - pptrans_normal(tmpdat_o-tmperr(n), ppcdf_o(ii,jj,:), ppzero_o(ii,jj))
!!            if (obserr_n < min_ppobserr) obserr_n = min_ppobserr
!!            tmperr(n) = 0.5d0 * (obserr_p + obserr_n)
!!          else if (opt_pptrans == 3) then ! Gaussian transformation with modified median zero rain
!!            obserr_p = pptrans_normal_mdzero(tmpdat_o+tmperr(n), ppcdf_o(ii,jj,:), ppzero_o(ii,jj), ppzero_m(ii,jj), zero_mem, ym, sigma) - tmpdat(n)
!!            if (obserr_p < min_ppobserr) obserr_p = min_ppobserr
!!            obserr_n = tmpdat(n) - pptrans_normal_mdzero(tmpdat_o-tmperr(n), ppcdf_o(ii,jj,:), ppzero_o(ii,jj), ppzero_m(ii,jj), zero_mem, ym, sigma)
!!            if (obserr_n < min_ppobserr) obserr_n = min_ppobserr
!!            tmperr(n) = 0.5d0 * (obserr_p + obserr_n)
!!          end if
!!        else if (opt_ppobserr == 2) then ! constant obserr
!!          tmperr(n) = const_ppobserr
!!        end if
!!      end if ! [ opt_pptrans == 1 ]

!!    end if ! [ tmpelm(n) == id_rain_obs ]
!!!###### end PRECIP assimilation ######



!!!###### RADAR assimilation ######
    if (obs(iof)%elm(iidx) == id_radar_ref_obs) then
      if (.not. USE_RADAR_REF) then
        obsda%qc(n) = iqc_otype
        cycle
      end if

      if (obs(iof)%dat(iidx) == undef) then
        obsda%qc(n) = iqc_obs_bad
        cycle
      end if

      !!! obsda%ensval: already converted to dBZ
      mem_ref = 0
      do i = 1, MEMBER
        if (obsda%ensval(i,n) > RADAR_REF_THRES_DBZ+1.0d-6) then
          mem_ref = mem_ref + 1
        end if
      end do
      if (obs(iof)%dat(iidx) > RADAR_REF_THRES_DBZ+1.0d-6) then
        if (mem_ref < MIN_RADAR_REF_MEMBER_OBSREF) then
          obsda%qc(n) = iqc_ref_mem
!          write (6,'(A)') '* Reflectivity does not fit assimilation criterion'
!          write (6,'(A,F6.2,A,F6.2,A,I6,A,F7.3)') &
!                '*  (lon,lat)=(',obs(iof)%lon(iidx),',',obs(iof)%lat(iidx),'), mem_ref=', &
!                mem_ref,', ref_obs=', obs(iof)%dat(iidx)
          cycle
        end if
      else
        if (mem_ref < MIN_RADAR_REF_MEMBER) then
          obsda%qc(n) = iqc_ref_mem
!          write (6,'(A)') '* Reflectivity does not fit assimilation criterion'
!          write (6,'(A,F6.2,A,F6.2,A,I6,A,F7.3)') &
!                '*  (lon,lat)=(',obs(iof)%lon(iidx),',',obs(iof)%lat(iidx),'), mem_ref=', &
!                mem_ref,', ref_obs=', obs(iof)%dat(iidx)
          cycle
        end if
      end if
    end if

    if (obs(iof)%elm(iidx) == id_radar_vr_obs) then
      if (.not. USE_RADAR_VR) then
        obsda%qc(n) = iqc_otype
        cycle
      end if
    end if

!    if (obs(iof)%elm(iidx) == id_radar_prh_obs) then
!      if (.not. USE_RADAR_PSEUDO_RH) then
!        obsda%qc(n) = iqc_otype
!        cycle
!      end if
!    end if
!!!###### end RADAR assimilation ######



#ifdef H08
!!!###### Himawari-8 assimilation ###### ! H08
    if (obs(iof)%elm(iidx) == id_H08IR_obs) then
      if (obs(iof)%dat(iidx) == undef) then
        obsda%qc(n) = iqc_obs_bad
        cycle
      end if

! -- reject Himawari-8 obs sensitivie above H08_LIMIT_LEV (Pa) ! H08 --
      if (obsda%lev(n) < H08_LIMIT_LEV) then
        obsda%qc(n) = iqc_obs_bad
        cycle
      endif

!
! -- Counting how many members have cloud.
! -- Cloudy members should have negative values.
!
      mem_ref = 0
      do i = 1, MEMBER
        if (obsda%ensval(i,n) < 0.0d0) then
          mem_ref = mem_ref + 1
          obsda%ensval(i,n) = obsda%ensval(i,n) * (-1.0d0)
        end if
      end do

!
! -- reject Band #11(ch=5) & #12(ch=6) of Himawari-8 obs ! H08
! -- because these channels are sensitive to chemical tracers
! NOTE!!
!    channel num of Himawari-8 obs is stored in obs%lev (T.Honda 11/04/2015)
!      if ((int(obs(iof)%elm(iidx)) == 11) .or. &
!          (int(obs(iof)%lev(iidx)) == 12)) then
!        obsda%qc(n) = iqc_obs_bad
!        cycle
!      endif
    endif
!!!###### end Himawari-8 assimilation ###### ! H08
#endif



    obsda%val(n) = obsda%ensval(1,n)
    DO i=2,MEMBER
      obsda%val(n) = obsda%val(n) + obsda%ensval(i,n)
    END DO
    obsda%val(n) = obsda%val(n) / REAL(MEMBER,r_size)
#ifdef H08
! Compute CA (cloud effect average, Okamoto et al. 2014QJRMS)
! CA is stored in obsda%val2

    obsda%val2(n) = (abs(obsda%val(n) - obsda%val2(n)) & ! CM
                   + abs(obs(iof)%dat(iidx) - obsda%val2(n)) &! CO
                   &) * 0.5d0
#endif
    DO i=1,MEMBER
      obsda%ensval(i,n) = obsda%ensval(i,n) - obsda%val(n) ! Hdx
    END DO
    obsda%val(n) = obs(iof)%dat(iidx) - obsda%val(n) ! y-Hx

!   compute sprd in obs space ! H08

!    hx_sprd(n) = 0.0d0 !H08
!    DO i=1,MEMBER
!      hx_sprd(n) = hx_sprd(n) + obsda%ensval(i,n) * obsda%ensval(i,n)
!    ENDDO
!    hx_sprd(n) = dsqrt(hx_sprd(n) / REAL(MEMBER,r_size))

    select case (obs(iof)%elm(iidx)) !gross error
    case (id_rain_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_RAIN * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case (id_radar_ref_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_RADAR_REF * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case (id_radar_vr_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_RADAR_VR * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case (id_radar_prh_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_RADAR_PRH * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case (id_H08IR_obs)
      ! Adaptive QC depending on the sky condition in the background.
      ! !!Not finished yet!!
      ! 
      ! In config.nml.obsope,
      !  H08_CLDSKY_THRS  < 0.0: turn off ! all members are diagnosed as cloudy.
      !  H08_CLDSKY_THRS  > 0.0: turn on
      !
      IF(mem_ref < H08_MIN_CLD_MEMBER)THEN ! Clear sky
        IF(ABS(obsda%val(n)) > 1.0d0 * obs(iof)%err(iidx)) THEN
          obsda%qc(n) = iqc_gross_err
        END IF
      ELSE ! Cloudy sky
        IF(ABS(obsda%val(n)) > GROSS_ERROR_H08 * obs(iof)%err(iidx)) THEN
          obsda%qc(n) = iqc_gross_err
        END IF
      END IF

      IF(obs(iof)%dat(iidx) < H08_BT_MIN)THEN
        obsda%qc(n) = iqc_gross_err
      ENDIF

!      IF(ABS(obsda%val(n)) > GROSS_ERROR_H08 * obs(iof)%err(iidx)) THEN
!        obsda%qc(n) = iqc_gross_err
!      END IF
    case (id_tclon_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_TCX * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case (id_tclat_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_TCY * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case (id_tcmip_obs)
      IF(ABS(obsda%val(n)) > GROSS_ERROR_TCP * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    case default
      IF(ABS(obsda%val(n)) > GROSS_ERROR * obs(iof)%err(iidx)) THEN
        obsda%qc(n) = iqc_gross_err
      END IF
    end select

    IF(obs(iof)%elm(iidx) == id_H08IR_obs)THEN

#ifdef H08
!
! Derived H08 obs height (based on the weighting function output from RTTOV fwd
! model) is substituted into obs%lev.
! Band num. is substituded into obsda%lev. This will be used in monit_obs.
!
      ch_num = obs(iof)%lev(iidx)

      IF(DEPARTURE_STAT_H08)THEN
!
! For obs err correlation statistics based on Desroziers et al. (2005, QJRMS).
!
        write(6, '(a,2I6,2F8.2,4F12.4,2I6,F10.4)')"H08-O-B", &
             obs(iof)%elm(iidx), &
             nint(ch_num), & ! obsda%lev includes band num.
             obs(iof)%lon(iidx), &
             obs(iof)%lat(iidx), &
             obsda%val(n),& ! O-B
             obsda%lev(n), & ! sensitive height
             obs(iof)%dat(iidx), &
             obs(iof)%err(iidx), &
             obsda%qc(n),        &
             mem_ref,  &  ! # of cloudy member
             obsda%val2(n)
      ELSE
        write(6, '(2I6,2F8.2,4F12.4,I3)') &
             obs(iof)%elm(iidx), & ! id
             nint(ch_num), & ! band num
             obs(iof)%lon(iidx), &
             obs(iof)%lat(n), &
             obsda%lev(iidx), & ! sensitive height
             obs(iof)%dat(iidx), &
             obs(iof)%err(iidx), &
             obsda%val(n), &
             obsda%qc(n)
      ENDIF !  [.not. DEPARTURE_STAT_H08]
#endif
    ELSE
      write (6, '(2I6,2F8.2,4F12.4,I3)') obs(iof)%elm(iidx), &
                                         obs(iof)%typ(iidx), &
                                         obs(iof)%lon(iidx), &
                                         obs(iof)%lat(iidx), &
                                         obs(iof)%lev(iidx), &
                                         obs(iof)%dat(iidx), &
                                         obs(iof)%err(iidx), &
                                         obsda%val(n), &
                                         obsda%qc(n)
    ENDIF


  END DO ! [ n = 1, obsda%nobs ]
!$OMP END PARALLEL DO

!!
!! output departure statistics
!!
!-------------------------------------------------------------------------------

  WRITE(6,'(A)') 'OBSERVATIONAL DEPARTURE STATISTICS (IN THIS SUBDOMAIN):'

  CALL monit_dep(obsda%nobs,tmpelm,obsda%val,obsda%qc,monit_nobs,bias,rmse)
  CALL monit_print(monit_nobs,bias,rmse)
  deallocate(tmpelm)

!!
!! temporal observation localization
!!
!!  DO n=1,nobs
!!    tmperr(n) = tmperr(n) * exp(0.25d0 * (tmpdif(n) / SIGMA_OBST)**2)
!!  END DO
!!
!! SELECT OBS IN THE NODE
!!

!  nobs = nn
!  WRITE(6,'(I10,A,I3.3)') nobs,' OBSERVATIONS TO BE ASSIMILATED IN MYRANK ',myrank



! Sorting
!-----------------------------------

  allocate ( obsda2(0:MEM_NP-1) )

  obsda2(PRC_myrank)%nobs = 0
  allocate ( nobsgrd(0:nlon,1:nlat,0:MEM_NP-1) )
  nobsgrd = 0
  do n = 1, obsda%nobs
    if (obsda%qc(n) == iqc_good) then
      obsda2(PRC_myrank)%nobs = obsda2(PRC_myrank)%nobs + 1
!      call rank_1d_2d(PRC_myrank, iproc, jproc)
!      i = ceiling(obsda%ri(n)-0.5) - IHALO - iproc * nlon
!      j = ceiling(obsda%rj(n)-0.5) - JHALO - jproc * nlat
      call ij_g2l(PRC_myrank, &
                  ceiling(obsda%ri(n)-0.5) - IHALO, &
                  ceiling(obsda%rj(n)-0.5) - JHALO, &
                  i, j)

!if (i < 1 .or. i > nlon .or. j < 1 .or. j > nlat) then
!  call rank_1d_2d(PRC_myrank, iproc, jproc)
!  print *, '%%%', myrank, PRC_myrank, obsda%ri(n) - real(IHALO + iproc * nlon, r_size), obsda%rj(n) - real(JHALO + jproc * nlat, r_size)
!end if

      if (i < 1) i = 1       ! Assume the process assignment was correct,
      if (i > nlon) i = nlon ! so this correction is only to remedy the round-off problem.
      if (j < 1) j = 1       !
      if (j > nlat) j = nlat !

      nobsgrd(i,j,PRC_myrank) = nobsgrd(i,j,PRC_myrank) + 1
    end if
  end do

!  write (6,'(I10,A)') obsda2(PRC_myrank)%nobs,' OBSERVATIONS TO BE ASSIMILATED'

! regional count
  do j = 1, nlat
    if (j > 1) then
      nobsgrd(0,j,PRC_myrank) = nobsgrd(nlon,j-1,PRC_myrank)
    end if
    do i = 1, nlon
      nobsgrd(i,j,PRC_myrank) = nobsgrd(i-1,j,PRC_myrank) + nobsgrd(i,j,PRC_myrank)
    end do
  end do

!do i = 0, MEM_NP-1
!do j = 1, nlat
!write (6,'(31I4)') nobsgrd(:,j,i)
!end do
!write (6,*)
!end do

!write (6, *) 'XXXXXX'



!do i = 0, MEM_NP-1
!do j = 1, nlat
!write (6,'(31I4)') nobsgrd(:,j,i)
!end do
!write (6,*)
!end do

!! global count
!    do j = 1, nlat
!      if (j > 1) then
!        nobsgrd(0,j) = nobsgrd(nlon,j-1)
!      end if
!      do i = 1, nlon
!        nobsgrd(i,j) = nobsgrd(i-1,j) + nobsgrd(i,j)
!      end do
!    end do

  allocate ( nnext (nlon,nlat) )
  nnext(1:nlon,:) = nobsgrd(0:nlon-1,:,PRC_myrank) + 1
  call obs_da_value_allocate(obsda2(PRC_myrank),MEMBER)
  do n = 1, obsda%nobs
    if (obsda%qc(n) == iqc_good) then
!      call rank_1d_2d(PRC_myrank, iproc, jproc)
!      i = ceiling(obsda%ri(n)-0.5) - IHALO - iproc * nlon
!      j = ceiling(obsda%rj(n)-0.5) - JHALO - jproc * nlat
      call ij_g2l(PRC_myrank, &
                  ceiling(obsda%ri(n)-0.5) - IHALO, &
                  ceiling(obsda%rj(n)-0.5) - JHALO, &
                  i, j)


      if (i < 1) i = 1       ! Assume the process assignment was correct,
      if (i > nlon) i = nlon ! so this correction is only to remedy the round-off problem.
      if (j < 1) j = 1       !
      if (j > nlat) j = nlat !


      obsda2(PRC_myrank)%set(nnext(i,j)) = obsda%set(n)
      obsda2(PRC_myrank)%idx(nnext(i,j)) = obsda%idx(n)
      obsda2(PRC_myrank)%val(nnext(i,j)) = obsda%val(n)
      obsda2(PRC_myrank)%ensval(:,nnext(i,j)) = obsda%ensval(:,n)
      obsda2(PRC_myrank)%qc(nnext(i,j)) = obsda%qc(n)
      obsda2(PRC_myrank)%ri(nnext(i,j)) = obsda%ri(n)
      obsda2(PRC_myrank)%rj(nnext(i,j)) = obsda%rj(n)
#ifdef H08
      obsda2(PRC_myrank)%lev(nnext(i,j)) = obsda%lev(n) ! H08
#endif

      nnext(i,j) = nnext(i,j) + 1
    end if
  end do
  deallocate (nnext)

  call obs_da_value_deallocate(obsda)



! Communication
!-----------------------------------

!#ifdef H08
! -- H08
!  if((obs(3)%nobs >= 1) .and. OBS_IN_NUM >= 3)then
!    allocate (bufr2(obs(3)%nobs))
!    bufr2 = 0.0d0
!    CALL MPI_BARRIER(MPI_COMM_d,ierr)
!    CALL MPI_ALLREDUCE(obs(3)%lev,bufr2,obs(3)%nobs,MPI_r_size,MPI_MAX,MPI_COMM_d,ierr)
!    obs(3)%lev = bufr2
!    deallocate(bufr2)
!  endif
! -- H08
!#endif

  allocate ( bufri2 (0:nlon,1:nlat,0:MEM_NP-1) )
  call MPI_ALLREDUCE(nobsgrd,bufri2,(nlon+1)*nlat*MEM_NP,MPI_INTEGER,MPI_SUM,MPI_COMM_d,ierr)
  nobsgrd(0:nlon,1:nlat,0:MEM_NP-1) = bufri2(0:nlon,1:nlat,0:MEM_NP-1)
  deallocate ( bufri2 )

  nobstotalg = sum(nobsgrd(nlon,nlat,:))
!  WRITE(6,'(A,I8)') 'Target observation numbers (global) : NOBS=',nobstotalg

  allocate ( nobsgrd2(0:nlon,1:nlat,0:MEM_NP-1) )
  nobsgrd2(:,:,PRC_myrank) = nobsgrd(:,:,PRC_myrank)

  allocate(nr(MEM_NP),nrt(MEM_NP))


!if (PRC_myrank == 0 .and. myrank_e == 0) then
!print *, nobsgrd(nlon,nlat,:)
!print *, maxval(nobsgrd(nlon,nlat,:))
!end if

  obsbufs%nobs = maxval(nobsgrd(nlon,nlat,:))
  call obs_da_value_allocate(obsbufs,MEMBER)

  allocate ( obsidx(maxval(nobsgrd(nlon,nlat,:))) )

  do ip = 0, MEM_NP-1
!    do ip = 0, 0

    call rank_1d_2d(ip, iproc, jproc)

!if (PRC_myrank == 0 .and. myrank_e == 0) then
!print *, PRC_NUM_X,nlon,iproc,jproc,SIGMA_OBS,DX,DY,dlon_zero,dlat_zero
!end if

    imin1 = max(1, iproc*nlon+1 - ceiling(dlon_zero))
    imax1 = min(PRC_NUM_X*nlon, (iproc+1)*nlon + ceiling(dlon_zero))
    jmin1 = max(1, jproc*nlat+1 - ceiling(dlat_zero))
    jmax1 = min(PRC_NUM_Y*nlat, (jproc+1)*nlat + ceiling(dlat_zero))

!if (myrank_e == 0) then
!print *, PRC_myrank,imin1,imax1,jmin1,jmax1
!end if

    ns = 0
    nr = 0
    nrt = 0

    do ip2 = 0, MEM_NP-1

      if (ip2 /= ip) then
        call rank_1d_2d(ip2, iproc, jproc)

        imin2 = max(1, imin1 - iproc*nlon)
        imax2 = min(nlon, imax1 - iproc*nlon)
        jmin2 = max(1, jmin1 - jproc*nlat)
        jmax2 = min(nlat, jmax1 - jproc*nlat)

!if (myrank_e == 0) then
!print *, PRC_myrank,ip2,imin2,imax2,jmin2,jmax2
!end if

        nr(ip2+1) = 0

        if (ip == PRC_myrank) then
          call obs_choose(imin2,imax2,jmin2,jmax2,ip2,nr(ip2+1),obsidx,.true.)
          obsda2(ip2)%nobs = nr(ip2+1)
          call obs_da_value_allocate(obsda2(ip2),MEMBER)
        else if (ip2 == PRC_myrank) then
          call obs_choose(imin2,imax2,jmin2,jmax2,ip2,nr(ip2+1),obsidx)
          ns = nr(ip2+1)
          do n = 1, ns
            obsbufs%set(n) = obsda2(PRC_myrank)%set(obsidx(n))
            obsbufs%idx(n) = obsda2(PRC_myrank)%idx(obsidx(n))
            obsbufs%val(n) = obsda2(PRC_myrank)%val(obsidx(n))
            obsbufs%ensval(:,n) = obsda2(PRC_myrank)%ensval(:,obsidx(n))
            obsbufs%qc(n) = obsda2(PRC_myrank)%qc(obsidx(n))
            obsbufs%ri(n) = obsda2(PRC_myrank)%ri(obsidx(n))
            obsbufs%rj(n) = obsda2(PRC_myrank)%rj(obsidx(n))
#ifdef H08
            obsbufs%lev(n) = obsda2(PRC_myrank)%lev(obsidx(n)) ! H08
#endif
          end do
        else
          call obs_choose(imin2,imax2,jmin2,jmax2,ip2,nr(ip2+1))
        end if


      end if

      if (ip2 > 0) then
        nrt(ip2+1) = nrt(ip2) + nr(ip2)
      end if

    end do ! [ ip2 = 0, MEM_NP-1 ]

    obsbufr%nobs = nrt(MEM_NP)+nr(MEM_NP)
    call obs_da_value_allocate(obsbufr,MEMBER)

    call MPI_GATHERV(obsbufs%set, ns, MPI_INTEGER, obsbufr%set, nr, nrt, MPI_INTEGER, ip, MPI_COMM_d, ierr)
    call MPI_GATHERV(obsbufs%idx, ns, MPI_INTEGER, obsbufr%idx, nr, nrt, MPI_INTEGER, ip, MPI_COMM_d, ierr)
    call MPI_GATHERV(obsbufs%val, ns, MPI_r_size, obsbufr%val, nr, nrt, MPI_r_size, ip, MPI_COMM_d, ierr)
    call MPI_GATHERV(obsbufs%ensval, ns*MEMBER, MPI_r_size, obsbufr%ensval, nr*MEMBER, nrt*MEMBER, MPI_r_size, ip, MPI_COMM_d, ierr)
    call MPI_GATHERV(obsbufs%qc, ns, MPI_INTEGER, obsbufr%qc, nr, nrt, MPI_INTEGER, ip, MPI_COMM_d, ierr)
    call MPI_GATHERV(obsbufs%ri, ns, MPI_r_size, obsbufr%ri, nr, nrt, MPI_r_size, ip, MPI_COMM_d, ierr)
    call MPI_GATHERV(obsbufs%rj, ns, MPI_r_size, obsbufr%rj, nr, nrt, MPI_r_size, ip, MPI_COMM_d, ierr)
#ifdef H08
    call MPI_GATHERV(obsbufs%lev, ns, MPI_r_size, obsbufr%lev, nr, nrt, MPI_r_size, ip, MPI_COMM_d, ierr) ! H08
#endif


    if (PRC_myrank == ip) then

      do ip2 = 0, MEM_NP-1
        if (ip2 /= ip) then
          obsda2(ip2)%set = obsbufr%set(nrt(ip2+1)+1:nrt(ip2+1)+nr(ip2+1))
          obsda2(ip2)%idx = obsbufr%idx(nrt(ip2+1)+1:nrt(ip2+1)+nr(ip2+1))
          obsda2(ip2)%val = obsbufr%val(nrt(ip2+1)+1:nrt(ip2+1)+nr(ip2+1))
          obsda2(ip2)%ensval = obsbufr%ensval(:,nrt(ip2+1)+1:nrt(ip2+1)+nr(ip2+1))
          obsda2(ip2)%qc = obsbufr%qc(nrt(ip2+1)+1:nrt(ip2+1)+nr(ip2+1))
          obsda2(ip2)%ri = obsbufr%ri(nrt(ip2+1)+1:nrt(ip2+1)+nr(ip2+1))
          obsda2(ip2)%rj = obsbufr%rj(nrt(ip2+1)+1:nrt(ip2+1)+nr(ip2+1))
#ifdef H08
          obsda2(ip2)%lev = obsbufr%lev(nrt(ip2+1)+1:nrt(ip2+1)+nr(ip2+1)) ! H08
#endif

!            write(6,*) obsda2(ip2)%idx

        end if
      end do

!      write(6,*) ns
!      write(6,*) nr(:)
!      write(6,*) nrt(:)
!!        write(6,*)
!!        write(6,*) nrt(MEM_NP)+nr(MEM_NP)
!!        write(6,*)
!!        write(6,*) obsbufr(:)

    end if

    call obs_da_value_deallocate(obsbufr)

  end do ! [ ip = 0, MEM_MP ]


  call obs_da_value_deallocate(obsbufs)
  deallocate(obsidx)

  nobsgrd = nobsgrd2

  deallocate (nobsgrd2)

  nobstotal = 0
  do ip = 0, MEM_NP-1
    nobstotal = nobstotal + obsda2(ip)%nobs


!write(6,'(A,I4,A,I10)') 'nobs(', ip+1, ') =', obsda2(ip)%nobs

  end do


!write(6,'(A,I10)') 'nobstotal =', nobstotal


  if (nobstotal /= sum(nobsgrd(nlon,nlat,:))) then

print *, myrank, nobstotalg, nobstotal, nobsgrd(nlon,nlat,:)

    print *, 'FFFFFF wrong!!'
    stop
  end if



!do i = 0, MEM_NP-1
!do j = 1, nlat
!write (6,'(31I4)') nobsgrd(:,j,i)
!end do
!write (6,*)
!end do


!do i = 0, MEM_NP-1
!  write(6,*) obsda2(i)%val(20:23)
!  write(6,*) obsda2(i)%ensval(:,20:23)
!end do

!    call unset_scalelib


  RETURN
END SUBROUTINE set_letkf_obs




SUBROUTINE obs_choose(imin,imax,jmin,jmax,proc,nn,nobs_use,nobsgrdout)

!  use scale_process, only: &
!    PRC_myrank


  implicit none

  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax,proc
  INTEGER,INTENT(INOUT) :: nn
  INTEGER,INTENT(INOUT),OPTIONAL :: nobs_use(:)
  logical,intent(in),optional :: nobsgrdout
  logical :: nobsgrdout_
  INTEGER :: j,ip

  nobsgrdout_ = .false.
  if (present(nobsgrdout)) nobsgrdout_ = nobsgrdout
  if (nobsgrdout_) then
    nobsgrd2(:,:,proc) = 0
  end if

  if (nobsgrd(nlon,nlat,proc) == 0) return

!  if (imin < 1 .or. imax > nlon) ....
  if (imin > imax .or. jmin > jmax) then
    return
  end if


  DO j = jmin, jmax

!if (myrank_e == 0 .and. PRC_myrank == 0) then
!print *, nobsgrd(imin-1,j,proc)+1, nobsgrd(imax,j,proc)
!end if


    DO ip = nobsgrd(imin-1,j,proc)+1, nobsgrd(imax,j,proc)
!      IF(nn > nobs) THEN
!        WRITE(6,*) 'FATALERROR, NN > NOBS', NN, NOBS
!      END IF
      nn = nn + 1
      if (present(nobs_use)) then
        nobs_use(nn) = ip
      end if
    END DO

    if (nobsgrdout_) then
      if (j > 1) then
        nobsgrd2(0:imax-1,j,proc) = nobsgrd2(nlon,j-1,proc)
      end if
      nobsgrd2(imin:imax,j,proc) = nobsgrd(imin:imax,j,proc) - nobsgrd(imin-1,j,proc) + nobsgrd2(imin-1,j,proc)
      if (imax < nlon) then
        nobsgrd2(imax+1:nlon,j,proc) = nobsgrd2(imax,j,proc)
      end if
    end if

  END DO

  if (nobsgrdout_) then
    if (jmax < nlat) then
      nobsgrd2(:,jmax+1:nlat,proc) = nobsgrd2(nlon,jmax,proc)
    end if
  end if

  RETURN
END SUBROUTINE obs_choose


!!-----------------------------------------------------------------------
!! Monitor observation diagnostics after the letkf update
!!  Moved from main program code to this subroutine, 2013/12/26 Guo-Yuan Lien
!!-----------------------------------------------------------------------
!SUBROUTINE monit_obs
!  IMPLICIT NONE
!  REAL(r_size),ALLOCATABLE :: ohx(:)
!  REAL(r_size),ALLOCATABLE :: dep(:)
!  INTEGER,ALLOCATABLE :: oqc(:)
!  INTEGER :: l,im
!  CHARACTER(10) :: ombfile='omb.dat'
!  CHARACTER(10) :: omafile='oma.dat'
!  CHARACTER(14) :: obsguesfile='obsguesNNN.dat'
!  CHARACTER(14) :: obsanalfile='obsanalNNN.dat'
!
!  IF(omb_output .AND. myrank == 0) THEN
!    ALLOCATE(ohx(nobs),oqc(nobs),dep(nobs))
!    CALL monit_output('gues',0,ohx,oqc)
!    dep = obsdat - ohx
!    WRITE(6,'(A)') 'OBSERVATIONAL DEPARTURE STATISTICS [gues mean]:'
!    CALL monit_dep(nobs,obselm,dep,oqc)
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',ombfile
!    CALL write_obs2(ombfile,nobs,obselm,obslon,obslat,obslev, &
!                    obsdat,obserr,obstyp,obsdif,dep,obsqc,0)
!  END IF
!  IF(oma_output .AND. myrank == 0) THEN
!    IF(.NOT. ALLOCATED(ohx)) ALLOCATE(ohx(nobs),oqc(nobs),dep(nobs))
!    CALL monit_output('anal',0,ohx,oqc)
!    dep = obsdat - ohx
!    WRITE(6,'(A)') 'OBSERVATIONAL DEPARTURE STATISTICS [anal mean]:'
!    CALL monit_dep(nobs,obselm,dep,oqc)
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',omafile
!    CALL write_obs2(omafile,nobs,obselm,obslon,obslat,obslev, &
!                    obsdat,obserr,obstyp,obsdif,dep,obsqc,0)
!  END IF

!  IF(obsgues_output) THEN
!    IF(.NOT. ALLOCATED(ohx)) ALLOCATE(ohx(nobs),oqc(nobs))
!    IF(myrank == 0) THEN
!      ohx = obsdat - obsdep
!      oqc = 1
!      WRITE(obsguesfile(8:10),'(A3)') '_me'
!      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsguesfile
!      CALL write_obs2(obsguesfile,nobs,obselm,obslon,obslat,obslev, &
!                      obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
!    END IF
!    l=0
!    DO
!      im = myrank+1 + nprocs * l
!      IF(im > MEMBER) EXIT
!      ohx(:) = obsdat(:) - obsdep(:) + obshdxf(:,im)
!      WRITE(obsguesfile(8:10),'(I3.3)') im
!      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsguesfile
!      CALL write_obs2(obsguesfile,nobs,obselm,obslon,obslat,obslev, &
!                      obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
!      l = l+1
!    END DO
!  END IF

!! This is not an accurate estimate of obsanal.
!! To obtain an accurate estimate of obsanal, use code in [letkf_tools:das_letkf]
!!------
!!  IF(obsanal_output) THEN
!!    IF(.NOT. ALLOCATED(ohx)) ALLOCATE(ohx(nobs),oqc(nobs))
!!    CALL monit_output('anal',0,ohx,oqc)
!!    WRITE(obsanalfile(8:10),'(A3)') '_me'
!!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsanalfile
!!    CALL write_obs2(obsanalfile,nobs,obselm,obslon,obslat,obslev, &
!!                    obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
!!    l=0
!!    DO
!!      im = myrank+1 + nprocs * l
!!      IF(im > MEMBER) EXIT
!!      CALL monit_output('anal',im,ohx,oqc)
!!      WRITE(obsanalfile(8:10),'(I3.3)') im
!!      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsanalfile
!!      CALL write_obs2(obsanalfile,nobs,obselm,obslon,obslat,obslev, &
!!                      obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
!!      l = l+1
!!    END DO
!!  END IF

!  IF(ALLOCATED(ohx)) DEALLOCATE(ohx,oqc)
!  IF(ALLOCATED(dep)) DEALLOCATE(dep)

!  RETURN
!END SUBROUTINE monit_obs
!-----------------------------------------------------------------------
! Monitor h(xb) or h(xa) from a LETKF output file
! Adopted from 'monit_mean' subroutine, 2013/12/24 Guo-Yuan Lien
!-----------------------------------------------------------------------
!  file: 'gues' or 'anal'
!  im:   member # (integer); 0 for ensmean (always called from myrank == 0)
!  ohx:  h(xb) or h(xa)
!  oqc:  quality of h(xb) or h(xa)
!-----------------------------------------------------------------------
!SUBROUTINE monit_output(file,im,ohx,oqc)
!  IMPLICIT NONE
!  CHARACTER(4),INTENT(IN) :: file
!  INTEGER,INTENT(IN) :: im
!  REAL(r_size),INTENT(OUT) :: ohx(nobs)
!  INTEGER,INTENT(OUT) :: oqc(nobs)
!  REAL(r_size) :: v3d(nlon,nlat,nlev,nv3dx)
!  REAL(r_size) :: v2d(nlon,nlat,nv2dx)
!  REAL(r_size) :: v3dtmp(nlon,nlat,nlev,nv3d)
!  REAL(r_size) :: v2dtmp(nlon,nlat,nv2d)
!  REAL(r_size) :: elem
!  REAL(r_size) :: bias_u,bias_v,bias_t,bias_ps,bias_q,bias_rh,bias_rain
!  REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_ps,rmse_q,rmse_rh,rmse_rain
!  REAL(r_size) :: hdxf,dep,ri,rj,rk
!  INTEGER :: n,iu,iv,it,iq,ips,irh,irain
!  CHARACTER(11) :: filename='filexxx.grd'
!  INTEGER :: iret
!  REAL(r_size) :: tmpps(nlon*nlat)
!  REAL(r_size) :: tmptv(nlon*nlat,nlev)
!  REAL(r_size) :: tmpp(nlon*nlat,nlev)



!  CHARACTER(4),INTENT(IN) :: file
!  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,MEMBER,nv3d)
!  REAL(r_size),INTENT(OUT) :: v2d(nij1,MEMBER,nv2d)
!  REAL(RP) :: v3dg(nlev,nlon,nlat,nv3d)
!  REAL(RP) :: v2dg(nlon,nlat,nv2d)
!  CHARACTER(9) :: filename='file.0000'
!  integer :: it,im,mstart,mend




!  IF(im == 0) THEN
!    WRITE(filename(1:7),'(A4,A3)') file,'_me'
!  ELSE
!    WRITE(filename(1:7),'(A4,I3.3)') file,im
!  END IF
!  CALL read_grd(filename,v3dtmp,v2dtmp,0) ! read ensemble mean into a temporary array
!  IF(im == 0) THEN
!    WRITE(filename(1:7),'(A4,I3.3)') 'gues',1 ! assume called from myrank == 0
!  ELSE
!    WRITE(filename(1:7),'(A4,I3.3)') 'gues',im
!  END IF
!  CALL read_grdx(filename,v3d,v2d) ! only the orography is used, P will be recalulated






!!-----------------------------------------------------------------------




!      WRITE(filename(1:4),'(A4)') file
!      WRITE(filename(6:9),'(I4.4)') im
!!      WRITE(6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',proc2mem(2,it,myrank+1),'.nc'
!      call read_restart(filename,v3dg,v2dg)






!!  v3d(:,:,:,iv3d_u) = v3dtmp(:,:,:,iv3d_u)
!!  v3d(:,:,:,iv3d_v) = v3dtmp(:,:,:,iv3d_v)
!!  v3d(:,:,:,iv3d_t) = v3dtmp(:,:,:,iv3d_t)
!!  v3d(:,:,:,iv3d_q) = v3dtmp(:,:,:,iv3d_q)
!!  v3d(:,:,:,iv3d_qc) = v3dtmp(:,:,:,iv3d_qc)
!!  v2d(:,:,iv2d_ps) = v2dtmp(:,:,iv2d_ps)
!!  tmpps = reshape(v2d(:,:,iv2d_ps),(/nlon*nlat/))
!!  tmptv = reshape(v3d(:,:,:,iv3d_t) * (1.0d0 + fvirt * v3d(:,:,:,iv3d_q)),(/nlon*nlat,nlev/))
!!  call sigio_modprd(nlon*nlat,nlon*nlat,nlev,gfs_nvcoord,gfs_idvc,gfs_idsl, &
!!                    gfs_vcoord,iret,tmpps,tmptv,pm=tmpp)
!!  v3d(:,:,:,iv3d_p) = reshape(tmpp,(/nlon,nlat,nlev/))

!!  oqc = 1
!!  DO n=1,nobs
!!    CALL phys2ijk(v3d(:,:,:,iv3d_p),obselm(n),obslon(n),obslat(n),obslev(n),ri,rj,rk)
!!    !
!!    ! For monitoring, don't skip any observation below or above model vertical extent.
!!    ! Just put bad QC but still get estimate.
!!    !
!!    IF(CEILING(rk) > nlev) THEN
!!      rk = REAL(nlev,r_size)
!!      oqc(n) = 0
!!    END IF
!!    IF(CEILING(rk) < 2 .AND. NINT(obselm(n)) /= id_ps_obs) THEN
!!      IF(NINT(obselm(n)) > 9999) THEN
!!        rk = 0.0d0
!!      ELSE IF(NINT(obselm(n)) == id_u_obs .OR. NINT(obselm(n)) == id_v_obs) THEN
!!        rk = 1.00001d0
!!      ELSE
!!        rk = 1.00001d0
!!        oqc(n) = 0
!!      END IF
!!    END IF
!!    IF(NINT(obselm(n)) == id_ps_obs) THEN
!!      CALL itpl_2d(v2d(:,:,iv2d_orog),ri,rj,rk)
!!      rk = obslev(n) - rk
!!    END IF
!!    IF(NINT(obselm(n)) == id_rain_obs) THEN ! No way to get the accumulated precipitation value
!!      ohx(n) = obsdat(n)
!!      oqc(n) = 0
!!    ELSE
!!      CALL Trans_XtoY(obselm(n),ri,rj,rk,v3d,v2d,ohx(n))
!!    END IF
!!  END DO

!  RETURN
!END SUBROUTINE monit_output
!!-----------------------------------------------------------------------
!! Read observation diagnostics for EFSO
!!  Adapted from Y.Ohta's EFSO code for SPEEDY-LETKF   2013/07/17 D.Hotta
!!  Modified, renamed from 'read_monit_obs' to 'set_efso_obs', 2013/12/26 Guo-Yuan Lien
!!-----------------------------------------------------------------------
!SUBROUTINE set_efso_obs
!  IMPLICIT NONE
!  REAL(r_size),ALLOCATABLE :: tmpdep(:)
!  INTEGER,ALLOCATABLE :: tmpqc0(:,:)
!  INTEGER,ALLOCATABLE :: tmpqc(:)
!  INTEGER :: nj(0:nlat-1)
!  INTEGER :: njs(1:nlat-1)
!  INTEGER :: l,im,i,j,n,nn,ierr
!  CHARACTER(14) :: obsguesfile='obsguesNNN.dat'
!  CHARACTER(14) :: obsanalfile='obsanalNNN.dat'

!  WRITE(6,'(A)') 'Hello from set_efso_obs'

!  dist_zero = SIGMA_OBS * SQRT(10.0d0/3.0d0) * 2.0d0
!  dist_zero_rain = SIGMA_OBS_RAIN * SQRT(10.0d0/3.0d0) * 2.0d0
!  dist_zerov = SIGMA_OBSV * SQRT(10.0d0/3.0d0) * 2.0d0
!  dist_zerov_rain = SIGMA_OBSV_RAIN * SQRT(10.0d0/3.0d0) * 2.0d0

!  CALL get_nobs_mpi(obsanalfile,10,nobs)
!  WRITE(6,'(I10,A)') nobs,' TOTAL OBSERVATIONS INPUT'
!  IF(nobs == 0) RETURN
!!
!! INITIALIZE GLOBAL VARIABLES
!!
!  ALLOCATE( obselm(nobs) )
!  ALLOCATE( obslon(nobs) )
!  ALLOCATE( obslat(nobs) )
!  ALLOCATE( obslev(nobs) )
!  ALLOCATE( obsdat(nobs) )
!  ALLOCATE( obserr(nobs) )
!  ALLOCATE( obstyp(nobs) )
!  ALLOCATE( obsdif(nobs) )
!  ALLOCATE( obsdep(nobs) )
!  ALLOCATE( obshdxf(nobs,MEMBER) )
!  ALLOCATE( obsqc(nobs) )
!  ALLOCATE( tmpdep(nobs) )
!  ALLOCATE( tmpqc0(nobs,MEMBER) )
!!
!! reading background observation data and compute departure
!!
!  IF(myrank == 0) THEN
!    WRITE(obsguesfile(8:10),'(A3)') '_me'
!    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',obsguesfile
!    CALL read_obs2(obsguesfile,nobs,obselm,obslon,obslat,obslev, &
!                   obsdat,obserr,obstyp,obsdif,obsdep,obsqc)
!    obsdep = obsdat - obsdep
!  END IF
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  CALL MPI_BCAST(obsdep,nobs,MPI_r_size,0,MPI_COMM_WORLD,ierr)
!  CALL MPI_BCAST(obsqc,nobs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!
!! reading ensemble analysis observations
!!
!  CALL read_obs2_mpi(obsanalfile,nobs,MEMBER,obselm,obslon,obslat,obslev, &
!                     obsdat,obserr,obstyp,obsdif,obshdxf,tmpqc0)

!!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i)
!  DO n=1,nobs
!    tmpdep(n) = obshdxf(n,1)
!    DO i=2,MEMBER
!      tmpdep(n) = tmpdep(n) + obshdxf(n,i)
!    END DO
!    tmpdep(n) = tmpdep(n) / REAL(MEMBER,r_size) ! mean
!    DO i=1,MEMBER
!      obshdxf(n,i) = obshdxf(n,i) - tmpdep(n)
!    END DO
!  END DO
!!$OMP END PARALLEL DO
!  DEALLOCATE(tmpdep,tmpqc0)
!!
!! Create observation box
!!
!  nobsgrd = 0
!  nj = 0
!!$OMP PARALLEL PRIVATE(i,j,n,nn)
!!$OMP DO SCHEDULE(DYNAMIC)
!  DO j=1,nlat-1
!    DO n=1,nobs
!      IF(obslat(n) < lat(j) .OR. lat(j+1) <= obslat(n)) CYCLE
!      nj(j) = nj(j) + 1
!    END DO
!  END DO
!!$OMP END DO
!!$OMP DO SCHEDULE(DYNAMIC)
!  DO j=1,nlat-1
!    njs(j) = SUM(nj(0:j-1))
!  END DO
!!$OMP END DO
!!$OMP DO SCHEDULE(DYNAMIC)
!  DO j=1,nlat-1
!    nn = 0
!    DO n=1,nobs
!      IF(obslat(n) < lat(j) .OR. lat(j+1) <= obslat(n)) CYCLE
!      nn = nn + 1
!    END DO
!  END DO
!!$OMP END DO
!!$OMP DO SCHEDULE(DYNAMIC)
!  DO j=1,nlat-1
!    IF(nj(j) == 0) THEN
!      nobsgrd(:,j) = njs(j)
!      CYCLE
!    END IF
!    nn = 0
!    DO i=1,nlon
!      DO n=njs(j)+1,njs(j)+nj(j)
!        IF(i < nlon) THEN
!          IF(obslon(n) < lon(i) .OR. lon(i+1) <= obslon(n)) CYCLE
!        ELSE
!          IF(obslon(n) < lon(nlon) .OR. 360.0d0 <= obslon(n)) CYCLE
!        END IF
!        nn = nn + 1
!      END DO
!      nobsgrd(i,j) = njs(j) + nn
!    END DO
!    IF(nn /= nj(j)) THEN
!!$OMP CRITICAL
!      WRITE(6,'(A,2I10)') 'OBS DATA SORT ERROR: ',nn,nj(j)
!      WRITE(6,'(F6.2,A,F6.2)') lat(j),'< LAT <',lat(j+1)
!      WRITE(6,'(F6.2,A,F6.2)') MINVAL(obslat(njs(j)+1:njs(j)+nj(j))),'< OBSLAT <',MAXVAL(obslat(njs(j)+1:njs(j)+nj(j)))
!!$OMP END CRITICAL
!    END IF
!  END DO
!!$OMP END DO
!!$OMP END PARALLEL

!  RETURN
!END SUBROUTINE set_efso_obs

END MODULE letkf_obs
