MODULE obsope_tools
!=======================================================================
!
! [PURPOSE:] Observation operator tools
!
! [HISTORY:]
!   November 2014  Guo-Yuan Lien  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_scale
  USE common_mpi_scale
  USE common_obs_scale

!  use common_scalelib

  use common_nml

  use scale_process, only: &
       PRC_myrank
!       MPI_COMM_d => LOCAL_COMM_WORLD

  use scale_grid_index, only: &
    KHALO, IHALO, JHALO
#ifdef H08
  use scale_grid, only: &
      DX, DY,           &
      BUFFER_DX,        &
      BUFFER_DY
#endif


  IMPLICIT NONE
  PUBLIC

!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------



CONTAINS
!!-----------------------------------------------------------------------
!! Read namelist for obsope
!!-----------------------------------------------------------------------
!subroutine read_nml_obsope
!  implicit none

!  call read_nml_letkf
!  call read_nml_letkf_prc

!  return
!end subroutine read_nml_obsope
!!-----------------------------------------------------------------------
!! Read namelist for obsmake
!!-----------------------------------------------------------------------
!subroutine read_nml_obsmake
!  implicit none

!  call read_nml_letkf
!  call read_nml_letkf_prc
!  call read_nml_letkf_obsmake

!  return
!end subroutine read_nml_obsmake
!-----------------------------------------------------------------------
! PARAM_LETKF_OBSMAKE
!-----------------------------------------------------------------------
!subroutine read_nml_letkf_obsmake
!  implicit none
!  integer :: ierr

!  namelist /PARAM_LETKF_OBSMAKE/ &
!    OBSERR_U, &
!    OBSERR_V, &
!    OBSERR_T, &
!    OBSERR_Q, &
!    OBSERR_RH, &
!    OBSERR_PS, &
!    OBSERR_RADAR_REF, &
!    OBSERR_RADAR_VR

!  rewind(IO_FID_CONF)
!  read(IO_FID_CONF,nml=PARAM_LETKF_OBSMAKE,iostat=ierr)
!  if (ierr < 0) then !--- missing
!    write(6,*) 'xxx Not found namelist. Check!'
!    stop
!  elseif (ierr > 0) then !--- fatal error
!    write(6,*) 'xxx Not appropriate names in namelist LETKF_PARAM_OBSMAKE. Check!'
!    stop
!  endif

!  return
!end subroutine read_nml_letkf_obsmake

!-----------------------------------------------------------------------
! Observation operator calculation
!-----------------------------------------------------------------------
SUBROUTINE obsope_cal(obs, obsda_return)
  IMPLICIT NONE

  TYPE(obs_info),INTENT(IN) :: obs(OBS_IN_NUM)
  type(obs_da_value),OPTIONAL,INTENT(INOUT) :: obsda_return
  type(obs_da_value) :: obsda
  REAL(r_size),ALLOCATABLE :: v3dg(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2dg(:,:,:)

  integer :: it,islot,proc,im,iof
  integer :: n,nn,nslot,nobs,nobs_0,nobs_slot,nobs_alldomain
!  real(r_size) :: rig,rjg,ri,rj,rk
  real(r_size) :: rig,rjg,rk
  real(r_size),allocatable :: ri(:),rj(:)
  real(r_size) :: ritmp,rjtmp
  real(r_size) :: slot_lb,slot_ub

#ifdef H08
! -- for Himawari-8 obs --
  INTEGER :: nallprof ! H08: Num of all profiles (entire domain) required by RTTOV
  INTEGER :: ns ! H08 obs count
  INTEGER :: nprof_H08 ! num of H08 obs
  REAL(r_size),ALLOCATABLE :: ri_H08(:),rj_H08(:)
  REAL(r_size),ALLOCATABLE :: lon_H08(:),lat_H08(:)
  REAL(r_size),ALLOCATABLE :: tmp_ri_H08(:),tmp_rj_H08(:)
  REAL(r_size),ALLOCATABLE :: tmp_lon_H08(:),tmp_lat_H08(:)

  REAL(r_size),ALLOCATABLE :: yobs_H08(:),plev_obs_H08(:)
  REAL(r_size),ALLOCATABLE :: yobs_H08_clr(:)
  INTEGER :: ch
  INTEGER,ALLOCATABLE :: qc_H08(:)

! -- Rejecting obs over the buffer regions. --
!
! bris: "ri" at the wetern end of the domain excluding buffer regions
! brie: "ri" at the eastern end of the domain excluding buffer regions
! bris: "rj" at the southern end of the domain excluding buffer regions
! bris: "rj" at the northern end of the domain excluding buffer regions
!
! e.g.,   ri:    ...bris...........brie...
!             buffer |  NOT buffer  | buffer
!
!
  REAL(r_size) :: bris, brie
  REAL(r_size) :: brjs, brje
#endif

! -- for TC vital assimilation --
  INTEGER :: obs_idx_TCX, obs_idx_TCY, obs_idx_TCP ! obs index
  INTEGER :: bTC_proc ! the process where the background TC is located.
! bTC: background TC in each subdomain
! bTC(1,:) : tcx (m), bTC(2,:): tcy (m), bTC(3,:): mslp (Pa)
  REAL(r_size),ALLOCATABLE :: bTC(:,:)
  REAL(r_size),ALLOCATABLE :: bufr(:,:)
  REAL(r_size) :: bTC_mslp

  character(filelenmax) :: obsdafile
  character(11) :: obsda_suffix = '.000000.dat'

!-----------------------------------------------------------------------

  integer :: ierr
  REAL(r_dble) :: rrtimer00,rrtimer

!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer00 = MPI_WTIME()

#ifdef H08
!  call phys2ij(MSLP_TC_LON,MSLP_TC_LAT,MSLP_TC_rig,MSLP_TC_rjg)
  bris = real(BUFFER_DX/DX,r_size) + real(IHALO,r_size) 
  brjs = real(BUFFER_DY/DY,r_size) + real(JHALO,r_size)
  brie = (real(nlong+2*IHALO,r_size) - bris)
  brje = (real(nlatg+2*JHALO,r_size) - brjs)
#endif

  nobs_alldomain = 0
  do iof = 1, OBS_IN_NUM
    if (OBSDA_RUN(iof)) then
      nobs_alldomain = nobs_alldomain + obs(iof)%nobs
    end if
  end do
  obsda%nobs = nobs_alldomain
  call obs_da_value_allocate(obsda,0)
  allocate ( ri(nobs_alldomain) )
  allocate ( rj(nobs_alldomain) )

  allocate ( v3dg (nlevh,nlonh,nlath,nv3dd) )
  allocate ( v2dg (nlonh,nlath,nv2dd) )

  do it = 1, nitmax
    im = proc2mem(1,it,myrank+1)
    if (im >= 1 .and. im <= MEMBER) then
      write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is processing member ', &
            im, ', subdomain id #', proc2mem(2,it,myrank+1)

!write(6,*) '%%%%%%', MPI_WTIME(), 0

      nobs = 0

      !!!!!!
      if (nobs_alldomain > 0) then
      !!!!!!

      obsda%qc = iqc_time

      do islot = SLOT_START, SLOT_END
        slot_lb = (real(islot-SLOT_BASE,r_size) - 0.5d0) * SLOT_TINTERVAL
        slot_ub = (real(islot-SLOT_BASE,r_size) + 0.5d0) * SLOT_TINTERVAL
        write (6,'(A,I3,A,F9.1,A,F9.1,A)') 'Slot #', islot-SLOT_START+1, ': time interval (', slot_lb, ',', slot_ub, '] sec'

        call read_ens_history_iter(HISTORY_IN_BASENAME,it,islot,v3dg,v2dg)
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer = MPI_WTIME()
  WRITE(6,'(A,I3,A,I3,A,4x,F15.7)') '###### obsope_cal:read_ens_history_iter:',it,':',islot,':',rrtimer-rrtimer00
  rrtimer00=rrtimer


        do iof = 1, OBS_IN_NUM

          if (.not. OBSDA_RUN(iof)) cycle

          obs_idx_TCX = -1
          obs_idx_TCY = -1
          obs_idx_TCP = -1

          nslot = 0
          nobs_slot = 0

!write(6,*) '%%%===', MPI_WTIME(), 'im:', im, 'islot:', islot, 'iof:', iof

          IF(OBS_IN_FORMAT(iof) /= 3)THEN ! except H08 obs ! H08

            ! do this small computation first, without OpenMP
            nobs_0 = nobs
            do n = 1, obs(iof)%nobs

              select case (obs(iof)%elm(n))
              case (id_tclon_obs)
                obs_idx_TCX = n
                cycle
              case (id_tclat_obs)
                obs_idx_TCY = n
                cycle
              case (id_tcmip_obs)
                obs_idx_TCP = n
                cycle
              end select

              if (obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub) then
                nslot = nslot + 1
                call phys2ij(obs(iof)%lon(n),obs(iof)%lat(n),rig,rjg)
                call rij_g2l_auto(proc,rig,rjg,ritmp,rjtmp)

                if (PRC_myrank == proc) then
                  nobs = nobs + 1
                  nobs_slot = nobs_slot + 1
                  obsda%set(nobs) = iof
                  obsda%idx(nobs) = n
                  obsda%ri(nobs) = rig
                  obsda%rj(nobs) = rjg
                  ri(nobs) = ritmp
                  rj(nobs) = rjtmp
                end if ! [ PRC_myrank == proc ]
              end if ! [ obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub ]
            end do ! [ n = 1, obs%nobs ]

#ifdef H08
          ELSEIF( OBS_IN_FORMAT(iof) == 3) THEN ! for H08 obs (OBS_IN_FORMAT(iof) = 3) ! H08

            nprof_H08 = 0
            nobs_0 = nobs
            nallprof = obs(iof)%nobs/nch

            ALLOCATE(tmp_ri_H08(nallprof))
            ALLOCATE(tmp_rj_H08(nallprof))
            ALLOCATE(tmp_lon_H08(nallprof))
            ALLOCATE(tmp_lat_H08(nallprof))

            do n = 1, nallprof
              ns = (n - 1) * nch + 1
              if (obs(iof)%dif(ns) > slot_lb .and. obs(iof)%dif(ns) <= slot_ub) then
                nslot = nslot + 1
                call phys2ij(obs(iof)%lon(ns),obs(iof)%lat(ns),rig,rjg)
                call rij_g2l_auto(proc,rig,rjg,ritmp,rjtmp)

                if (PRC_myrank == proc) then
                  nprof_H08 = nprof_H08 + 1 ! num of prof in myrank node
                  tmp_ri_H08(nprof_H08) = ritmp
                  tmp_rj_H08(nprof_H08) = rjtmp
                  tmp_lon_H08(nprof_H08) = obs(iof)%lon(ns)
                  tmp_lat_H08(nprof_H08) = obs(iof)%lat(ns)

                  nobs = nobs + nch
                  nobs_slot = nobs_slot + 1
                  obsda%set(nobs-nch+1:nobs) = iof
                  obsda%ri(nobs-nch+1:nobs) = rig
                  obsda%rj(nobs-nch+1:nobs) = rjg
                  ri(nobs-nch+1:nobs) = ritmp
                  rj(nobs-nch+1:nobs) = rjtmp
                  do ch = 1, nch
                    obsda%idx(nobs-nch+ch) = ns + ch - 1
                  enddo

                end if ! [ PRC_myrank == proc ]
              end if ! [ obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub ]
            end do ! [ n = 1, nallprof ]

            IF(nprof_H08 >=1)THEN
              ALLOCATE(ri_H08(nprof_H08))
              ALLOCATE(rj_H08(nprof_H08))
              ALLOCATE(lon_H08(nprof_H08))
              ALLOCATE(lat_H08(nprof_H08))

              ri_H08 = tmp_ri_H08(1:nprof_H08)
              rj_H08 = tmp_rj_H08(1:nprof_H08)
              lon_H08 = tmp_lon_H08(1:nprof_H08)
              lat_H08 = tmp_lat_H08(1:nprof_H08)

            ENDIF

            DEALLOCATE(tmp_ri_H08,tmp_rj_H08)
            DEALLOCATE(tmp_lon_H08,tmp_lat_H08)

#endif
          ENDIF ! end of nobs count [if (OBS_IN_FORMAT(iof) = 3)]



!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer = MPI_WTIME()
  WRITE(6,'(A,I3,A,I3,A,I3,A,F15.7)') '###### obsope_cal:obsope_step_1:        ',it,':',islot,':',iof,':',rrtimer-rrtimer00
  rrtimer00=rrtimer


          ! then do this heavy computation with OpenMP

          IF(OBS_IN_FORMAT(iof) /= 3)THEN ! H08


!write(6,*) '%%%===', MPI_WTIME(), nobs_0 + 1, nobs

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(nn,n,rk)
            do nn = nobs_0 + 1, nobs
              n = obsda%idx(nn)

!if (mod(nn,50) == 0) then
!  write(6,*) '%%%%%%', MPI_WTIME(), nn
!end if

              if (obs(iof)%elm(n) == id_radar_ref_obs .or. obs(iof)%elm(n) == id_radar_vr_obs) then
                if (obs(iof)%lev(n) > RADAR_ZMAX) then
                  obsda%qc(nn) = iqc_radar_vhi
                  write(6,'(A,F8.1,A,I5)') 'warning: radar observation is too high: lev=', obs(iof)%lev(n), ', elem=', obs(iof)%elm(n)
                else
                  call phys2ijkz(v3dg(:,:,:,iv3dd_hgt),ri(nn),rj(nn),obs(iof)%lev(n),rk,obsda%qc(nn))
                end if
              else
                call phys2ijk(v3dg(:,:,:,iv3dd_p),obs(iof)%elm(n),ri(nn),rj(nn),obs(iof)%lev(n),rk,obsda%qc(nn))
              end if


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,I3,A,I3,A,I3,A,I8,A,F15.7)') '###### obsope_cal:obsope_step_2_phys2ijkz:',it,':',islot,':',iof,':',nn,':',rrtimer-rrtimer00
!  rrtimer00=rrtimer


              if (obsda%qc(nn) == iqc_good) then
                select case (OBS_IN_FORMAT(iof))
                case (1)
                  call Trans_XtoY(obs(iof)%elm(n),ri(nn),rj(nn),rk, &
                                  obs(iof)%lon(n),obs(iof)%lat(n),v3dg,v2dg,obsda%val(nn),obsda%qc(nn))
                case (2)
                  call Trans_XtoY_radar(obs(iof)%elm(n),obs(iof)%meta(1),obs(iof)%meta(2),obs(iof)%meta(3),ri(nn),rj(nn),rk, &
                                        obs(iof)%lon(n),obs(iof)%lat(n),obs(iof)%lev(n),v3dg,v2dg,obsda%val(nn),obsda%qc(nn))
                  if (obsda%qc(nn) == iqc_ref_low) obsda%qc(nn) = iqc_good ! when process the observation operator, we don't care if reflectivity is too small

                !!!!!! may not need to do this at this stage...
                !if (obs(iof)%elm(n) == id_radar_ref_obs) then
                !  obsda%val(nn) = 10.0d0 * log10(obsda%val(nn))
                !end if
                !!!!!!

                end select
              end if


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
!  rrtimer = MPI_WTIME()
!  WRITE(6,'(A,I3,A,I3,A,I3,A,I8,A,F15.7)') '###### obsope_cal:obsope_step_2_Trans_XtoY_radar:',it,':',islot,':',iof,':',nn,':',rrtimer-rrtimer00
!  rrtimer00=rrtimer


            end do ! [ nn = nobs_0 + 1, nobs ]
!$OMP END PARALLEL DO

#ifdef H08
          ELSEIF((OBS_IN_FORMAT(iof) == 3).and.(nprof_H08 >=1 ))THEN ! H08
! -- Note: Trans_XtoY_H08 is called without OpenMP but it can use a parallel (with OpenMP) RTTOV routine
!

            ALLOCATE(yobs_H08(nprof_H08*nch))
            ALLOCATE(yobs_H08_clr(nprof_H08*nch))
            ALLOCATE(plev_obs_H08(nprof_H08*nch))
            ALLOCATE(qc_H08(nprof_H08*nch))

            CALL Trans_XtoY_H08(nprof_H08,ri_H08,rj_H08,&
                                lon_H08,lat_H08,v3dg,v2dg,&
                                yobs_H08,plev_obs_H08,&
                                qc_H08,yobs_H08_clr=yobs_H08_clr)

! Clear sky yobs(>0)
! Cloudy sky yobs(<0)

            obsda%qc(nobs_0+1:nobs) = iqc_obs_bad

            ns = 0
            DO nn = nobs_0 + 1, nobs
              ns = ns + 1

              obsda%val(nn) = yobs_H08(ns)
              obsda%qc(nn) = qc_H08(ns)

              if(obsda%qc(nn) == iqc_good)then
                rig = obsda%ri(nn)
                rjg = obsda%rj(nn)

! -- tentative treatment around the TC center --
!                dist_MSLP_TC = sqrt(((rig - MSLP_TC_rig) * DX)**2&
!                                   +((rjg - MSLP_TC_rjg) * DY)**2)

!                if(dist_MSLP_TC <= dist_MSLP_TC_MIN)then
!                  obsda%qc(nn) = iqc_obs_bad
!                endif

! -- Rejecting Himawari-8 obs over the buffer regions. --
                if((rig <= bris) .or. (rig >= brie) .or.&
                   (rjg <= brjs) .or. (rjg >= brje))then
                  obsda%qc(nn) = iqc_obs_bad
                endif
              endif

!
!  NOTE: T.Honda (10/16/2015)
!  The original H08 obs does not inlcude the level information.
!  However, we have the level information derived by RTTOV (plev_obs_H08) here, 
!  so that we substitute the level information into obsda%lev.  
!  The substituted level information is used in letkf_tools.f90
!
              obsda%lev(nn) = plev_obs_H08(ns)
              obsda%val2(nn) = yobs_H08_clr(ns)

!              write(6,'(a,f12.1,i9)')'H08 debug_plev',obsda%lev(nn),nn

            END DO ! [ nn = nobs_0 + 1, nobs ]

            DEALLOCATE(ri_H08, rj_H08)
            DEALLOCATE(lon_H08, lat_H08)
            DEALLOCATE(yobs_H08, plev_obs_H08)
            DEALLOCATE(yobs_H08_clr)
            DEALLOCATE(qc_H08)

#endif
          ENDIF ! H08

          write (6,'(3A,I10)') ' -- [', trim(OBS_IN_NAME(iof)), '] nobs in the slot = ', nslot
          write (6,'(3A,I6,A,I10)') ' -- [', trim(OBS_IN_NAME(iof)), '] nobs in the slot and processed by rank ' &
                                    , myrank, ' = ', nobs_slot


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer = MPI_WTIME()
  WRITE(6,'(A,I3,A,I3,A,I3,A,F15.7)') '###### obsope_cal:obsope_step_2:        ',it,':',islot,':',iof,':',rrtimer-rrtimer00
  rrtimer00=rrtimer



! ###  -- TC vital assimilation -- ###
          if (obs_idx_TCX > 0 .and. obs_idx_TCY > 0 .and. obs_idx_TCP > 0 .and. &
              obs(iof)%dif(obs_idx_TCX) == obs(iof)%dif(obs_idx_TCY) .and. &
              obs(iof)%dif(obs_idx_TCY) == obs(iof)%dif(obs_idx_TCP)) then
           
            if (obs(iof)%dif(obs_idx_TCX) > slot_lb .and. &
              obs(iof)%dif(obs_idx_TCX) <= slot_ub) then
              nslot = nslot + 3 ! TC vital obs should have 3 data (i.e., lon, lat, and MSLP)

              !!! bTC(1,:) : lon, bTC(2,:): lat, bTC(3,:): mslp
              ! bTC(1,:) : tcx (m), bTC(2,:): tcy (m), bTC(3,:): mslp
              allocate(bTC(3,0:MEM_NP-1))
              allocate(bufr(3,0:MEM_NP-1))

              bTC = 9.99d33
              bufr = 9.99d33

              ! Note: obs(iof)%dat(obs_idx_TCX) is not longitude (deg) but X (m).
              !       Units of the original TC vital position are converted in
              !       subroutine read_obs in common_obs_scale.f90.
              !
              call phys2ij(obs(iof)%lon(obs_idx_TCX),obs(iof)%lat(obs_idx_TCX),rig,rjg) 
              call rij_g2l_auto(proc,rig,rjg,ritmp,rjtmp)  
              call search_tc_subdom(rig,rjg,v2dg,bTC(1,PRC_myrank),bTC(2,PRC_myrank),bTC(3,PRC_myrank))
  
              CALL MPI_BARRIER(MPI_COMM_d,ierr)
              CALL MPI_ALLREDUCE(bTC,bufr,3*MEM_NP,MPI_r_size,MPI_MIN,MPI_COMM_d,ierr)
              bTC = bufr

              deallocate(bufr)

              ! Assume MSLP of background TC is lower than 1100 (hPa). 
              bTC_mslp = 1100.0d2
              do n = 0, MEM_NP - 1
                write(6,'(3e20.5)')bTC(1,n),bTC(2,n),bTC(3,n) ! debug
                if (bTC(3,n) < bTC_mslp ) then
                  bTC_mslp = bTC(3,n)
                  bTC_proc = n
                endif
              enddo ! [ n = 0, MEM_NP - 1]

              if (PRC_myrank == proc) then
                do n = 1, 3
                  nobs = nobs + 1
                  nobs_slot = nobs_slot + 1
                  obsda%set(nobs) = iof
                  if(n==1) obsda%idx(nobs) = obs_idx_TCX
                  if(n==2) obsda%idx(nobs) = obs_idx_TCY
                  if(n==3) obsda%idx(nobs) = obs_idx_TCP
                  obsda%ri(nobs) = rig
                  obsda%rj(nobs) = rjg
                  ri(nobs) = ritmp
                  rj(nobs) = rjtmp

                  obsda%val(nobs) = bTC(n,bTC_proc)
                  obsda%qc(nobs) = iqc_good
                enddo ! [ n = 1, 3 ]

              endif
              deallocate(bTC)

            endif ! [ obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub ]
          endif ! [ obs_idx_TCX > 0 ...]

        end do ! [ do iof = 1, OBS_IN_NUM ]

 
!      IF(NINT(elem(n)) == id_ps_obs .AND. odat(n) < -100.0d0) THEN
!        CYCLE
!      END IF
!      IF(NINT(elem(n)) == id_ps_obs) THEN
!        CALL itpl_2d(v2d(:,:,iv2d_orog),ri,rj,dz)
!        rk = rlev(n) - dz
!        IF(ABS(rk) > threshold_dz) THEN ! pressure adjustment threshold
!!          WRITE(6,'(A)') '* PS obs vertical adjustment beyond threshold'
!!          WRITE(6,'(A,F10.2,A,F6.2,A,F6.2,A)') '*   dz=',rk,&
!!           & ', (lon,lat)=(',elon(n),',',elat(n),')'
!          CYCLE
!        END IF
!      END IF

      end do ! [ islot = SLOT_START, SLOT_END ]

!write(6,*) '%%%%%%', MPI_WTIME(), nobs

      !!!!!!
      end if ! [ nobs_alldomain > 0 ]
      !!!!!!

      if (it == 1) then
        obsda%nobs = nobs
      else if (nobs /= obsda%nobs) then
        write (6, '(A)') '[Error] numbers of observations found are different among members.'
        stop
      end if

      write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' finishes processing member ', &
            im, ', subdomain id #', proc2mem(2,it,myrank+1)
      write (6,'(A,I8,A)') ' -- ', nobs, ' observations found'


      if (OBSDA_OUT) then
        write (6,'(A,I6.6,A,I4.4,A,I6.6)') 'MYRANK ',myrank,' is writing observations for member ', &
              im, ', subdomain id #', proc2mem(2,it,myrank+1)
        call file_member_replace(im, OBSDA_OUT_BASENAME, obsdafile)
        write (obsda_suffix(2:7),'(I6.6)') proc2mem(2,it,myrank+1)
        call write_obs_da(trim(obsdafile)//obsda_suffix,obsda,0)
      end if

      if (present(obsda_return)) then
        ! variables without an ensemble dimension
        if (it == 1) then
          obsda_return%nobs = nobs + obsda_return%nobs ! obsda_return%nobs: additional space for externally processed observations
          call obs_da_value_allocate(obsda_return,MEMBER)
        end if
        if (nobs > 0) then
          if (it == 1) then
            obsda_return%set(1:nobs) = obsda%set(1:nobs)
            obsda_return%idx(1:nobs) = obsda%idx(1:nobs)
            obsda_return%ri(1:nobs) = obsda%ri(1:nobs)
            obsda_return%rj(1:nobs) = obsda%rj(1:nobs)
            obsda_return%qc(1:nobs) = obsda%qc(1:nobs)
#ifdef H08
            obsda_return%lev(1:nobs) = obsda%lev(1:nobs)
            obsda_return%val2(1:nobs) = obsda%val2(1:nobs)
#endif
          else
            obsda_return%qc(1:nobs) = max(obsda_return%qc(1:nobs), obsda%qc(1:nobs))
#ifdef H08
            obsda_return%lev(1:nobs) = obsda_return%lev(1:nobs) + obsda%lev(1:nobs)
            obsda_return%val2(1:nobs) = obsda_return%val2(1:nobs) + obsda%val2(1:nobs)
#endif
          end if

          ! variables with an ensemble dimension
          obsda_return%ensval(im,1:nobs) = obsda%val(1:nobs)
        end if ! [ nobs > 0 ]
      end if ! [ present(obsda_return) ]


!  CALL MPI_BARRIER(MPI_COMM_a,ierr)
  rrtimer = MPI_WTIME()
  WRITE(6,'(A,I3,A,8x,F15.7)') '###### obsope_cal:write_obs_da:         ',it,':',rrtimer-rrtimer00
  rrtimer00=rrtimer


    end if ! [ im >= 1 .and. im <= MEMBER ]

  end do ! [ it = 1, nitmax ]

  call obs_da_value_deallocate(obsda)

  deallocate ( ri, rj, v3dg, v2dg )

!write(6,*) ri(1),rj(1)
!write(6,*) '$$$$ 0'
!!  deallocate ( ri )
!write(6,*) '$$$$ 1'
!!  deallocate ( rj )
!write(6,*) '$$$$ 2'
!!  deallocate ( v3dg )
!write(6,*) '$$$$ 3'
!!  deallocate ( v2dg )
!write(6,*) '$$$$ 4'

end subroutine obsope_cal

!-----------------------------------------------------------------------
! Observation generator calculation
!-----------------------------------------------------------------------
SUBROUTINE obsmake_cal(obs)
  IMPLICIT NONE

  TYPE(obs_info),INTENT(INOUT) :: obs(OBS_IN_NUM)
  REAL(r_size),ALLOCATABLE :: v3dg(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2dg(:,:,:)

  integer :: islot,proc
  integer :: n,nslot,nobs,nobs_slot,ierr,iqc,iof
  integer :: nobsmax,nobsall
  real(r_size) :: rig,rjg,ri,rj,rk
  real(r_size) :: slot_lb,slot_ub
  real(r_size),allocatable :: bufr(:)
  real(r_size),allocatable :: error(:)

  CHARACTER(10) :: obsoutfile = 'obsout.dat'
  INTEGER :: ns 
#ifdef H08
! obsmake for H08 is not available !! (03/17/2016) T.Honda
! -- for Himawari-8 obs --
  INTEGER :: nallprof ! H08: Num of all profiles (entire domain) required by RTTOV
  INTEGER :: nprof_H08 ! num of H08 obs
  REAL(r_size),ALLOCATABLE :: ri_H08(:),rj_H08(:)
  REAL(r_size),ALLOCATABLE :: lon_H08(:),lat_H08(:)
  REAL(r_size),ALLOCATABLE :: tmp_ri_H08(:),tmp_rj_H08(:)
  REAL(r_size),ALLOCATABLE :: tmp_lon_H08(:),tmp_lat_H08(:)

  REAL(r_size),ALLOCATABLE :: yobs_H08(:),plev_obs_H08(:)
  INTEGER,ALLOCATABLE :: qc_H08(:)
  INTEGER,ALLOCATABLE :: idx_H08(:) ! index array
  INTEGER :: ich
#endif

!-----------------------------------------------------------------------

  write (6,'(A,I6.6,A,I6.6)') 'MYRANK ',myrank,' is processing subdomain id #', proc2mem(2,1,myrank+1)

  allocate ( v3dg (nlevh,nlonh,nlath,nv3dd) )
  allocate ( v2dg (nlonh,nlath,nv2dd) )

  do iof = 1, OBS_IN_NUM
    obs(iof)%dat = 0.0d0
  end do

  nobs = 0
  do islot = SLOT_START, SLOT_END
    slot_lb = (real(islot-SLOT_BASE,r_size) - 0.5d0) * SLOT_TINTERVAL
    slot_ub = (real(islot-SLOT_BASE,r_size) + 0.5d0) * SLOT_TINTERVAL
    write (6,'(A,I3,A,F9.1,A,F9.1,A)') 'Slot #', islot-SLOT_START+1, ': time interval (', slot_lb, ',', slot_ub, '] sec'

    call read_ens_history_iter(HISTORY_IN_BASENAME,1,islot,v3dg,v2dg)

    do iof = 1, OBS_IN_NUM
      IF(OBS_IN_FORMAT(iof) /= 3)THEN ! except H08 obs
        nslot = 0
        nobs_slot = 0
        do n = 1, obs(iof)%nobs

          if (obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub) then
            nslot = nslot + 1

            call phys2ij(obs(iof)%lon(n),obs(iof)%lat(n),rig,rjg)
            call rij_g2l_auto(proc,rig,rjg,ri,rj)

  !          if (PRC_myrank == 0) then
  !            print *, proc, rig, rjg, ri, rj
  !          end if

            if (proc < 0 .and. PRC_myrank == 0) then ! if outside of the global domain, processed by PRC_myrank == 0
              obs(iof)%dat(n) = undef
            end if

            if (PRC_myrank == proc) then
              nobs = nobs + 1
              nobs_slot = nobs_slot + 1

  !IF(NINT(elem(n)) == id_ps_obs) THEN
  !  CALL itpl_2d(v2d(:,:,iv2d_orog),ri,rj,dz)
  !  rk = rlev(n) - dz
  !  IF(ABS(rk) > threshold_dz) THEN ! pressure adjustment threshold
  !    ! WRITE(6,'(A)') '* PS obs vertical adjustment beyond threshold'
  !    ! WRITE(6,'(A,F10.2,A,F6.2,A,F6.2,A)') '* dz=',rk,&
  !    ! & ', (lon,lat)=(',elon(n),',',elat(n),')'
  !    CYCLE
  !  END IF
  !END IF

              if (obs(iof)%elm(n) == id_radar_ref_obs .or. obs(iof)%elm(n) == id_radar_vr_obs) then
                call phys2ijkz(v3dg(:,:,:,iv3dd_hgt),ri,rj,obs(iof)%lev(n),rk,iqc)
              else
                call phys2ijk(v3dg(:,:,:,iv3dd_p),obs(iof)%elm(n),ri,rj,obs(iof)%lev(n),rk,iqc)
              end if

              if (iqc /= iqc_good) then
                obs(iof)%dat(n) = undef
              else
                select case (OBS_IN_FORMAT(iof))
                case (1)
                  call Trans_XtoY(obs(iof)%elm(n),ri,rj,rk, &
                                  obs(iof)%lon(n),obs(iof)%lat(n),v3dg,v2dg,obs(iof)%dat(n),iqc)
                case (2)
                  call Trans_XtoY_radar(obs(iof)%elm(n),obs(iof)%meta(1),obs(iof)%meta(2),obs(iof)%meta(3),ri,rj,rk, &
                                        obs(iof)%lon(n),obs(iof)%lat(n),obs(iof)%lev(n),v3dg,v2dg,obs(iof)%dat(n),iqc)
                end select

 !!! For radar observation, when reflectivity value is too low, do not generate ref/vr observations
 !!! No consideration of the terrain blocking effects.....

                if (iqc /= iqc_good) then
                  obs(iof)%dat(n) = undef
                end if
              end if

            end if ! [ PRC_myrank == proc ]

          end if ! [ obs%dif(n) > slot_lb .and. obs%dif(n) <= slot_ub ]

        end do ! [ n = 1, obs%nobs ]

#ifdef H08
! -- H08 part --
      ELSEIF(OBS_IN_FORMAT(iof) == 3)THEN ! H08
        nslot = 0
        nobs_slot = 0
        nprof_H08 = 0

        nallprof = obs(iof)%nobs/nch

        ALLOCATE(tmp_ri_H08(nallprof))
        ALLOCATE(tmp_rj_H08(nallprof))
        ALLOCATE(tmp_lon_H08(nallprof))
        ALLOCATE(tmp_lat_H08(nallprof))
        ALLOCATE(idx_H08(nallprof))

        do n = 1, nallprof
          ns = (n - 1) * nch + 1
          if (obs(iof)%dif(n) > slot_lb .and. obs(iof)%dif(n) <= slot_ub) then
            nslot = nslot + 1

            call phys2ij(obs(iof)%lon(ns),obs(iof)%lat(ns),rig,rjg)
            call rij_g2l_auto(proc,rig,rjg,ri,rj)


            if (proc < 0 .and. PRC_myrank == 0) then ! if outside of the global domain, processed by PRC_myrank == 0
              obs(iof)%dat(ns:ns+nch-1) = undef
            end if

            if (PRC_myrank == proc) then
              nprof_H08 = nprof_H08 + 1 ! num of prof in myrank node
              idx_H08(nprof_H08) = ns ! idx of prof in myrank node
              tmp_ri_H08(nprof_H08) = ri
              tmp_rj_H08(nprof_H08) = rj
              tmp_lon_H08(nprof_H08) = obs(iof)%lon(ns)
              tmp_lat_H08(nprof_H08) = obs(iof)%lat(ns)

              nobs = nobs + nch
              nobs_slot = nobs_slot + nch

            end if ! [ PRC_myrank == proc ]

          end if ! [ obs%dif(n) > slot_lb .and. obs%dif(n) <= slot_ub ]

        end do ! [ n = 1, nallprof ]

        IF(nprof_H08 >=1)THEN
          ALLOCATE(ri_H08(nprof_H08))
          ALLOCATE(rj_H08(nprof_H08))
          ALLOCATE(lon_H08(nprof_H08))
          ALLOCATE(lat_H08(nprof_H08))

          ri_H08 = tmp_ri_H08(1:nprof_H08)
          rj_H08 = tmp_rj_H08(1:nprof_H08)
          lon_H08 = tmp_lon_H08(1:nprof_H08)
          lat_H08 = tmp_lat_H08(1:nprof_H08)

          ALLOCATE(yobs_H08(nprof_H08*nch))
          ALLOCATE(plev_obs_H08(nprof_H08*nch))
          ALLOCATE(qc_H08(nprof_H08*nch))

          CALL Trans_XtoY_H08(nprof_H08,ri_H08,rj_H08,&
                              lon_H08,lat_H08,v3dg,v2dg,&
                              yobs_H08,plev_obs_H08,&
                              qc_H08)

          DO n = 1, nprof_H08
            ns = idx_H08(n)

            obs(iof)%lon(ns:ns+nch-1)=lon_H08(n:n+nch-1)
            obs(iof)%lat(ns:ns+nch-1)=lat_H08(n:n+nch-1)

            DO ich = 1, nch-1
              IF(qc_H08(n+ich-1) == iqc_good)THEN
                obs(iof)%dat(ns+ich-1)=undef
              ELSE
                obs(iof)%dat(ns+ich-1)=yobs_H08(n+ich-1)
              ENDIF
            ENDDO
          ENDDO

        ENDIF

        DEALLOCATE(tmp_ri_H08,tmp_rj_H08)
        DEALLOCATE(tmp_lon_H08,tmp_lat_H08)


! -- end of H08 part --
#endif
      ENDIF

    end do ! [ iof = 1, OBS_IN_NUM ]

    write (6,'(3A,I10)') ' -- [', trim(OBS_IN_NAME(iof)), '] nobs in the slot = ', nslot
    write (6,'(3A,I6,A,I10)') ' -- [', trim(OBS_IN_NAME(iof)), '] nobs in the slot and processed by rank ', myrank, ' = ', nobs_slot

  end do ! [ islot = SLOT_START, SLOT_END ]

  deallocate ( v3dg, v2dg )

  if (PRC_myrank == 0) then
    nobsmax = 0
    nobsall = 0
    do iof = 1, OBS_IN_NUM
      if (obs(iof)%nobs > nobsmax) nobsmax = obs(iof)%nobs
      nobsall = nobsall + obs(iof)%nobs
    end do

    allocate ( bufr(nobsmax) )
    allocate ( error(nobsall) )

    call com_randn(nobsall, error) ! generate all random numbers at the same time
    ns = 0
  end if

  do iof = 1, OBS_IN_NUM

    call MPI_REDUCE(obs(iof)%dat,bufr(1:obs(iof)%nobs),obs(iof)%nobs,MPI_r_size,MPI_SUM,0,MPI_COMM_d,ierr)

    if (PRC_myrank == 0) then
      obs(iof)%dat = bufr(1:obs(iof)%nobs)

      do n = 1, obs(iof)%nobs
        select case(obs(iof)%elm(n))
        case(id_u_obs)
          obs(iof)%err(n) = OBSERR_U
        case(id_v_obs)
          obs(iof)%err(n) = OBSERR_V
        case(id_t_obs,id_tv_obs)
          obs(iof)%err(n) = OBSERR_T
        case(id_q_obs)
          obs(iof)%err(n) = OBSERR_Q
        case(id_rh_obs)
          obs(iof)%err(n) = OBSERR_RH
        case(id_ps_obs)
          obs(iof)%err(n) = OBSERR_PS
        case(id_radar_ref_obs)
          obs(iof)%err(n) = OBSERR_RADAR_REF
        case(id_radar_vr_obs)
          obs(iof)%err(n) = OBSERR_RADAR_VR
!
! -- Not available (02/09/2015)
!        case(id_H08IR_obs) ! H08
!          obs(iof)%err(n) = OBSERR_H08(ch) !H08
!        case default
          write(6,'(A)') 'warning: skip assigning observation error (unsupported observation type)' 
        end select

        if (obs(iof)%dat(n) /= undef .and. obs(iof)%err(n) /= undef) then
          obs(iof)%dat(n) = obs(iof)%dat(n) + obs(iof)%err(n) * error(ns+n)
        end if

!print *, '######', obs%elm(n), obs%dat(n)
      end do ! [ n = 1, obs(iof)%nobs ]

      ns = ns + obs(iof)%nobs
    end if ! [ PRC_myrank == 0 ]

  end do ! [ iof = 1, OBS_IN_NUM ]

  if (PRC_myrank == 0) then
    deallocate ( bufr )
    deallocate ( error )

    call write_obs_all(obs, missing=.false., file_suffix='.out') ! only at the head node
  end if

end subroutine obsmake_cal
!=======================================================================

END MODULE obsope_tools
