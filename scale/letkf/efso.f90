PROGRAM efso
!=======================================================================
!
! [PURPOSE:] Main program of forecast sensitivity to observations using LETKF
!
! [HISTORY:]
!   09/29/2011 Yoichiro Ohta     created from main program of LETKF
!   07/01/2013 Daisuke Hotta     ported to GFS-LEKTF system
!   12/19/2013 Guo-Yuan Lien     merged to GFS-LETKF main development
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_gfs
  USE common_mpi
  USE common_mpi_gfs
  USE common_letkf
  USE efso_nml
  USE efso_tools
  USE letkf_obs
  USE letkf_tools

  IMPLICIT NONE
  REAL(r_size),ALLOCATABLE :: gues3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: gues2d(:,:)
  REAL(r_size),ALLOCATABLE :: fcst3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: fcst2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: fcer3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: fcer2d(:,:)
  REAL(r_size),ALLOCATABLE :: work3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: work2d(:,:)
  REAL(r_sngl),ALLOCATABLE :: work3dg(:,:,:,:)
  REAL(r_sngl),ALLOCATABLE :: work2dg(:,:,:)
  REAL(r_size),ALLOCATABLE :: uadf(:,:), vadf(:,:)
  REAL(r_size),ALLOCATABLE :: uada(:,:), vada(:,:)

  REAL(r_sngl) :: rtimer00,rtimer
  INTEGER :: n,ierr,ilev
  CHARACTER(8) :: stdoutf='NOUT-000'
  CHARACTER(4) :: fcstf='fc01'
  CHARACTER(9) :: analf='anal0.grd'
  CHARACTER(9) :: gmeanf='gmean.grd'
  CHARACTER(9) :: fmean0='fme00.grd'
  CHARACTER(9) :: fmean6='fme06.grd'
  CHARACTER(9) :: ameanf='amean.grd'
!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------
  CALL CPU_TIME(rtimer00)
  CALL read_namelist
  CALL initialize_mpi
!
  WRITE(stdoutf(6:8), '(I3.3)') myrank
  WRITE(6,'(3A,I3.3)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  OPEN(6,FILE=stdoutf)
  WRITE(6,'(A,I3.3,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf
!
  CALL set_common_gfs
  CALL set_common_mpi_gfs
!
  ALLOCATE(gues3d(nij1,nlev,nv3d))
  ALLOCATE(gues2d(nij1,nv2d))
  ALLOCATE(fcst3d(nij1,nlev,nbv,nv3d))
  ALLOCATE(fcst2d(nij1,nbv,nv2d))
  ALLOCATE(fcer3d(nij1,nlev,nv3d))
  ALLOCATE(fcer2d(nij1,nv2d))
  ALLOCATE(work3d(nij1,nlev,nv3d))
  ALLOCATE(work2d(nij1,nv2d))
  ALLOCATE(work3dg(nlon,nlat,nlev,nv3d))
  ALLOCATE(work2dg(nlon,nlat,nv2d))
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(INITIALIZE):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Read observation diagnostics
!-----------------------------------------------------------------------
  CALL set_efso_obs
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_OBS):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Read model data
!-----------------------------------------------------------------------
  !
  ! Forecast ensemble
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL read_ens_mpi(fcstf,nbv,fcst3d,fcst2d)
  !!! fcst3d,fcst2d: (xmean+X)^f_t [Eq.(6), Ota et al. 2013]
  !
  ! Forecast error at evaluation time
  !
  IF(myrank == 0) THEN
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',fmean0
    CALL read_grd4(fmean0,work3dg,work2dg,1)
  END IF
  CALL scatter_grd_mpi(0,work3dg,work2dg,fcer3d,fcer2d)
  IF(myrank == 0) THEN
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',fmean6
    CALL read_grd4(fmean6,work3dg,work2dg,1)
  END IF
  CALL scatter_grd_mpi(0,work3dg,work2dg,work3d,work2d)
  fcer3d(:,:,:) = 0.5_r_size * (fcer3d(:,:,:) + work3d(:,:,:))
  fcer2d(:,:) = 0.5_r_size * (fcer2d(:,:) + work2d(:,:))
  IF(myrank == 0) THEN
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',ameanf
    CALL read_grd4(ameanf,work3dg,work2dg,1)
  END IF
  CALL scatter_grd_mpi(0,work3dg,work2dg,work3d,work2d)
  fcer3d(:,:,:) = (fcer3d(:,:,:) - work3d(:,:,:)) / REAL(nbv-1,r_size)
  fcer2d(:,:) = (fcer2d(:,:) - work2d(:,:)) / REAL(nbv-1,r_size)
  !!! fcer3d,fcer2d: [1/2(K-1)](e^f_t+e^g_t) [Eq.(6), Ota et al. 2013]
  !
  ! Norm
  !
  CALL lnorm(fcst3d,fcst2d,fcer3d,fcer2d)
  !!! fcst3d,fcst2d: C^(1/2)*X^f_t
  !!! fcer3d,fcer2d: C^(1/2)*[1/2(K-1)](e^f_t+e^g_t)
  !
  ! Guess mean for full-level pressure computation
  !
  IF(myrank == 0) THEN
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',gmeanf
    CALL read_grd4(gmeanf,work3dg,work2dg,0)
  END IF
  CALL scatter_grd_mpi(0,work3dg,work2dg,gues3d,gues2d)
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_FORECAST):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Winds for advection
!-----------------------------------------------------------------------
  IF(ABS(locadv_rate) > TINY(locadv_rate)) THEN
    ALLOCATE(uadf(nij1,nlev))
    ALLOCATE(vadf(nij1,nlev))
    ALLOCATE(uada(nij1,nlev))
    ALLOCATE(vada(nij1,nlev))
    uadf(:,:) = work3d(:,:,iv3d_u)
    vadf(:,:) = work3d(:,:,iv3d_v)
    IF(myrank == 0) THEN
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',analf
      CALL read_grd4(analf,work3dg,work2dg,0)
    END IF
    CALL scatter_grd_mpi(0,work3dg,work2dg,work3d,work2d)
    uada(:,:) = work3d(:,:,iv3d_u)
    vada(:,:) = work3d(:,:,iv3d_v)
    CALL loc_advection(uada,vada,uadf,vadf) ! ADVECTION for FSO
    DEALLOCATE(uadf,vadf,uada,vada)
!
    CALL CPU_TIME(rtimer)
    WRITE(6,'(A,2F10.2)') '### TIMER(WIND_ADVECTION):',rtimer,rtimer-rtimer00
    rtimer00=rtimer
  END IF
  DEALLOCATE(work3d,work2d)
!-----------------------------------------------------------------------
! EFSO computation
!-----------------------------------------------------------------------
  CALL init_obsense()
!
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL das_efso(gues3d,gues2d,fcst3d,fcst2d,fcer3d,fcer2d)
  DEALLOCATE(gues3d,gues2d,fcst3d,fcst2d,fcer3d,fcer2d)
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(DAS_EFSO):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! EFSO output
!-----------------------------------------------------------------------
  IF(myrank == 0) CALL print_obsense()
  CALL destroy_obsense()
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(EFSO_OUTPUT):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi

  STOP
END PROGRAM efso
