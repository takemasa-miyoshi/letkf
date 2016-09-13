PROGRAM mean 
!=======================================================================
!
! [PURPOSE:] Compute ensemble mean and spread
!
! [HISTORY:]
!   01/07/2013 Guo-Yuan Lien  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_gfs
  USE common_mpi_gfs
  USE common_letkf, only: nbv
!  USE letkf_tools

  IMPLICIT NONE
  REAL(r_size),ALLOCATABLE :: gues3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: gues2d(:,:,:)
  REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr
  CHARACTER(8) :: stdoutf='NOUT-000'
!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------
  CALL CPU_TIME(rtimer00)
  CALL initialize_mpi
!
  WRITE(stdoutf(6:8), '(I3.3)') myrank
  WRITE(6,'(3A,I3.3)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  OPEN(6,FILE=stdoutf)
  WRITE(6,'(A,I3.3,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf
!
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '  LOCAL ENSEMBLE TRANSFORM KALMAN FILTERING  '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '   LL      EEEEEE  TTTTTT  KK  KK  FFFFFF    '
  WRITE(6,'(A)') '   LL      EE        TT    KK KK   FF        '
  WRITE(6,'(A)') '   LL      EEEEE     TT    KKK     FFFFF     '
  WRITE(6,'(A)') '   LL      EE        TT    KK KK   FF        '
  WRITE(6,'(A)') '   LLLLLL  EEEEEE    TT    KK  KK  FF        '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '             WITHOUT LOCAL PATCH             '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '          Coded by Takemasa Miyoshi          '
  WRITE(6,'(A)') '  Based on Ott et al (2004) and Hunt (2005)  '
  WRITE(6,'(A)') '  Tested by Miyoshi and Yamane (2006)        '
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '              LETKF PARAMETERS               '
  WRITE(6,'(A)') ' ------------------------------------------- '
  WRITE(6,'(A,I15)')   '  nbv          :',nbv
  WRITE(6,'(A)') '============================================='
  CALL set_common_gfs
  CALL set_common_mpi_gfs
  ALLOCATE(gues3d(nij1,nlev,nbv,nv3d))
  ALLOCATE(gues2d(nij1,nbv,nv2d))
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(INITIALIZE):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! First guess ensemble
!-----------------------------------------------------------------------
  !
  ! READ GUES
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL read_ens_mpi('gues',nbv,gues3d,gues2d)
  !
  ! WRITE ENS MEAN and SPRD
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL write_ensmspr_mpi('gues',nbv,gues3d,gues2d)
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_GUES):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi

  STOP
END PROGRAM mean
