MODULE common_mpi
!=======================================================================
!
! [PURPOSE:] General MPI procedures
!
! [HISTORY:]
!   09/06/2005 Takemasa MIYOSHI  created
!
!=======================================================================
  USE common, only: r_size, r_dble, r_sngl
  IMPLICIT NONE
  PUBLIC
  INCLUDE 'mpif.h'

  INTEGER,SAVE :: nprocs
  INTEGER,SAVE :: myrank
  INTEGER,SAVE :: MPI_r_size

CONTAINS
SUBROUTINE initialize_mpi
  IMPLICIT NONE
  INTEGER :: ierr
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
  WRITE(6,'(A,I6.6,A,I6.6)') 'Hello from MYRANK ',myrank,'/',nprocs-1
  IF(r_size == r_dble) THEN
    MPI_r_size = MPI_DOUBLE_PRECISION
  ELSE IF(r_size == r_sngl) THEN
    MPI_r_size = MPI_REAL
  END IF

  RETURN
END SUBROUTINE initialize_mpi

SUBROUTINE finalize_mpi
  IMPLICIT NONE
  INTEGER :: ierr
  CALL MPI_FINALIZE(ierr)

  RETURN
END SUBROUTINE finalize_mpi

END MODULE common_mpi
