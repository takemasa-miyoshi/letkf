MODULE common_mpi
!=======================================================================
!
! [PURPOSE:] General MPI procedures
!
! [HISTORY:]
!   09/06/2005 Takemasa MIYOSHI  created
!
!=======================================================================
  IMPLICIT NONE
  PUBLIC
  INCLUDE 'mpif.h'

  INTEGER,SAVE :: nprocs
  INTEGER,SAVE :: myrank

CONTAINS
SUBROUTINE initialize_mpi
  IMPLICIT NONE
  INTEGER :: ierr
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
  WRITE(6,'(A,I3.3,A,I3.3)') 'Hello from MYRANK ',myrank,'/',nprocs-1

  RETURN
END SUBROUTINE initialize_mpi

SUBROUTINE finalize_mpi
  IMPLICIT NONE
  INTEGER :: ierr
  CALL MPI_FINALIZE(ierr)

  RETURN
END SUBROUTINE finalize_mpi

END MODULE common_mpi
