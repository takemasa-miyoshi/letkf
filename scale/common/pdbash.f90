program pdbash
!=======================================================================
!
! [PURPOSE:] Execute parallel distributed bash script
!
!=======================================================================
  implicit none
  include 'mpif.h'

  integer, parameter :: maxlen = 6400
  integer :: nprocs, myrank, ierr, pos
  character(len=maxlen) :: cmd
  character(len=10) :: myranks

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

  call get_command(cmd)
  pos = index(cmd, ' ')
  cmd = cmd(pos+1:maxlen)
  pos = index(cmd, ' ')

  write (myranks, '(I10)') myrank

!     print *, 'bash ' // trim(cmd(1:pos-1)) // trim(myranks) // trim(cmd(pos:maxlen))
  call system('bash ' // trim(cmd(1:pos-1)) // trim(myranks) // trim(cmd(pos:maxlen)))

  call MPI_FINALIZE(ierr)

end program pdbash
