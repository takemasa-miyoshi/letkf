!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Compute ensemble mean from the GFS sig/sfc files
!
!  September 2013, Created Guo-Yuan Lien, University of Maryland
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program gfsmean

  use sigio_module
  use sfcio_module
  use common, only: r_size
  use common_mpi
  use common_letkf, only: nbv
  use sscal_tools

  implicit none

!-------------------------------------------------------------------------------

  character(len=6) :: siginfile = 'sigNNN'
  character(len=6) :: sfcinfile = 'sfcNNN'
  character(len=6) :: sigoutfile = 'sig_me'
  character(len=6) :: sfcoutfile = 'sfc_me'

  type(sigio_head) :: headsig
  type(sfcio_head) :: headsfc
  type(sigio_data) :: datasig, tmpdatasig
  type(sfcio_data) :: datasfc, tmpdatasfc

  character(len=8) :: stdoutf='NOUT-000'

  integer :: l, im, n, iret, ierr
  integer :: MPI_C, MPI_G, MPI_G_WORLD
  integer :: useranks(nbv)
  real(r_size) :: rtimer00, rtimer
  logical :: data_allocated

!-------------------------------------------------------------------------------
! Initialize
!-------------------------------------------------------------------------------

  call CPU_TIME(rtimer00)
  call initialize_mpi

  if (nprocs > nbv) then
    do n = 1, nbv
      useranks(n) = n-1
    end do
    call MPI_COMM_GROUP(MPI_COMM_WORLD, MPI_G_WORLD, ierr)
    call MPI_GROUP_INCL(MPI_G_WORLD, nbv, useranks, MPI_G, ierr)
    call MPI_COMM_CREATE(MPI_COMM_WORLD, MPI_G, MPI_C, ierr)
  else
    MPI_C = MPI_COMM_WORLD
  end if

  write (stdoutf(6:8), '(I3.3)') myrank
  write (6, '(3A,I3.3)') 'STDOUT goes to ', stdoutf, ' for MYRANK ', myrank
  open (6, file=stdoutf)
  write (6, '(A,I3.3,2A)') 'MYRANK=', myrank, ', STDOUTF=', stdoutf

  data_allocated = .false.

  if (myrank+1 > nbv) then
    write (6, '(A)') 'THIS PROCESSOR IS NOT USED.'
  else

    call CPU_TIME(rtimer)
    write (6, '(A,2F10.2)') '### TIMER(INITIALIZE):', rtimer, rtimer-rtimer00
    rtimer00 = rtimer

!-------------------------------------------------------------------------------
! Read sig/sfc files and do local summation
!-------------------------------------------------------------------------------

    l=0
    do
      im = myrank+1 + nprocs * l
      if (im > nbv) exit
      write (siginfile(4:6), '(I3.3)') im
      write (sfcinfile(4:6), '(I3.3)') im
      write (6, '(A,I3.3,4A)') 'MYRANK ', myrank, ' is reading a file ', siginfile, '/', sfcinfile

      if (l == 0) then
        call sigio_srohdc(21, siginfile, headsig, tmpdatasig, iret)
        if (iret /= 0) then
          write (0, '(2A)') '[Error] Open and read sigma file: ', siginfile
          stop
        endif
        call sfcio_srohdc(22, sfcinfile, headsfc, tmpdatasfc, iret)
        if (iret /= 0) then
          write (0, '(2A)') '[Error] Open and read sigma file: ', siginfile
          stop
        endif
      else
        call sigio_srohdc(21, siginfile, headsig, datasig, iret)
        if (iret /= 0) then
          write (0, '(2A)') '[Error] Open and read sigma file: ', siginfile
          stop
        endif
        call sfcio_srohdc(22, sfcinfile, headsfc, datasfc, iret)
        if (iret /= 0) then
          write (0, '(2A)') '[Error] Open and read sigma file: ', siginfile
          stop
        endif
        call sigio_data_add(tmpdatasig, datasig)
        call sfcio_data_add(tmpdatasfc, datasfc)
        data_allocated = .true.
      end if

      l = l + 1
    end do

    call CPU_TIME(rtimer)
    write (6, '(A,2F10.2)') '### TIMER(READ FILES):', rtimer, rtimer-rtimer00
    rtimer00 = rtimer

!-------------------------------------------------------------------------------
! MPI_REDUCE to compute the grand summation, then average
!-------------------------------------------------------------------------------

    if (myrank == 0) then
      if (.not. data_allocated) then
        call sigio_aldata(headsig, datasig, iret)
        call sfcio_aldata(headsfc, datasfc, iret)
        datasfc%slmsk  = tmpdatasfc%slmsk
        datasfc%vtype  = tmpdatasfc%vtype
        datasfc%stype  = tmpdatasfc%stype
        datasfc%srflag = tmpdatasfc%srflag
      end if
    end if

    call MPI_BARRIER(MPI_C, ierr)
    call sigio_data_mpi_reduce(tmpdatasig, datasig, MPI_SUM, 0, MPI_C)
    call MPI_BARRIER(MPI_C, ierr)
    call sfcio_data_mpi_reduce(tmpdatasfc, datasfc, MPI_SUM, 0, MPI_C)

    if (myrank == 0) then
      call sigio_data_scalar_mul(datasig, 1.0d0/real(nbv,r_size))
      call sfcio_data_scalar_mul(datasfc, 1.0d0/real(nbv,r_size))

      write (6, '(A,I3.3,4A)') 'MYRANK ', myrank, ' is writing a file ', sigoutfile, '/', sfcoutfile
      call sigio_swohdc(21, sigoutfile, headsig, datasig, iret)
      if (iret /= 0) then
        write (0, '(2A)') '[Error] Open and write sigma file: ', sigoutfile
        stop
      endif
      call sfcio_swohdc(22, sfcoutfile, headsfc, datasfc, iret)
      if (iret /= 0) then
        write (0, '(2A)') '[Error] Open and write sigma file: ', sfcoutfile
        stop
      endif
    end if

    call CPU_TIME(rtimer)
    write (6, '(A,2F10.2)') '### TIMER(MPI_REDUCE):', rtimer, rtimer-rtimer00
    rtimer00 = rtimer

  end if ! [ myrank+1 > nbv ]

!-------------------------------------------------------------------------------
! Finalize
!-------------------------------------------------------------------------------

  call finalize_mpi

!-------------------------------------------------------------------------------

end program gfsmean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
