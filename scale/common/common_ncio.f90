module common_ncio
!=======================================================================
!
! [PURPOSE:] NetCDF I/O
!
! [HISTORY:]
!   07/24/2014 Guo-Yuan Lien  created
!
!=======================================================================
  use netcdf
  use common, only: r_size, r_dble, r_sngl

  implicit none
  public

  interface ncio_read
    module procedure ncio_read_1d_r4
    module procedure ncio_read_2d_r4
    module procedure ncio_read_3d_r4
    module procedure ncio_read_1d_r8
    module procedure ncio_read_2d_r8
    module procedure ncio_read_3d_r8
  end interface ncio_read

  interface ncio_write
    module procedure ncio_write_1d_r4
    module procedure ncio_write_2d_r4
    module procedure ncio_write_3d_r4
    module procedure ncio_write_1d_r8
    module procedure ncio_write_2d_r8
    module procedure ncio_write_3d_r8
  end interface ncio_write

contains

!-----------------------------------------------------------------------
! Check the status
!-----------------------------------------------------------------------
subroutine ncio_check(status)
  implicit none
  integer, intent(in) :: status

  if (status /= nf90_noerr) then
    write (6,*) trim(nf90_strerror(status))
    stop 10
  end if
end subroutine ncio_check
!-----------------------------------------------------------------------
! Open a netcdf file
!-----------------------------------------------------------------------
subroutine ncio_open(filename, mode, ncid)
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: mode
  integer, intent(out) :: ncid

  call ncio_check(nf90_open(filename, mode, ncid))
end subroutine ncio_open
!-----------------------------------------------------------------------
! Close a netcdf file
!-----------------------------------------------------------------------
subroutine ncio_close(ncid)
  implicit none
  integer, intent(in) :: ncid

  call ncio_check(nf90_close(ncid))
end subroutine ncio_close
!-----------------------------------------------------------------------
! Read netcdf dimension
!-----------------------------------------------------------------------
subroutine ncio_read_dim(ncid, dimname, dimlen)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: dimname
  integer, intent(out) :: dimlen
  integer :: dimid

  call ncio_check(nf90_inq_dimid(ncid, dimname, dimid))
  call ncio_check(nf90_inquire_dimension(ncid, dimid, len=dimlen))
end subroutine ncio_read_dim
!-----------------------------------------------------------------------
! Read netcdf integer global attribute
!-----------------------------------------------------------------------
!subroutine ncio_read_gattr_i(ncid, attrname, attr)
!  implicit none
!  integer, intent(in) :: ncid
!  character(len=*), intent(in) :: attrname
!  integer, intent(out) :: attr

!  call ncio_check(nf90_get_att(ncid, nf90_global, attrname, attr))
!end subroutine ncio_read_gattr_i
!-----------------------------------------------------------------------
! Read netcdf single-precision global attribute
!-----------------------------------------------------------------------
!subroutine ncio_read_gattr_r4(ncid, attrname, attr)
!  implicit none
!  integer, intent(in) :: ncid
!  character(len=*), intent(in)  :: attrname
!  real(r_sngl), intent(out) :: attr

!  call ncio_check(nf90_get_att(ncid, nf90_global, attrname, attr))
!end subroutine ncio_read_gattr_r4
!-----------------------------------------------------------------------
! Read netcdf single-precision 1-D variable
!-----------------------------------------------------------------------
subroutine ncio_read_1d_r4(ncid, varname, dim1, t, var)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: dim1, t
  real(r_sngl), intent(out) :: var(dim1)
  integer :: varid

  call ncio_check(nf90_inq_varid(ncid, varname, varid))
  call ncio_check(nf90_get_var(ncid, varid, var,   &
                               start = (/ 1, t /), &
                               count = (/ dim1, 1 /)))
end subroutine ncio_read_1d_r4
!-----------------------------------------------------------------------
! Read netcdf single-precision 2-D variable
!-----------------------------------------------------------------------
subroutine ncio_read_2d_r4(ncid, varname, dim1, dim2, t, var)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: dim1, dim2, t
  real(r_sngl), intent(out) :: var(dim1,dim2)
  integer :: varid

  call ncio_check(nf90_inq_varid(ncid, varname, varid))
  call ncio_check(nf90_get_var(ncid, varid, var,      &
                               start = (/ 1, 1, t /), &
                               count = (/ dim1, dim2, 1 /)))
end subroutine ncio_read_2d_r4
!-----------------------------------------------------------------------
! Read netcdf single-precision 3-D variable
!-----------------------------------------------------------------------
subroutine ncio_read_3d_r4(ncid, varname, dim1, dim2, dim3, t, var)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: dim1, dim2, dim3, t
  real(r_sngl), intent(out) :: var(dim1,dim2,dim3)
  integer :: varid

  call ncio_check(nf90_inq_varid(ncid, varname, varid))
  call ncio_check(nf90_get_var(ncid, varid, var,         &
                               start = (/ 1, 1, 1, t /), &
                               count = (/ dim1, dim2, dim3, 1 /)))
end subroutine ncio_read_3d_r4
!-----------------------------------------------------------------------
! Read netcdf double-precision 1-D variable
!-----------------------------------------------------------------------
subroutine ncio_read_1d_r8(ncid, varname, dim1, t, var)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: dim1, t
  real(r_dble), intent(out) :: var(dim1)
  integer :: varid

  call ncio_check(nf90_inq_varid(ncid, varname, varid))
  call ncio_check(nf90_get_var(ncid, varid, var,   &
                               start = (/ 1, t /), &
                               count = (/ dim1, 1 /)))
end subroutine ncio_read_1d_r8
!-----------------------------------------------------------------------
! Read netcdf double-precision 2-D variable
!-----------------------------------------------------------------------
subroutine ncio_read_2d_r8(ncid, varname, dim1, dim2, t, var)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: dim1, dim2, t
  real(r_dble), intent(out) :: var(dim1,dim2)
  integer :: varid

  call ncio_check(nf90_inq_varid(ncid, varname, varid))
  call ncio_check(nf90_get_var(ncid, varid, var,      &
                               start = (/ 1, 1, t /), &
                               count = (/ dim1, dim2, 1 /)))
end subroutine ncio_read_2d_r8
!-----------------------------------------------------------------------
! Read netcdf double-precision 3-D variable
!-----------------------------------------------------------------------
subroutine ncio_read_3d_r8(ncid, varname, dim1, dim2, dim3, t, var)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: dim1, dim2, dim3, t
  real(r_dble), intent(out) :: var(dim1,dim2,dim3)
  integer :: varid

  call ncio_check(nf90_inq_varid(ncid, varname, varid))
  call ncio_check(nf90_get_var(ncid, varid, var,         &
                               start = (/ 1, 1, 1, t /), &
                               count = (/ dim1, dim2, dim3, 1 /)))
end subroutine ncio_read_3d_r8
!-----------------------------------------------------------------------
! Read netcdf text variable
!-----------------------------------------------------------------------
subroutine ncio_read_text(ncid, varname, length, t, var)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: length, t
  character(len=length), intent(out) :: var
  integer :: varid

  call ncio_check(nf90_inq_varid(ncid, varname, varid))
  call ncio_check(nf90_get_var(ncid, varid, var,   &
                               start = (/ 1, t /), &
                               count = (/ length, 1 /)))
end subroutine ncio_read_text
!-----------------------------------------------------------------------
! Write netcdf single-precision 1-D variable
!-----------------------------------------------------------------------
subroutine ncio_write_1d_r4(ncid, varname, dim1, t, var)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: dim1, t
  real(r_sngl), intent(in) :: var(dim1)
  integer :: varid

  call ncio_check(nf90_inq_varid(ncid, varname, varid))
  call ncio_check(nf90_put_var(ncid, varid, var,   &
                               start = (/ 1, t /), &
                               count = (/ dim1, 1 /)))
end subroutine ncio_write_1d_r4
!-----------------------------------------------------------------------
! Write netcdf single-precision 2-D variable
!-----------------------------------------------------------------------
subroutine ncio_write_2d_r4(ncid, varname, dim1, dim2, t, var)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: dim1, dim2, t
  real(r_sngl), intent(in) :: var(dim1,dim2)
  integer :: varid

  call ncio_check(nf90_inq_varid(ncid, varname, varid))
  call ncio_check(nf90_put_var(ncid, varid, var,      &
                               start = (/ 1, 1, t /), &
                               count = (/ dim1, dim2, 1 /)))
end subroutine ncio_write_2d_r4
!-----------------------------------------------------------------------
! Write netcdf single-precision 3-D variable
!-----------------------------------------------------------------------
subroutine ncio_write_3d_r4(ncid, varname, dim1, dim2, dim3, t, var)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: dim1, dim2, dim3, t
  real(r_sngl), intent(in) :: var(dim1,dim2,dim3)
  integer :: varid

  call ncio_check(nf90_inq_varid(ncid, varname, varid))
  call ncio_check(nf90_put_var(ncid, varid, var,         &
                               start = (/ 1, 1, 1, t /), &
                               count = (/ dim1, dim2, dim3, 1 /)))
end subroutine ncio_write_3d_r4
!-----------------------------------------------------------------------
! Write netcdf double-precision 1-D variable
!-----------------------------------------------------------------------
subroutine ncio_write_1d_r8(ncid, varname, dim1, t, var)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: dim1, t
  real(r_dble), intent(in) :: var(dim1)
  integer :: varid

  call ncio_check(nf90_inq_varid(ncid, varname, varid))
  call ncio_check(nf90_put_var(ncid, varid, var,   &
                               start = (/ 1, t /), &
                               count = (/ dim1, 1 /)))
end subroutine ncio_write_1d_r8
!-----------------------------------------------------------------------
! Write netcdf double-precision 2-D variable
!-----------------------------------------------------------------------
subroutine ncio_write_2d_r8(ncid, varname, dim1, dim2, t, var)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: dim1, dim2, t
  real(r_dble), intent(in) :: var(dim1,dim2)
  integer :: varid

  call ncio_check(nf90_inq_varid(ncid, varname, varid))
  call ncio_check(nf90_put_var(ncid, varid, var,      &
                               start = (/ 1, 1, t /), &
                               count = (/ dim1, dim2, 1 /)))
end subroutine ncio_write_2d_r8
!-----------------------------------------------------------------------
! Write netcdf double-precision 3-D variable
!-----------------------------------------------------------------------
subroutine ncio_write_3d_r8(ncid, varname, dim1, dim2, dim3, t, var)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: dim1, dim2, dim3, t
  real(r_dble), intent(in) :: var(dim1,dim2,dim3)
  integer :: varid

  call ncio_check(nf90_inq_varid(ncid, varname, varid))
  call ncio_check(nf90_put_var(ncid, varid, var,         &
                               start = (/ 1, 1, 1, t /), &
                               count = (/ dim1, dim2, dim3, 1 /)))
end subroutine ncio_write_3d_r8
!-----------------------------------------------------------------------
! Write netcdf text variable
!-----------------------------------------------------------------------
subroutine ncio_write_text (ncid, varname, length, t, var)
  implicit none
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: length, t
  character(len=length), intent(in) :: var
  integer :: varid

  call ncio_check(nf90_inq_varid(ncid, varname, varid))
  call ncio_check(nf90_put_var(ncid, varid, var,   &
                               start = (/ 1, t /), &
                               count = (/ length, 1 /)))
end subroutine ncio_write_text

end module common_ncio
