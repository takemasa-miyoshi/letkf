!===============================================================================
module sscal_tools
!-------------------------------------------------------------------------------

  use sigio_module
  use sfcio_module
  use common, only: r_size
  use common_mpi

  implicit none

contains

!===============================================================================
subroutine sigio_data_add (datasig, datasig2)
!-------------------------------------------------------------------------------

  implicit none
  type(sigio_data), intent(inout) :: datasig
  type(sigio_data), intent(in)    :: datasig2

  datasig%hs  = datasig%hs  + datasig2%hs
  datasig%ps  = datasig%ps  + datasig2%ps
  datasig%t   = datasig%t   + datasig2%t
  datasig%d   = datasig%d   + datasig2%d
  datasig%z   = datasig%z   + datasig2%z
  datasig%q   = datasig%q   + datasig2%q
  datasig%xgr = datasig%xgr + datasig2%xgr
  datasig%xss = datasig%xss + datasig2%xss

!-------------------------------------------------------------------------------
end subroutine sigio_data_add
!===============================================================================

!===============================================================================
subroutine sigio_data_scalar_mul (datasig, scalar)
!-------------------------------------------------------------------------------

  implicit none
  type(sigio_data), intent(inout) :: datasig
  real(r_size), intent(in) :: scalar

  datasig%hs  = datasig%hs  * scalar
  datasig%ps  = datasig%ps  * scalar
  datasig%t   = datasig%t   * scalar
  datasig%d   = datasig%d   * scalar
  datasig%z   = datasig%z   * scalar
  datasig%q   = datasig%q   * scalar
  datasig%xgr = datasig%xgr * scalar
  datasig%xss = datasig%xss * scalar

!-------------------------------------------------------------------------------
end subroutine sigio_data_scalar_mul
!===============================================================================

!===============================================================================
subroutine sigio_data_mpi_reduce (datasig, datasig_res, OP, ROOT, MPI_COMM)
!-------------------------------------------------------------------------------

  implicit none
  type(sigio_data), intent(in)    :: datasig
  type(sigio_data), intent(inout) :: datasig_res
  integer, intent(in) :: OP
  integer, intent(in) :: ROOT
  integer, intent(in) :: MPI_COMM
  integer :: ierr

  call MPI_REDUCE(datasig%hs , datasig_res%hs , size(datasig%hs ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasig%ps , datasig_res%ps , size(datasig%ps ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasig%t  , datasig_res%t  , size(datasig%t  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasig%d  , datasig_res%d  , size(datasig%d  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasig%z  , datasig_res%z  , size(datasig%z  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasig%q  , datasig_res%q  , size(datasig%q  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasig%xgr, datasig_res%xgr, size(datasig%xgr), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasig%xss, datasig_res%xss, size(datasig%xss), MPI_REAL, OP, ROOT, MPI_COMM, ierr)

!-------------------------------------------------------------------------------
end subroutine sigio_data_mpi_reduce
!===============================================================================

!===============================================================================
subroutine sfcio_data_add (datasfc, datasfc2)
!-------------------------------------------------------------------------------

  implicit none
  type(sfcio_data), intent(inout) :: datasfc
  type(sfcio_data), intent(in)    :: datasfc2

  datasfc%tsea   = datasfc%tsea   + datasfc2%tsea
  datasfc%smc    = datasfc%smc    + datasfc2%smc
  datasfc%sheleg = datasfc%sheleg + datasfc2%sheleg
  datasfc%stc    = datasfc%stc    + datasfc2%stc
  datasfc%tg3    = datasfc%tg3    + datasfc2%tg3
  datasfc%zorl   = datasfc%zorl   + datasfc2%zorl
  datasfc%cv     = datasfc%cv     + datasfc2%cv
  datasfc%cvb    = datasfc%cvb    + datasfc2%cvb
  datasfc%cvt    = datasfc%cvt    + datasfc2%cvt
  datasfc%alvsf  = datasfc%alvsf  + datasfc2%alvsf
  datasfc%alvwf  = datasfc%alvwf  + datasfc2%alvwf
  datasfc%alnsf  = datasfc%alnsf  + datasfc2%alnsf
  datasfc%alnwf  = datasfc%alnwf  + datasfc2%alnwf
!  datasfc%slmsk  = datasfc%slmsk  + datasfc2%slmsk
  datasfc%vfrac  = datasfc%vfrac  + datasfc2%vfrac
  datasfc%canopy = datasfc%canopy + datasfc2%canopy
  datasfc%f10m   = datasfc%f10m   + datasfc2%f10m
  datasfc%t2m    = datasfc%t2m    + datasfc2%t2m
  datasfc%q2m    = datasfc%q2m    + datasfc2%q2m
!  datasfc%vtype  = datasfc%vtype  + datasfc2%vtype
!  datasfc%stype  = datasfc%stype  + datasfc2%stype
  datasfc%facsf  = datasfc%facsf  + datasfc2%facsf
  datasfc%facwf  = datasfc%facwf  + datasfc2%facwf
  datasfc%uustar = datasfc%uustar + datasfc2%uustar
  datasfc%ffmm   = datasfc%ffmm   + datasfc2%ffmm
  datasfc%ffhh   = datasfc%ffhh   + datasfc2%ffhh
  datasfc%hice   = datasfc%hice   + datasfc2%hice
  datasfc%fice   = datasfc%fice   + datasfc2%fice
  datasfc%tisfc  = datasfc%tisfc  + datasfc2%tisfc
  datasfc%tprcp  = datasfc%tprcp  + datasfc2%tprcp
!  datasfc%srflag = datasfc%srflag + datasfc2%srflag
  datasfc%snwdph = datasfc%snwdph + datasfc2%snwdph
  datasfc%slc    = datasfc%slc    + datasfc2%slc
  datasfc%shdmin = datasfc%shdmin + datasfc2%shdmin
  datasfc%shdmax = datasfc%shdmax + datasfc2%shdmax
  datasfc%slope  = datasfc%slope  + datasfc2%slope
  datasfc%snoalb = datasfc%snoalb + datasfc2%snoalb
  datasfc%orog   = datasfc%orog   + datasfc2%orog

!-------------------------------------------------------------------------------
end subroutine sfcio_data_add
!===============================================================================

!===============================================================================
subroutine sfcio_data_scalar_mul (datasfc, scalar)
!-------------------------------------------------------------------------------

  implicit none
  type(sfcio_data), intent(inout) :: datasfc
  real(r_size), intent(in) :: scalar

  datasfc%tsea   = datasfc%tsea   * scalar
  datasfc%smc    = datasfc%smc    * scalar
  datasfc%sheleg = datasfc%sheleg * scalar
  datasfc%stc    = datasfc%stc    * scalar
  datasfc%tg3    = datasfc%tg3    * scalar
  datasfc%zorl   = datasfc%zorl   * scalar
  datasfc%cv     = datasfc%cv     * scalar
  datasfc%cvb    = datasfc%cvb    * scalar
  datasfc%cvt    = datasfc%cvt    * scalar
  datasfc%alvsf  = datasfc%alvsf  * scalar
  datasfc%alvwf  = datasfc%alvwf  * scalar
  datasfc%alnsf  = datasfc%alnsf  * scalar
  datasfc%alnwf  = datasfc%alnwf  * scalar
!  datasfc%slmsk  = datasfc%slmsk  * scalar
  datasfc%vfrac  = datasfc%vfrac  * scalar
  datasfc%canopy = datasfc%canopy * scalar
  datasfc%f10m   = datasfc%f10m   * scalar
  datasfc%t2m    = datasfc%t2m    * scalar
  datasfc%q2m    = datasfc%q2m    * scalar
!  datasfc%vtype  = datasfc%vtype  * scalar
!  datasfc%stype  = datasfc%stype  * scalar
  datasfc%facsf  = datasfc%facsf  * scalar
  datasfc%facwf  = datasfc%facwf  * scalar
  datasfc%uustar = datasfc%uustar * scalar
  datasfc%ffmm   = datasfc%ffmm   * scalar
  datasfc%ffhh   = datasfc%ffhh   * scalar
  datasfc%hice   = datasfc%hice   * scalar
  datasfc%fice   = datasfc%fice   * scalar
  datasfc%tisfc  = datasfc%tisfc  * scalar
  datasfc%tprcp  = datasfc%tprcp  * scalar
!  datasfc%srflag = datasfc%srflag * scalar
  datasfc%snwdph = datasfc%snwdph * scalar
  datasfc%slc    = datasfc%slc    * scalar
  datasfc%shdmin = datasfc%shdmin * scalar
  datasfc%shdmax = datasfc%shdmax * scalar
  datasfc%slope  = datasfc%slope  * scalar
  datasfc%snoalb = datasfc%snoalb * scalar
  datasfc%orog   = datasfc%orog   * scalar

!-------------------------------------------------------------------------------
end subroutine sfcio_data_scalar_mul
!===============================================================================

!===============================================================================
subroutine sfcio_data_mpi_reduce (datasfc, datasfc_res, OP, ROOT, MPI_COMM)
!-------------------------------------------------------------------------------

  implicit none
  type(sfcio_data), intent(in)    :: datasfc
  type(sfcio_data), intent(inout) :: datasfc_res
  integer, intent(in) :: OP
  integer, intent(in) :: ROOT
  integer, intent(in) :: MPI_COMM
  integer :: ierr

  call MPI_REDUCE(datasfc%tsea   , datasfc_res%tsea   , size(datasfc%tsea   ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%smc    , datasfc_res%smc    , size(datasfc%smc    ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%sheleg , datasfc_res%sheleg , size(datasfc%sheleg ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%stc    , datasfc_res%stc    , size(datasfc%stc    ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%tg3    , datasfc_res%tg3    , size(datasfc%tg3    ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%zorl   , datasfc_res%zorl   , size(datasfc%zorl   ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%cv     , datasfc_res%cv     , size(datasfc%cv     ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%cvb    , datasfc_res%cvb    , size(datasfc%cvb    ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%cvt    , datasfc_res%cvt    , size(datasfc%cvt    ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%alvsf  , datasfc_res%alvsf  , size(datasfc%alvsf  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%alvwf  , datasfc_res%alvwf  , size(datasfc%alvwf  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%alnsf  , datasfc_res%alnsf  , size(datasfc%alnsf  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%alnwf  , datasfc_res%alnwf  , size(datasfc%alnwf  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
!  call MPI_REDUCE(datasfc%slmsk  , datasfc_res%slmsk  , size(datasfc%slmsk  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%vfrac  , datasfc_res%vfrac  , size(datasfc%vfrac  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%canopy , datasfc_res%canopy , size(datasfc%canopy ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%f10m   , datasfc_res%f10m   , size(datasfc%f10m   ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%t2m    , datasfc_res%t2m    , size(datasfc%t2m    ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%q2m    , datasfc_res%q2m    , size(datasfc%q2m    ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
!  call MPI_REDUCE(datasfc%vtype  , datasfc_res%vtype  , size(datasfc%vtype  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
!  call MPI_REDUCE(datasfc%stype  , datasfc_res%stype  , size(datasfc%stype  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%facsf  , datasfc_res%facsf  , size(datasfc%facsf  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%facwf  , datasfc_res%facwf  , size(datasfc%facwf  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%uustar , datasfc_res%uustar , size(datasfc%uustar ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%ffmm   , datasfc_res%ffmm   , size(datasfc%ffmm   ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%ffhh   , datasfc_res%ffhh   , size(datasfc%ffhh   ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%hice   , datasfc_res%hice   , size(datasfc%hice   ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%fice   , datasfc_res%fice   , size(datasfc%fice   ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%tisfc  , datasfc_res%tisfc  , size(datasfc%tisfc  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%tprcp  , datasfc_res%tprcp  , size(datasfc%tprcp  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
!  call MPI_REDUCE(datasfc%srflag , datasfc_res%srflag , size(datasfc%srflag ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%snwdph , datasfc_res%snwdph , size(datasfc%snwdph ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%slc    , datasfc_res%slc    , size(datasfc%slc    ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%shdmin , datasfc_res%shdmin , size(datasfc%shdmin ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%shdmax , datasfc_res%shdmax , size(datasfc%shdmax ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%slope  , datasfc_res%slope  , size(datasfc%slope  ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%snoalb , datasfc_res%snoalb , size(datasfc%snoalb ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)
  call MPI_REDUCE(datasfc%orog   , datasfc_res%orog   , size(datasfc%orog   ), MPI_REAL, OP, ROOT, MPI_COMM, ierr)

!-------------------------------------------------------------------------------
end subroutine sfcio_data_mpi_reduce
!===============================================================================

end module sscal_tools
!===============================================================================
