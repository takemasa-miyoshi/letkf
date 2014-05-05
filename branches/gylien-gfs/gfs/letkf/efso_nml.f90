MODULE efso_nml
!=======================================================================
!
! [PURPOSE:] Read Namelist
!
! [HISTORY:]
!   04/12/2009 Yoichiro Ohta  created
!   04/27/2010 Yoichiro Ohta  Adaptive inflation
!   07/27/2011 Yoichiro Ohta  Modified for observation sensitivity
!   07/01/2013 Daisuke Hotta  Ported to GFS-LETKF system
!   12/19/2013 Guo-Yuan Lien  merged to GFS-LETKF main development
!
!=======================================================================
  USE common

  IMPLICIT NONE
  PUBLIC

  REAL(r_size),SAVE :: wmoist          ! Weight for moist total energy term
  REAL(r_size),SAVE :: eft             ! Evaluation FT for FSO
  REAL(r_size),SAVE :: locadv_rate     ! Localization advection rate
  REAL(r_size),SAVE :: tar_minlon      ! Target area (longitude)
  REAL(r_size),SAVE :: tar_maxlon      ! Target area (longitude)
  REAL(r_size),SAVE :: tar_minlat      ! Target area (latitude)
  REAL(r_size),SAVE :: tar_maxlat      ! Target area (latitude)
  INTEGER,SAVE :: tar_minlev           ! Target area (level)
  INTEGER,SAVE :: tar_maxlev           ! Target area (level)

CONTAINS
SUBROUTINE read_namelist
  IMPLICIT NONE
  CHARACTER(len=13) :: namelistfile='namelist.efso'
  INTEGER :: ierr
  !
  ! Namelist
  !
  NAMELIST /EFSOPRM/ wmoist,eft,locadv_rate, &
                   & tar_minlon,tar_maxlon,tar_minlat,tar_maxlat, &
                   & tar_minlev,tar_maxlev
  !
  ! Default Setting
  !
  wmoist=0.0d0; eft=0.0d0
  locadv_rate=0.0d0
  tar_minlon=0.0d0; tar_maxlon=360.0d0
  tar_minlat=-90.0d0; tar_maxlat=90.0d0
  tar_minlev=1; tar_maxlev=64
  !
  ! Read Namelist
  !
  OPEN(7,FILE=namelistfile,STATUS='old',FORM='formatted',iostat=ierr)
  READ(7,EFSOPRM)
  CLOSE(7)

  RETURN

END SUBROUTINE read_namelist
END MODULE efso_nml
