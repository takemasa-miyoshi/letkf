!!!!!  ==========================================================  !!!!!
!!!!!                 module  'physcons' description               !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!   this module contains some the most frequently used math and        !
!   physics constatns for gcm models.                                  !
!                                                                      !
!   references:                                                        !
!     as set in NMC handbook from Smithsonian tables.                  !
!                                                                      !
!   modification history:                                              !
!                                                                      !
!     1990-04-30  g and rd are made consistent with NWS usage          !
!     2001-10-22  g made consistent with SI usage                      !
!     2005-04-13  added molicular weights for gases          - y-t hou !
!                                                                      !
!   external modules referenced:                                       !
!                                                                      !
!       'module machine'                    in 'machine.f'             !
!                                                                      !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!

!========================================!
          module physcons                !
!........................................!
!
  use common, only: r_size
!
  implicit none
!
  public

!  --- ...  Math constants

  real(r_size),parameter:: con_pi     =3.1415926535897931 ! pi
  real(r_size),parameter:: con_sqrt2  =1.414214e+0        ! square root of 2
  real(r_size),parameter:: con_sqrt3  =1.732051e+0        ! square root of 3

!  --- ...  Geophysics/Astronomy constants

  real(r_size),parameter:: con_rerth  =6.3712e+6      ! radius of earth   (m)
  real(r_size),parameter:: con_g      =9.80665e+0     ! gravity           (m/s2)
  real(r_size),parameter:: con_omega  =7.2921e-5      ! ang vel of earth  (1/s)
  real(r_size),parameter:: con_p0     =1.01325e5      ! std atms pressure (pa)
! real(r_size),parameter:: con_solr   =1.3533e+3      ! solar constant    (W/m2)
! real(r_size),parameter:: con_solr   =1.36822e+3     ! solar constant    (W/m2)-aer(2001)
  real(r_size),parameter:: con_solr   =1.3660e+3      ! solar constant    (W/m2)-liu(2002)
! real(r_size),parameter:: con_solr   =1.36742732e+3  ! solar constant    (W/m2)-gfdl(1989)

!  --- ...  Thermodynamics constants

  real(r_size),parameter:: con_rgas   =8.314472       ! molar gas constant  (J/mol/K)
  real(r_size),parameter:: con_rd     =2.8705e+2      ! gas constant air    (J/kg/K)
  real(r_size),parameter:: con_rv     =4.6150e+2      ! gas constant H2O    (J/kg/K)
  real(r_size),parameter:: con_cp     =1.0046e+3      ! spec heat air @p    (J/kg/K)
  real(r_size),parameter:: con_cv     =7.1760e+2      ! spec heat air @v    (J/kg/K)
  real(r_size),parameter:: con_cvap   =1.8460e+3      ! spec heat H2O gas   (J/kg/K)
  real(r_size),parameter:: con_cliq   =4.1855e+3      ! spec heat H2O liq   (J/kg/K)
  real(r_size),parameter:: con_csol   =2.1060e+3      ! spec heat H2O ice   (J/kg/K)
  real(r_size),parameter:: con_hvap   =2.5000e+6      ! lat heat H2O cond   (J/kg)
  real(r_size),parameter:: con_hfus   =3.3358e+5      ! lat heat H2O fusion (J/kg)
  real(r_size),parameter:: con_psat   =6.1078e+2      ! pres at H2O 3pt     (Pa)  
  real(r_size),parameter:: con_t0c    =2.7315e+2      ! temp at 0C          (K)
  real(r_size),parameter:: con_ttp    =2.7316e+2      ! temp at H2O 3pt     (K)
  real(r_size),parameter:: con_jcal   =4.1855E+0      ! joules per calorie  ()

!  Secondary constants

  real(r_size),parameter:: con_rocp   =con_rd/con_cp
  real(r_size),parameter:: con_cpor   =con_cp/con_rd
  real(r_size),parameter:: con_rog    =con_rd/con_g
  real(r_size),parameter:: con_fvirt  =con_rv/con_rd-1.
  real(r_size),parameter:: con_eps    =con_rd/con_rv
  real(r_size),parameter:: con_epsm1  =con_rd/con_rv-1.
  real(r_size),parameter:: con_dldt   =con_cvap-con_cliq
  real(r_size),parameter:: con_xpona  =-con_dldt/con_rv
  real(r_size),parameter:: con_xponb  =-con_dldt/con_rv+con_hvap/(con_rv*con_ttp)

!  --- ...  Other Physics/Chemistry constants

  real(r_size),parameter:: con_plnk   =6.6260693e-34  ! planck constatn   (J/s) -nist(2002)
! real(r_size),parameter:: con_sbc    =5.6730e-8      ! stefan-boltzmann    (W/m2/K4)
  real(r_size),parameter:: con_sbc    =5.670400e-8    ! stefan-boltzmann    (W/m2/K4) -nist(2002)
! real(r_size),parameter:: con_avgd   =6.02214e23     ! avogadro constant   (1/mol) -aer
  real(r_size),parameter:: con_avgd   =6.0221415e23   ! avogadro constant   (1/mol) -nist(2002)
  real(r_size),parameter:: con_gasv   =22413.996e-6   ! vol of ideal gas at 273.15k, 101.325kpa (m3/mol) -nist(2002)
! real(r_size),parameter:: con_amd    =28.970         ! molecular wght of dry air (g/mol)
  real(r_size),parameter:: con_amd    =28.9644        ! molecular wght of dry air (g/mol)
  real(r_size),parameter:: con_amw    =18.0154        ! molecular wght of water vapor (g/mol)
  real(r_size),parameter:: con_amo3   =47.9982        ! molecular wght of o3  (g/mol)
! real(r_size),parameter:: con_amo3   =48.0           ! molecular wght of o3  (g/mol)
  real(r_size),parameter:: con_amco2  =44.011         ! molecular wght of co2 (g/mol)
  real(r_size),parameter:: con_amo2   =31.9999        ! molecular wght of o2  (g/mol)
  real(r_size),parameter:: con_amch4  =16.043         ! molecular wght of ch4 (g/mol)
  real(r_size),parameter:: con_amn2o  =44.013         ! molecular wght of n2o (g/mol)


!........................................!
      end module physcons                !
!========================================!
