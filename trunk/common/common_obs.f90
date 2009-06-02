MODULE common_obs
!=======================================================================
!
! [PURPOSE:] General procedures related to observational data
!
! [HISTORY:]
!   10/14/2005 Takemasa MIYOSHI  created
!
!=======================================================================
  USE common

  IMPLICIT NONE
  PUBLIC

  INTEGER :: nobs
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331
  INTEGER,PARAMETER :: id_ps_obs=14593
  INTEGER,PARAMETER :: id_z_obs=2567
  INTEGER,PARAMETER :: id_s_obs=3332
  INTEGER,PARAMETER :: id_rain_obs=9999

CONTAINS
SUBROUTINE get_nobs(iunit)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: iunit
  REAL(r_sngl) :: wk(6)
  INTEGER :: ios
  INTEGER :: iu,iv,it,iq,irh,ips,iz,is

  nobs = 0
  iu = 0
  iv = 0
  it = 0
  iq = 0
  irh = 0
  ips = 0
  iz = 0
  is = 0
  DO
    READ(iunit,IOSTAT=ios) wk
    IF(ios /= 0) EXIT
    SELECT CASE(INT(wk(1)))
    CASE(id_u_obs)
      iu = iu + 1
    CASE(id_v_obs)
      iv = iv + 1
    CASE(id_t_obs)
      it = it + 1
    CASE(id_q_obs)
      iq = iq + 1
    CASE(id_rh_obs)
      irh = irh + 1
    CASE(id_ps_obs)
      ips = ips + 1
    CASE(id_z_obs)
      iz = iz + 1
    CASE(id_s_obs)
      is = is + 1
    END SELECT
    nobs = nobs + 1
  END DO
  WRITE(6,'(I10,A)') nobs,' OBSERVATIONS INPUT'
  WRITE(6,'(A12,I10)') '          U:',iu
  WRITE(6,'(A12,I10)') '          V:',iv
  WRITE(6,'(A12,I10)') '          T:',it
  WRITE(6,'(A12,I10)') '          Q:',iq
  WRITE(6,'(A12,I10)') '         RH:',irh
  WRITE(6,'(A12,I10)') '         PS:',ips
  WRITE(6,'(A12,I10)') '          Z:',iz
  WRITE(6,'(A12,I10)') '       SALT:',is
  REWIND(iunit)

  RETURN
END SUBROUTINE get_nobs

SUBROUTINE read_obs(iunit,elem,rlon,rlat,rlev,odat,oerr)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: iunit
  REAL(r_size),INTENT(OUT) :: elem(nobs) ! element number
  REAL(r_size),INTENT(OUT) :: rlon(nobs)
  REAL(r_size),INTENT(OUT) :: rlat(nobs)
  REAL(r_size),INTENT(OUT) :: rlev(nobs)
  REAL(r_size),INTENT(OUT) :: odat(nobs)
  REAL(r_size),INTENT(OUT) :: oerr(nobs)
  REAL(r_sngl) :: wk(6)
  INTEGER :: n

  n = 1
  DO
    READ(iunit) wk
    SELECT CASE(NINT(wk(1)))
    CASE(id_u_obs)
      wk(4) = wk(4) * 100.0d0 ! hPa -> Pa
    CASE(id_v_obs)
      wk(4) = wk(4) * 100.0d0 ! hPa -> Pa
    CASE(id_t_obs)
      wk(4) = wk(4) * 100.0d0 ! hPa -> Pa
    CASE(id_q_obs)
      wk(4) = wk(4) * 100.0d0 ! hPa -> Pa
    CASE(id_ps_obs)
      wk(5) = wk(5) * 100.0d0 ! hPa -> Pa
      wk(6) = wk(6) * 100.0d0 ! hPa -> Pa
    CASE(id_rh_obs)
      wk(4) = wk(4) * 100.0d0 ! hPa -> Pa
      wk(5) = wk(5) * 0.01d0 ! percent
      wk(6) = wk(6) * 0.01d0 ! percent
    CASE(id_z_obs)
      wk(4) = wk(4) * 100.0d0 ! hPa -> Pa
    END SELECT
    elem(n) = REAL(wk(1),r_size)
    rlon(n) = REAL(wk(2),r_size)
    rlat(n) = REAL(wk(3),r_size)
    rlev(n) = REAL(wk(4),r_size)
    odat(n) = REAL(wk(5),r_size)
    oerr(n) = REAL(wk(6),r_size)
    n = n + 1
    IF(n > nobs) EXIT
  END DO
  REWIND(iunit)

  RETURN
END SUBROUTINE read_obs

END MODULE common_obs
