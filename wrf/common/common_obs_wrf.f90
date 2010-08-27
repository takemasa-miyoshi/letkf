MODULE common_obs_wrf
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_wrf

  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331
  INTEGER,PARAMETER :: id_ps_obs=14593
  INTEGER,PARAMETER :: id_rain_obs=9999

CONTAINS
!-----------------------------------------------------------------------
! Transformation from model variables to an observation
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY(elm,ri,rj,rk,v3d,v2d,yobs)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rk
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: yobs
  REAL(r_size) :: rh(nlon,nlat,nlev)
  REAL(r_size) :: tg,qg
  INTEGER :: i,j,k
  INTEGER :: is,ie,js,je,ks,ke
  ie = CEILING( ri )
  is = ie-1
  je = CEILING( rj )
  js = je-1
  ke = CEILING( rk )
  ks = ke-1

  SELECT CASE (NINT(elm))
  CASE(id_u_obs)  ! U
    CALL itpl_3d(v3d(:,:,:,iv3d_u),ri,rj,rk,yobs)
  CASE(id_v_obs)  ! V
    CALL itpl_3d(v3d(:,:,:,iv3d_v),ri,rj,rk,yobs)
  CASE(id_t_obs)  ! T
    CALL itpl_3d(v3d(:,:,:,iv3d_t),ri,rj,rk,yobs)
  CASE(id_q_obs)  ! Q
    CALL itpl_3d(v3d(:,:,:,iv3d_qv),ri,rj,rk,yobs)
  CASE(id_ps_obs) ! PS
    CALL itpl_2d(v2d(:,:,iv2d_t2),ri,rj,tg)
    CALL itpl_2d(v2d(:,:,iv2d_q2),ri,rj,qg)
    CALL itpl_2d(v2d(:,:,iv2d_ps),ri,rj,yobs)
    CALL prsadj(yobs,rk,tg,qg)
  CASE(id_rain_obs) ! RAIN
    CALL itpl_2d(v2d(:,:,iv2d_rain),ri,rj,yobs)
  CASE(id_rh_obs) ! RH
    DO k=ks,ke
      DO j=js,je
        IF(ie <= nlon ) THEN
          CALL calc_rh(v3d(is,j,k,iv3d_t),v3d(is,j,k,iv3d_qv),&
            & v3d(is,j,k,iv3d_p),rh(is,j,k))
          CALL calc_rh(v3d(ie,j,k,iv3d_t),v3d(ie,j,k,iv3d_qv),&
            & v3d(ie,j,k,iv3d_p),rh(ie,j,k))
        ELSE
          CALL calc_rh(v3d(is,j,k,iv3d_t),v3d(is,j,k,iv3d_qv),&
            & v3d(is,j,k,iv3d_p),rh(is,j,k))
          CALL calc_rh(v3d( 1,j,k,iv3d_t),v3d( 1,j,k,iv3d_qv),&
            & v3d( 1,j,k,iv3d_p),rh( 1,j,k))
        END IF
      END DO
    END DO
    CALL itpl_3d(rh,ri,rj,rk,yobs)
  END SELECT

  RETURN
END SUBROUTINE Trans_XtoY
!-----------------------------------------------------------------------
! Compute relative humidity (RH)
!-----------------------------------------------------------------------
SUBROUTINE calc_rh(t,q,p,rh)
  IMPLICIT NONE
  REAL(r_size),PARAMETER :: t0=273.15d0
  REAL(r_size),PARAMETER :: e0c=6.11d0
  REAL(r_size),PARAMETER :: al=17.3d0
  REAL(r_size),PARAMETER :: bl=237.3d0
  REAL(r_size),PARAMETER :: e0i=6.1121d0
  REAL(r_size),PARAMETER :: ai=22.587d0
  REAL(r_size),PARAMETER :: bi=273.86d0
  REAL(r_size),INTENT(IN) :: t,q,p
  REAL(r_size),INTENT(OUT) :: rh
  REAL(r_size) :: e,es,tc

  e = q * p * 0.01d0 / (0.378d0 * q + 0.622d0)

  tc = t-t0
  IF(tc >= 0.0d0) THEN
    es = e0c * exp(al*tc/(bl+tc))
  ELSE IF(tc <= -15.d0) THEN
    es = e0i * exp(ai*tc/(bi+tc))
  ELSE
    es = e0c * exp(al*tc/(bl+tc)) * (15.0d0+tc)/15.0d0 &
       + e0i * exp(ai*tc/(bi+tc)) * (-tc) / 15.0d0
  END IF

  rh = e/es

  RETURN
END SUBROUTINE calc_rh
!-----------------------------------------------------------------------
! Pressure adjustment for a different height level
!-----------------------------------------------------------------------
SUBROUTINE prsadj(p,dz,t,q)
  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: p
  REAL(r_size),INTENT(IN) :: dz ! height difference (target - original) [m]
  REAL(r_size),INTENT(IN) :: t  ! temperature [K] at original level
  REAL(r_size),INTENT(IN) :: q  ! humidity [kg/kg] at original level
  REAL(r_size),PARAMETER :: gamma=5.0d-3 ! lapse rate [K/m]
  REAL(r_size) :: tv

  IF(dz /= 0) THEN
    tv = t * (1.0d0 + 0.608d0 * q)
    p = p * ((-gamma*dz+tv)/tv)**(gg/(gamma*rd)) !tv is at original level
!    p = p * (tv/(tv+gamma*dz))**(gg/(gamma*rd)) !tv is at target level
  END IF

  RETURN
END SUBROUTINE prsadj
!-----------------------------------------------------------------------
! Horizontal coordinate conversion
!-----------------------------------------------------------------------
SUBROUTINE ll2ij(elem,rlon,rlat,ri,rj)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: elem
  REAL(r_size),INTENT(IN) :: rlon(1)
  REAL(r_size),INTENT(IN) :: rlat(1)
  REAL(r_size),INTENT(OUT) :: ri(1)
  REAL(r_size),INTENT(OUT) :: rj(1)
  INTEGER :: i,j
  !
  ! rlon, rlat --> ri, rj
  !
  SELECT CASE(NINT(elem))
  CASE(id_u_obs)
    i=nlon
    j=nlat-1
    CALL com_pos2ij(2,i,j,lonu(1:i,1:j),latu(1:i,1:j),1,rlon,rlat,ri,rj)
  CASE(id_v_obs)
    i=nlon-1
    j=nlat
    CALL com_pos2ij(2,i,j,lonv(1:i,1:j),latv(1:i,1:j),1,rlon,rlat,ri,rj)
  CASE DEFAULT
    i=nlon-1
    j=nlat-1
    CALL com_pos2ij(2,i,j,lon(1:i,1:j),lat(1:i,1:j),1,rlon,rlat,ri,rj)
  END SELECT

  RETURN
END SUBROUTINE ll2ij
!-----------------------------------------------------------------------
! Vertical coordinate conversion
!-----------------------------------------------------------------------
SUBROUTINE p2k(p_full,elem,ri,rj,rlev,rk)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: p_full(nlon,nlat,nlev)
  REAL(r_size),INTENT(IN) :: elem
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rlev ! pressure levels
  REAL(r_size),INTENT(OUT) :: rk
  REAL(r_size) :: ak
  REAL(r_size) :: lnps(nlon,nlat)
  REAL(r_size) :: plev(nlev)
  INTEGER :: i,j,k,n
  !
  ! rlev --> rk
  !
  IF(NINT(elem) == id_ps_obs) THEN ! surface pressure observation
    rk = 0.0d0
  ELSE
    !
    ! horizontal interpolation
    !
    i = CEILING(ri)
    j = CEILING(rj)
    IF(p_full(i,j,1) < 0.1) THEN
      rk = -1.0d0
      RETURN
    END IF
    DO k=1,nlev-1
      lnps(i-1:i,j-1:j) = LOG(p_full(i-1:i,j-1:j,k))
      CALL itpl_2d(lnps,ri,rj,plev(k))
    END DO
    !
    ! Log pressure
    !
    rk = LOG(rlev)
    !
    ! find rk
    !
    DO k=2,nlev-1
      IF(plev(k) < rk) EXIT ! assuming descending order of plev
    END DO
    ak = (rk - plev(k-1)) / (plev(k) - plev(k-1))
    rk = REAL(k-1,r_size) + ak
  END IF

  RETURN
END SUBROUTINE p2k
!-----------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------
SUBROUTINE itpl_2d(var,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlon,nlat)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  IF(i <= nlon) THEN
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(i  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(i  ,j  ) *    ai  *    aj
  ELSE
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(1  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(1  ,j  ) *    ai  *    aj
  END IF

  RETURN
END SUBROUTINE itpl_2d

SUBROUTINE itpl_3d(var,ri,rj,rk,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlon,nlat,nlev)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rk
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj,ak
  INTEGER :: i,j,k

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)
  k = CEILING(rk)
  ak = rk - REAL(k-1,r_size)

  IF(i <= nlon) THEN
    var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(i  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
       & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
       & + var(i  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
       & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
       & + var(i  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
       & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
       & + var(i  ,j  ,k  ) *    ai  *    aj  *    ak
  ELSE
    var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(1  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
       & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
       & + var(1  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
       & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
       & + var(1  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
       & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
       & + var(1  ,j  ,k  ) *    ai  *    aj  *    ak
  END IF

  RETURN
END SUBROUTINE itpl_3d
!-----------------------------------------------------------------------
! Monitor departure
!-----------------------------------------------------------------------
SUBROUTINE monit_dep(nn,elm,dep,qc)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elm(nn)
  REAL(r_size),INTENT(IN) :: dep(nn)
  INTEGER,INTENT(IN) :: qc(nn)
  REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_q,rmse_ps,rmse_rh
  REAL(r_size) :: bias_u,bias_v,bias_t,bias_q,bias_ps,bias_rh
  INTEGER :: n,iu,iv,it,iq,ips,irh

  rmse_u = 0.0d0
  rmse_v = 0.0d0
  rmse_t = 0.0d0
  rmse_q = 0.0d0
  rmse_ps = 0.0d0
  rmse_rh = 0.0d0
  bias_u = 0.0d0
  bias_v = 0.0d0
  bias_t = 0.0d0
  bias_q = 0.0d0
  bias_ps = 0.0d0
  bias_rh = 0.0d0
  iu = 0
  iv = 0
  it = 0
  iq = 0
  ips = 0
  irh = 0
  DO n=1,nn
    IF(qc(n) /= 1) CYCLE
    SELECT CASE(NINT(elm(n)))
    CASE(id_u_obs)
      rmse_u = rmse_u + dep(n)**2
      bias_u = bias_u + dep(n)
      iu = iu + 1
    CASE(id_v_obs)
      rmse_v = rmse_v + dep(n)**2
      bias_v = bias_v + dep(n)
      iv = iv + 1
    CASE(id_t_obs)
      rmse_t = rmse_t + dep(n)**2
      bias_t = bias_t + dep(n)
      it = it + 1
    CASE(id_q_obs)
      rmse_q = rmse_q + dep(n)**2
      bias_q = bias_q + dep(n)
      iq = iq + 1
    CASE(id_ps_obs)
      rmse_ps = rmse_ps + dep(n)**2
      bias_ps = bias_ps + dep(n)
      ips = ips + 1
    CASE(id_rh_obs)
      rmse_rh = rmse_rh + dep(n)**2
      bias_rh = bias_rh + dep(n)
      irh = irh + 1
    END SELECT
  END DO
  IF(iu == 0) THEN
    rmse_u = undef
    bias_u = undef
  ELSE
    rmse_u = SQRT(rmse_u / REAL(iu,r_size))
    bias_u = bias_u / REAL(iu,r_size)
  END IF
  IF(iv == 0) THEN
    rmse_v = undef
    bias_v = undef
  ELSE
    rmse_v = SQRT(rmse_v / REAL(iv,r_size))
    bias_v = bias_v / REAL(iv,r_size)
  END IF
  IF(it == 0) THEN
    rmse_t = undef
    bias_t = undef
  ELSE
    rmse_t = SQRT(rmse_t / REAL(it,r_size))
    bias_t = bias_t / REAL(it,r_size)
  END IF
  IF(iq == 0) THEN
    rmse_q = undef
    bias_q = undef
  ELSE
    rmse_q = SQRT(rmse_q / REAL(iq,r_size))
    bias_q = bias_q / REAL(iq,r_size)
  END IF
  IF(ips == 0) THEN
    rmse_ps = undef
    bias_ps = undef
  ELSE
    rmse_ps = SQRT(rmse_ps / REAL(ips,r_size))
    bias_ps = bias_ps / REAL(ips,r_size)
  END IF
  IF(irh == 0) THEN
    rmse_rh = undef
    bias_rh = undef
  ELSE
    rmse_rh = SQRT(rmse_rh / REAL(irh,r_size))
    bias_rh = bias_rh / REAL(irh,r_size)
  END IF

  WRITE(6,'(A)') '== OBSERVATIONAL DEPARTURE ============================================='
  WRITE(6,'(6A12)') 'U','V','T','Q','PS','RH'
  WRITE(6,'(6ES12.3)') bias_u,bias_v,bias_t,bias_q,bias_ps,bias_rh
  WRITE(6,'(6ES12.3)') rmse_u,rmse_v,rmse_t,rmse_q,rmse_ps,rmse_rh
  WRITE(6,'(A)') '== NUMBER OF OBSERVATIONS TO BE ASSIMILATED ============================'
  WRITE(6,'(6A12)') 'U','V','T','Q','PS','RH'
  WRITE(6,'(6I12)') iu,iv,it,iq,ips,irh
  WRITE(6,'(A)') '========================================================================'

  RETURN
END SUBROUTINE monit_dep
!-----------------------------------------------------------------------
! Basic modules for observation input
!-----------------------------------------------------------------------
SUBROUTINE get_nobs(cfile,nn)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(OUT) :: nn
  REAL(r_sngl) :: wk(7)
  INTEGER :: ios
  INTEGER :: iu,iv,it,iq,irh,ips
  INTEGER :: iunit
  LOGICAL :: ex

  nn = 0
  iu = 0
  iv = 0
  it = 0
  iq = 0
  irh = 0
  ips = 0
  iunit=91
  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
    DO
      READ(iunit,IOSTAT=ios) wk
      IF(ios /= 0) EXIT
      SELECT CASE(NINT(wk(1)))
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
      END SELECT
      nn = nn + 1
    END DO
    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
    WRITE(6,'(A12,I10)') '          U:',iu
    WRITE(6,'(A12,I10)') '          V:',iv
    WRITE(6,'(A12,I10)') '          T:',it
    WRITE(6,'(A12,I10)') '          Q:',iq
    WRITE(6,'(A12,I10)') '         RH:',irh
    WRITE(6,'(A12,I10)') '         Ps:',ips
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF

  RETURN
END SUBROUTINE get_nobs

SUBROUTINE read_obs(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,otyp)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(OUT) :: elem(nn) ! element number
  REAL(r_size),INTENT(OUT) :: rlon(nn)
  REAL(r_size),INTENT(OUT) :: rlat(nn)
  REAL(r_size),INTENT(OUT) :: rlev(nn)
  REAL(r_size),INTENT(OUT) :: odat(nn)
  REAL(r_size),INTENT(OUT) :: oerr(nn)
  REAL(r_size),INTENT(OUT) :: otyp(nn)
  REAL(r_sngl) :: wk(7)
  INTEGER :: n,iunit

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,nn
    READ(iunit) wk
    SELECT CASE(NINT(wk(1)))
    CASE(id_u_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_v_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_t_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_q_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_ps_obs)
      wk(5) = wk(5) * 100.0 ! hPa -> Pa
      wk(6) = wk(6) * 100.0 ! hPa -> Pa
    CASE(id_rh_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = wk(5) * 0.01 ! percent input
      wk(6) = wk(6) * 0.01 ! percent input
    END SELECT
    elem(n) = REAL(wk(1),r_size)
    rlon(n) = REAL(wk(2),r_size)
    rlat(n) = REAL(wk(3),r_size)
    rlev(n) = REAL(wk(4),r_size)
    odat(n) = REAL(wk(5),r_size)
    oerr(n) = REAL(wk(6),r_size)
    otyp(n) = REAL(wk(7),r_size)
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs

SUBROUTINE monit_obs(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,otyp,omb,oma,slot)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elem(nn) ! element number
  REAL(r_size),INTENT(IN) :: rlon(nn)
  REAL(r_size),INTENT(IN) :: rlat(nn)
  REAL(r_size),INTENT(IN) :: rlev(nn)
  REAL(r_size),INTENT(IN) :: odat(nn)
  REAL(r_size),INTENT(IN) :: oerr(nn)
  REAL(r_size),INTENT(IN) :: otyp(nn)
  REAL(r_size),INTENT(IN) :: omb(nn)
  REAL(r_size),INTENT(IN) :: oma(nn)
  REAL(r_size),INTENT(IN) :: slot(nn)
  REAL(r_sngl) :: wk(10)
  INTEGER :: n,iunit

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,nn
    wk(1) = REAL(elem(n),r_sngl)
    wk(2) = REAL(rlon(n),r_sngl)
    wk(3) = REAL(rlat(n),r_sngl)
    wk(4) = REAL(rlev(n),r_sngl)
    wk(5) = REAL(odat(n),r_sngl)
    wk(6) = REAL(oerr(n),r_sngl)
    wk(7) = REAL(otyp(n),r_sngl)
    wk(8) = REAL(omb(n),r_sngl)
    wk(9) = REAL(oma(n),r_sngl)
    wk(10) = REAL(slot(n),r_sngl)
    SELECT CASE(NINT(wk(1)))
    CASE(id_u_obs)
      wk(4) = wk(4) / 100.0 ! Pa -> hPa
    CASE(id_v_obs)
      wk(4) = wk(4) / 100.0 ! Pa -> hPa
    CASE(id_t_obs)
      wk(4) = wk(4) / 100.0 ! Pa -> hPa
    CASE(id_q_obs)
      wk(4) = wk(4) / 100.0 ! Pa -> hPa
    CASE(id_ps_obs)
      wk(5) = wk(5) / 100.0 ! Pa -> hPa
      wk(6) = wk(6) / 100.0 ! Pa -> hPa
    CASE(id_rh_obs)
      wk(4) = wk(4) / 100.0 ! Pa -> hPa
      wk(5) = wk(5) / 0.01 ! percent output
      wk(6) = wk(6) / 0.01 ! percent output
    END SELECT
    WRITE(iunit) wk
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE monit_obs

END MODULE common_obs_wrf
