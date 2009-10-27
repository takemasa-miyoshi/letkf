PROGRAM dec_prepbufr
!
! NOTE: output goes to fort.90
!
  USE common
  USE common_obs

  IMPLICIT NONE

  INTEGER,PARAMETER :: maxlev = 255     !Maximum number of BUFR levels
  INTEGER,PARAMETER :: maxevn = 10      !Maximum number of BUFR event sequences
!  CHARACTER(MXSTRL) :: head = 'SID XOB YOB DHR ELV TYP T29 ITP'
!  CHARACTER(MXSTRL) :: ostr(MXR8VT) = (/'POB PQM PPC PRC PFC PAN CAT',&
!                                    &   'QOB QQM QPC QRC QFC QAN CAT',&
!                                    &   'TOB TQM TPC TRC TFC TAN CAT',&
!                                    &   'ZOB ZQM ZPC ZRC ZFC ZAN CAT',&
!                                    &   'UOB WQM WPC WRC UFC UAN CAT',&
!                                    &   'VOB WQM WPC WRC VFC VAN CAT'/)
  CHARACTER(8) :: obtype
!  CHARACTER :: var(8) = (/'P','Q','T','Z','U','V'/)
  INTEGER,PARAMETER :: nobtype = 20
  INTEGER :: IREADNS
  INTEGER :: idate,idummy,nlev,n
  INTEGER :: ilev,ievn
  INTEGER :: iobs(nobtype+1),iobs_out(nobtype+1)
  CHARACTER(6) :: obtypelist(nobtype)
  REAL(r_dble) :: station(5)
  CHARACTER(8) :: cs
  REAL(r_dble) :: prs(3,maxlev,maxevn)
  REAL(r_dble) :: obs(3,maxlev,maxevn)
  REAL(r_sngl) :: wk(6)
  INTEGER :: iunit
  !
  ! Open the input file
  !
  OPEN(11,FILE='prepbufr.in',FORM='unformatted')
  CALL OPENBF(11,'IN',11)
  CALL DATELEN(10)
  !
  ! Main loop
  !
  obtypelist = (/'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR',&
               & 'VADWND', 'SATEMP', 'ADPSFC', 'SFCSHP', 'SFCBOG',&
               & 'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND', 'QKSWND',&
               & 'MSONET', 'GPSIPW', 'RASSDA', 'WDSATR', 'ASCATW'/)
  iobs = 0
  iobs_out = 0
  DO
    !
    ! next record
    !
    IF(IREADNS(11,obtype,idate) /= 0) EXIT
!    print *,obtype,idate ! for debugging
!    IF(n < 5) CYCLE      ! for debugging
!    READ(*,*)            ! for debugging
    DO n=1,nobtype+1
      IF(obtype == obtypelist(n) .OR. n > nobtype) THEN
        iobs(n) = iobs(n)+1
        EXIT
      END IF
    END DO
    !
    ! station location (lon, lat)
    !
    CALL UFBINT(11,station,5,1,idummy,'SID XOB YOB ELV DHR')
    WRITE(cs(1:8),'(A8)') station(1)
!    IF(n < 2) CYCLE                                       ! for debugging
!    print '(A,4F10.3,I)',cs,station(2:5),NINT(station(5)) ! for debugging
!    READ(*,*)                                             ! for debugging
    wk(2:3) = station(2:3)
    wk(4) = station(4)
    IF(NINT(station(5)) < -3 .OR. 3 < NINT(station(5))) CYCLE
    iunit = 90+NINT(station(5))
    !
    ! obs
    !
    CALL UFBEVN(11,prs,3,maxlev,maxevn,nlev,'POB POE PQM')
    IF(obtype == 'ADPSFC' .OR.&
     & obtype == 'SFCSHP' .OR.&
     & obtype == 'SFCBOG') CALL output_ps ! surface pressure report
    CALL UFBEVN(11,obs,3,maxlev,maxevn,ilev,'QOB QOE QQM')
    CALL output(id_q_obs)
    CALL UFBEVN(11,obs,3,maxlev,maxevn,ilev,'TOB TOE TQM')
    CALL output(id_t_obs)
    CALL UFBEVN(11,obs,3,maxlev,maxevn,ilev,'UOB WOE WQM')
    CALL output(id_u_obs)
    CALL UFBEVN(11,obs,3,maxlev,maxevn,ilev,'VOB WOE WQM')
    CALL output(id_v_obs)
  END DO

  PRINT '(A)','================================================================================'
  PRINT '(A)','                  SYNOPSIS OF PREPBUFR DECODER using BUFRLIB'
  PRINT '(A)','                              by TAKEMASA MIYOSHI'
  PRINT '(A)','--------------------------------------------------------------------------------'
  PRINT '(A,I10)',' TOTAL NUMBER OF READ-IN RECORDS:',SUM(iobs)
  PRINT '(A,I10,A)',' TOTAL NUMBER OF WRITTEN RECORDS:',SUM(iobs_out),' (levels/variables separately)'
  PRINT '(A)','--------------------------------------------------------------------------------'
  PRINT '(10(2X,A))',obtypelist(1:10)
  PRINT '(10I8)',iobs(1:10)
  PRINT '(10I8)',iobs_out(1:10)
  PRINT '(A)','--------------------------------------------------------------------------------'
  PRINT '(10(2X,A))',obtypelist(11:20)
  PRINT '(10I8)',iobs(11:20)
  PRINT '(10I8)',iobs_out(11:20)
  PRINT '(A)','--------------------------------------------------------------------------------'
  PRINT '(2X,A)','OTHERS'
  PRINT '(I8)',iobs(21)
  PRINT '(I8)',iobs_out(21)
  PRINT '(A)','================================================================================'

  STOP
CONTAINS
SUBROUTINE output(id)
  INTEGER,INTENT(IN) :: id
  INTEGER :: iqm

  IF(ilev /= nlev) THEN
    PRINT *,'FATAL ERROR, nlev /= ilev',nlev,ilev
    STOP
  END IF
  wk(1) = id
  DO ilev=1,nlev
    wk(4) = prs(1,ilev,1) * 100.0d0 ! hPa -> Pa
    iqm = NINT(prs(3,ilev,1))
    IF(iqm < 0 .OR. 2 < iqm) CYCLE
    wk(5:6) = obs(1:2,ilev,1)
    IF(id == id_q_obs) wk(5:6) = wk(5:6) * 1.E-6 ! Mg/kg -> g/kg
    iqm = NINT(obs(3,ilev,1))
    IF(iqm < 0 .OR. 2 < iqm) CYCLE
    IF(wk(6) > 1.E10) CYCLE
    WRITE(iunit) wk
    iobs_out(n) = iobs_out(n) + 1
  END DO

  RETURN
END SUBROUTINE output

SUBROUTINE output_ps
  INTEGER :: iqm

  wk(1) = id_ps_obs
  iqm = NINT(prs(3,1,1))
  IF(iqm < 0 .OR. 2 < iqm) RETURN
  wk(5:6) = prs(1:2,1,1)
  IF(wk(6) > 1.E10) RETURN
  WRITE(iunit) wk
  iobs_out(n) = iobs_out(n) + 1

  RETURN
END SUBROUTINE output_ps

END PROGRAM dec_prepbufr
