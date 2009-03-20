
      SUBROUTINE IOGRID (IMODE)
C--
C--   SUBROUTINE IOGRID (IMODE)
C--   Created by TM
C--
C--   Purpose : read or write a gridded file in sigma coordinate
C--   Input :   IMODE = 1 : read model variables from a gridded file (sigma)
C--                   = 2 : write model variables  to a gridded file (p)
C--                   = 3 : write a GrADS control file (for p)
C--                   = 4 : write model variables  to a gridded file (sigma)
C--                   = 5 : write a GrADS control file (for sigma)
C--   Initialized common blocks (if IMODE = 1) : DYNSP1, SFCANOM 
C--

      include "atparam.h"
      include "atparam1.h"

      PARAMETER (NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT)

      include "com_physcon.h"
      include "com_physvar.h"
      include "com_dynvar.h"
      include "com_dyncon1.h"
      include "com_anomvar.h"

      include "com_date.h"
      include "com_tsteps.h"

      COMPLEX UCOSTMP(MX,NX),VCOSTMP(MX,NX)
      REAL UGR(NGP,KX),VGR(NGP,KX),TGR(NGP,KX),QGR(NGP,KX),
     &  PHIGR(NGP,KX),PSGR(NGP)
      REAL UGR1(NGP,KX),VGR1(NGP,KX),TGR1(NGP,KX),QGR1(NGP,KX),
     &  PHIGR1(NGP,KX)
      REAL(4) UGR4(NGP,KX),VGR4(NGP,KX),TGR4(NGP,KX),
     &  QGR4(NGP,KX),PHIGR4(NGP,KX),PSGR4(NGP)

C-- For vertical interpolation !adapted from ppo_tminc.f
      INTEGER K0(ngp)
      REAL W0(ngp), ZOUT(ngp), ZINP(nlev), RDZINP(nlev)

C-- File names etc.
      CHARACTER(LEN=14) :: filename='yyyymmddhh.grd'
      CHARACTER(LEN=14) :: ctlname='yyyymmddhh.ctl'
      CHARACTER(LEN=16) :: filenamep='yyyymmddhh_p.grd'
      CHARACTER(LEN=16) :: ctlnamep='yyyymmddhh_p.ctl'
      CHARACTER(LEN=3) :: CMON3='JAN'
      INTEGER :: irec
      INTEGER :: iitest=0

      IF (IMODE.EQ.1) THEN

C--   1. Read the gridded dataset (Sigma-level)

         READ (2,*,END=200) IYEAR
         READ (2,*,END=200) IMONTH
         READ (2,*,END=200) IDATE
         READ (2,*,END=200) IHOUR

         PRINT '(A,I4.4,A,I2.2,A,I2.2,A,I2.2)',
     &     'Read gridded dataset for year/month/date/hour: ',
     &     IYEAR,'/',IMONTH,'/',IDATE,'/',IHOUR

         OPEN (90,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NGP)
         irec=1
         DO K=KX,1,-1
           READ (90,REC=irec) (UGR4(j,k),j=1,NGP)
           irec=irec+1
         END DO
         DO K=KX,1,-1
           READ (90,REC=irec) (VGR4(j,k),j=1,NGP)
           irec=irec+1
         END DO
         DO K=KX,1,-1
           READ (90,REC=irec) (TGR4(j,k),j=1,NGP)
           irec=irec+1
         END DO
         DO K=KX,1,-1
           READ (90,REC=irec) (QGR4(j,k),j=1,NGP)
           irec=irec+1
         END DO
         READ (90,REC=irec) (PSGR4(j),j=1,NGP)
         CLOSE (90)

         UGR = UGR4
         VGR = VGR4
         TGR = TGR4
         QGR = QGR4 *1.0d3
         PSGR = PSGR4
         PSGR = log(PSGR/P0)
         if(iitest==1) PRINT *,' UGR  :',MINVAL(UGR),MAXVAL(UGR)
         if(iitest==1) PRINT *,' VGR  :',MINVAL(VGR),MAXVAL(VGR)
         if(iitest==1) PRINT *,' TGR  :',MINVAL(TGR),MAXVAL(TGR)
         if(iitest==1) PRINT *,' QGR  :',MINVAL(QGR),MAXVAL(QGR)
         if(iitest==1) PRINT *,' PSGR :',MINVAL(PSGR),MAXVAL(PSGR)

!------------TRASH OF VERTICAL INTERPOLATION FROM P TO SIGMA
!------------IN ORDER TO MAKE P INPUT POSSIBLE BUT FAILED..
!         ZINP(1) = -log(POUT(1))
!         DO k=2,nlev
!           ZINP(k) = -log(POUT(k))
!           RDZINP(k) = 1.0d0/(ZINP(k-1)-ZINP(k))
!         END DO
!
!         DO k=1,KX
!
!           DO j=1,ngp
!             ZOUT(j) = -PSGR(j) - SIGL(k)
!           END DO
!
!           CALL SETVIN(ZINP,RDZINP,ZOUT,ngp,kx,K0,W0)
!
!           CALL VERINT(TGR1(1,k),TGR,ngp,kx,K0,W0)
!
!           CALL VERINT(UGR1(1,k),UGR,ngp,kx,K0,W0)
!           CALL VERINT(VGR1(1,k),VGR,ngp,kx,K0,W0)
!           CALL VERINT(QGR1(1,k),QGR,ngp,kx,K0,W0)
!
!           DO j=1,ngp
!             phi1 = PHIGR(j,K0(j)) +
!     &         0.5*RD*(TGR1(j,k)+TGR(j,K0(j)))*(ZOUT(j)-ZINP(K0(j)))
!             phi2 = PHIGR1(j,K0(j)-1) +
!     &         0.5*RD*(TGR1(j,k)+TGR(j,K0(j)-1))*(ZOUT(j)-ZINP(K0(j)-1))
!             PHIGR1(j,k) = phi1 + W0(j)*(phi2-phi1)
!           END DO
!
!         END DO
!
!         if(iitest==1) PRINT *,' UGR1  :',MINVAL(UGR1),MAXVAL(UGR1)
!         if(iitest==1) PRINT *,' VGR1  :',MINVAL(VGR1),MAXVAL(VGR1)
!         if(iitest==1) PRINT *,' TGR1  :',MINVAL(TGR1),MAXVAL(TGR1)
!         if(iitest==1) PRINT *,' QGR1  :',MINVAL(QGR1),MAXVAL(QGR1)
!         if(iitest==1) PRINT *,' PHIGR1:',MINVAL(PHIGR1),MAXVAL(PHIGR1)
!         if(iitest==1) PRINT *,' PSGR1 :',MINVAL(PSGR1),MAXVAL(PSGR1)
!------------END OF TRASH

C---- Conversion from gridded variable to spectral variable
         DO K=1,KX
           CALL VDSPEC(UGR(1,K),VGR(1,K),VOR(1,1,K,1),DIV(1,1,K,1),2)
           IF(IX.EQ.IY*4) THEN
             CALL TRUNCT(VOR(1,1,K,1))
             CALL TRUNCT(DIV(1,1,K,1))
           END IF
           CALL SPEC(TGR(1,K),T(1,1,K,1))
           CALL SPEC(QGR(1,K),TR(1,1,K,1,1))
         END DO
         CALL SPEC(PSGR(1),PS(1,1,1))

      ELSE IF (IMODE.EQ.2.OR.IMODE.EQ.4) THEN

C--   2. Write date and model variables to the gridded file (2:P,4:sigma)

C---- Conversion from spectral model variable to gridded variable

         DO K=1,KX
           CALL UVSPEC(VOR(1,1,K,1),DIV(1,1,K,1),UCOSTMP,VCOSTMP)
           CALL GRID(UCOSTMP,UGR(1,K),2)
           CALL GRID(VCOSTMP,VGR(1,K),2)
         END DO

         DO K=1,KX
           CALL GRID(T(1,1,K,1),TGR(1,K),1)
           CALL GRID(TR(1,1,K,1,1),QGR(1,K),1)
           CALL GRID(PHI(1,1,K),PHIGR(1,K),1)
         END DO

         CALL GRID(PS(1,1,1),PSGR(1),1)

C---- Vertical interpolation from sigma level to pressure level (ppo_tminc.f)
         IF (IMODE.EQ.2) THEN ! p-level output

         ZINP(1) = -SIGL(1)
         DO k=2,nlev
           ZINP(k) = -SIGL(k)
           RDZINP(k) = 1.0d0/(ZINP(k-1)-ZINP(k))
         END DO

         DO k=1,KX

           DO j=1,ngp
             ZOUT(j) = PSGR(j) - log(POUT(k))
           END DO

           CALL SETVIN(ZINP,RDZINP,ZOUT,ngp,kx,K0,W0)

           CALL VERINT(TGR1(1,k),TGR,ngp,kx,K0,W0)

           DO j=1,ngp
             IF(ZOUT(j).LT.ZINP(nlev)) THEN
               textr = MAX(TGR1(j,k),TGR(j,nlev))
               aref = RD*0.006d0/GG * (ZINP(nlev)-ZOUT(j))
               tref = TGR(j,nlev)*(1.0d0+aref+0.5*aref*aref)
               TGR1(j,k) = textr + 0.7d0*(tref-textr)
             END IF
           END DO

          DO j=1,ngp
            W0(j) = MAX(W0(j),0.0d0)
          END DO

           CALL VERINT(UGR1(1,k),UGR,ngp,kx,K0,W0)
           CALL VERINT(VGR1(1,k),VGR,ngp,kx,K0,W0)
           CALL VERINT(QGR1(1,k),QGR,ngp,kx,K0,W0)

           DO j=1,ngp
             phi1 = PHIGR(j,K0(j)) +
     &         0.5*RD*(TGR1(j,k)+TGR(j,K0(j)))*(ZOUT(j)-ZINP(K0(j)))
             phi2 = PHIGR(j,K0(j)-1) +
     &         0.5*RD*(TGR1(j,k)+TGR(j,K0(j)-1))*(ZOUT(j)-ZINP(K0(j)-1))
             PHIGR1(j,k) = phi1 + W0(j)*(phi2-phi1)
           END DO

         END DO

         ELSE  ! sigma-level output
           UGR1 = UGR
           VGR1 = VGR
           TGR1 = TGR
           QGR1 = QGR
           PHIGR1 = PHIGR
         END IF
C---- Output

         PRINT '(A,I4.4,A,I2.2,A,I2.2,A,I2.2)',
     &     'Write gridded dataset for year/month/date/hour: ',
     &     IYEAR,'/',IMONTH,'/',IDATE,'/',IHOUR

         UGR4=UGR1
         VGR4=VGR1
         TGR4=TGR1
         QGR4=QGR1*1.0d-3 ! kg/kg
         PHIGR4=PHIGR1/GG   ! m
         PSGR4=P0*exp(PSGR)! Pa

         IF (IMODE.EQ.2) THEN
           WRITE (filenamep(1:4),'(I4.4)') IYEAR
           WRITE (filenamep(5:6),'(I2.2)') IMONTH
           WRITE (filenamep(7:8),'(I2.2)') IDATE
           WRITE (filenamep(9:10),'(I2.2)') IHOUR
           OPEN (99,FILE=filenamep,FORM='UNFORMATTED',ACCESS='DIRECT',
     &       RECL=4*IX*IL)
         ELSE
           WRITE (filename(1:4),'(I4.4)') IYEAR
           WRITE (filename(5:6),'(I2.2)') IMONTH
           WRITE (filename(7:8),'(I2.2)') IDATE
           WRITE (filename(9:10),'(I2.2)') IHOUR
           OPEN (99,FILE=filename,FORM='UNFORMATTED',ACCESS='DIRECT',
     &       RECL=4*IX*IL)
         END IF
         irec=1
         DO k=KX,1,-1
           WRITE (99,REC=irec) (UGR4(j,k),j=1,NGP)
           irec=irec+1
         END DO
         DO k=KX,1,-1
           WRITE (99,REC=irec) (VGR4(j,k),j=1,NGP)
           irec=irec+1
         END DO
         DO k=KX,1,-1
           WRITE (99,REC=irec) (TGR4(j,k),j=1,NGP)
           irec=irec+1
         END DO
         DO k=KX,1,-1
           WRITE (99,REC=irec) (QGR4(j,k),j=1,NGP)
           irec=irec+1
         END DO
         IF (IMODE.EQ.2) THEN !Z output is only for p-level
           DO k=KX,1,-1
             WRITE (99,REC=irec) (PHIGR4(j,k),j=1,NGP)
             irec=irec+1
           END DO
         END IF
         WRITE (99,REC=irec) (PSGR4(j),j=1,NGP)
         CLOSE (99)
         if(iitest==1) PRINT *,' UGR  :',MINVAL(UGR4),MAXVAL(UGR4)
         if(iitest==1) PRINT *,' VGR  :',MINVAL(VGR4),MAXVAL(VGR4)
         if(iitest==1) PRINT *,' TGR  :',MINVAL(TGR4),MAXVAL(TGR4)
         if(iitest==1) PRINT *,' QGR  :',MINVAL(QGR4),MAXVAL(QGR4)
         if(iitest==1) PRINT *,' PHIGR:',MINVAL(PHIGR4),MAXVAL(PHIGR4)
         if(iitest==1) PRINT *,' PSGR :',MINVAL(PSGR4),MAXVAL(PSGR4)

      ELSE IF (IMODE.EQ.3.OR.IMODE.EQ.5) THEN

C--   3. Write a GrADS control file (3:p,5:sigma)

         IF (IMONTH.EQ.1) THEN
           CMON3='JAN'
         ELSE IF (IMONTH.EQ.2) THEN
           CMON3='FEB'
         ELSE IF (IMONTH.EQ.3) THEN
           CMON3='MAR'
         ELSE IF (IMONTH.EQ.4) THEN
           CMON3='APR'
         ELSE IF (IMONTH.EQ.5) THEN
           CMON3='MAY'
         ELSE IF (IMONTH.EQ.6) THEN
           CMON3='JUN'
         ELSE IF (IMONTH.EQ.7) THEN
           CMON3='JUL'
         ELSE IF (IMONTH.EQ.8) THEN
           CMON3='AUG'
         ELSE IF (IMONTH.EQ.9) THEN
           CMON3='SEP'
         ELSE IF (IMONTH.EQ.10) THEN
           CMON3='OCT'
         ELSE IF (IMONTH.EQ.11) THEN
           CMON3='NOV'
         ELSE IF (IMONTH.EQ.12) THEN
           CMON3='DEC'
         END IF

         IF (IMODE.EQ.3) THEN !p-level
           WRITE (ctlnamep(1:4),'(I4.4)') IYEAR
           WRITE (ctlnamep(5:6),'(I2.2)') IMONTH
           WRITE (ctlnamep(7:8),'(I2.2)') IDATE
           WRITE (ctlnamep(9:10),'(I2.2)') IHOUR
           OPEN (11,FILE=ctlnamep,FORM='FORMATTED')
           WRITE (11,'(A)') 'DSET ^%y4%m2%d2%h2_p.grd'
         ELSE !sigma-level
           WRITE (ctlname(1:4),'(I4.4)') IYEAR
           WRITE (ctlname(5:6),'(I2.2)') IMONTH
           WRITE (ctlname(7:8),'(I2.2)') IDATE
           WRITE (ctlname(9:10),'(I2.2)') IHOUR
           OPEN (11,FILE=ctlname,FORM='FORMATTED')
           WRITE (11,'(A)') 'DSET ^%y4%m2%d2%h2.grd'
         END IF
         WRITE (11,'(A)') 'TITLE SPEEDY MODEL OUTPUT'
         WRITE (11,'(A)') 'UNDEF -9.99E33'
         WRITE (11,'(A)') 'OPTIONS template big_endian'
         WRITE (11,'(A)') 'XDEF 96 LINEAR 0.0 3.75'
         WRITE (11,'(A,48F8.3)') 'YDEF 48 LEVELS ',
     &     (RADANG(J)*90.0d0/ASIN(1.0d0),J=1,48)
         IF (IMODE.EQ.3) THEN
           WRITE (11,'(A)') 'ZDEF 7 LEVELS 925 850 700 500 300 200 100'
         ELSE
           WRITE (11,'(A,7F6.3)') 'ZDEF 7 LEVELS ',(SIG(K),K=7,1,-1)
         END IF
         IF (NDAYSL.NE.0) THEN
           WRITE (11,'(A,I4,A,I2.2,A,I2.2,A,I4.4,A)') 'TDEF ',
     &       NDAYSL*4+1,' LINEAR ',IHOUR,'Z',IDATE,CMON3,IYEAR,' 6HR'
         ELSE
           WRITE (11,'(A,I4,A,I2.2,A,I2.2,A,I4.4,A)') 'TDEF ',
     &       2,' LINEAR ',IHOUR,'Z',IDATE,CMON3,IYEAR,' 6HR'
         END IF
         IF (IMODE.EQ.3) THEN !p-level
           WRITE (11,'(A)') 'VARS 6'
         ELSE !sigma-level
           WRITE (11,'(A)') 'VARS 5'
         END IF
         WRITE (11,'(A)') 'U 7 99 U-wind [m/s]'
         WRITE (11,'(A)') 'V 7 99 V-wind [m/s]'
         WRITE (11,'(A)') 'T 7 99 Temperature [K]'
         WRITE (11,'(A)') 'Q 7 99 Specific Humidity [kg/kg]'
         IF (IMODE.EQ.3) THEN
           WRITE (11,'(A)') 'Z 7 99 Geopotential Height [m]'
         END IF
         WRITE (11,'(A)') 'PS 0 99 Surface Pressure [Pa]'
         WRITE (11,'(A)') 'ENDVARS'
         CLOSE (11)

      ELSE
        PRINT *,'Hey, look at the usage! (IOGRID)'
        STOP
      ENDIF
C--
      RETURN

C--   4. Stop integration if gridded file is not found

  200 CONTINUE

      print*, ' Hey, what are you doing?',
     &  ' fort.2 should contain time setting'

      STOP 'invalid gridded data input'

C--
      END

