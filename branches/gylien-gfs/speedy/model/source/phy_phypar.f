      SUBROUTINE PHYPAR (VOR1,DIV1,T1,Q1,PHI1,PSL1,
     &                   UTEND,VTEND,TTEND,QTEND)
C--
C--   SUBROUTINE PHYPAR (VOR1,DIV1,T1,Q1,PHI1,PSL1,
C--  &                   UTEND,VTEND,TTEND,QTEND)
C--
C--   Purpose: compute physical parametrization tendencies for U, V, T, Q 
C--   and add them to dynamical grid-point tendencies
C--   Input-only  arguments:   VOR1   : vorticity (sp)
C--                            DIV1   : divergence (sp)
C--                            T1     : temperature (sp)
C--                            Q1     : specific humidity (sp)
C--                            PHI1   : geopotential (sp)
C--                            PSL1   : log of sfc pressure (sp)
C--   Input-output arguments:  UTEND  : u-wind tendency (gp)
C--                            VTEND  : v-wind tendency (gp)
C--                            TTEND  : temp. tendency (gp)
C--                            QTEND  : spec. hum. tendency (gp)
C--   Modified common blocks:  PHYGR1, PHYGR2, PHYGR3, PHYTEN, FLUXES
C--

C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Constants + functions of sigma and latitude
      include "com_physcon.h"

C     Model variables, tendencies and fluxes on gaussian grid
      include "com_physvar.h"

C     Surface forcing fields (time-inv. or functions of seasonal cycle)
      include "com_forcing.h"

C--   Random diabatic forcing 
      include "com_randfor.h"

C     Logical flags
      include "com_lflags.h"

      COMPLEX VOR1(MX,NX,NLEV), DIV1(MX,NX,NLEV), T1(MX,NX,NLEV),
     &          Q1(MX,NX,NLEV), PHI1(MX,NX,NLEV), PSL1(MX,NX),
     &          UCOS(MX,NX), VCOS(MX,NX)

      REAL UTEND(NGP,NLEV), VTEND(NGP,NLEV), TTEND(NGP,NLEV),
     &     QTEND(NGP,NLEV)

      INTEGER IPTOP(NGP), ICLTOP(NGP), ICNV(NGP)
      REAL    RPS(NGP), ST4S(NGP)

      iitest=0

C--   1. Compute grid-point fields

C     1.1 Convert model spectral variables to grid-point variables

      if (iitest.eq.1) print *, ' 1.1 in PHYPAR'

      DO K=1,NLEV

        CALL UVSPEC (VOR1(1,1,K),DIV1(1,1,K),UCOS,VCOS)
        CALL GRID   (UCOS,UG1(1,K),2)
        CALL GRID   (VCOS,VG1(1,K),2)

      ENDDO

      DO K=1,NLEV

        CALL GRID   (T1(1,1,K),  TG1(1,K),  1)
        CALL GRID   (Q1(1,1,K),  QG1(1,K),  1)
        CALL GRID   (PHI1(1,1,K),PHIG1(1,K),1)

      ENDDO

      CALL GRID (PSL1,PSLG1,1)

C     Remove negative humidity values
C     CALL QNEG (QG1)

C     1.2 Compute thermodynamic variables

      if (iitest.eq.1) print *, ' 1.2 in PHYPAR'

      DO J=1,NGP
       PSG(J)=EXP(PSLG1(J))
       RPS(J)=1./PSG(J)
      ENDDO

      DO K=1,NLEV
       DO J=1,NGP
c       remove when QNEG is implemented
	qg1(j,k)=max(qg1(j,k),0.)
        SE(J,K)=CP*TG1(J,K)+PHIG1(J,K)
       ENDDO
      ENDDO

      DO K=1,NLEV
       CALL SHTORH (1,NGP,TG1(1,K),PSG,SIG(K),QG1(1,K),
     &              RH(1,K),QSAT(1,K))
      ENDDO

C--   2. Precipitation 

C     2.1 Deep convection

      CALL CONVMF (PSG,SE,QG1,QSAT,
     &             IPTOP,CBMF,PRECNV,TT_CNV,QT_CNV)

      DO K=2,NLEV
       DO J=1,NGP
        TT_CNV(J,K)=TT_CNV(J,K)*RPS(J)*GRDSCP(K)
        QT_CNV(J,K)=QT_CNV(J,K)*RPS(J)*GRDSIG(K)
       ENDDO
      ENDDO

      DO J=1,NGP
        ICNV(J)=NLEV-IPTOP(J)
      ENDDO

C     2.2 Large-scale condensation

      CALL LSCOND (PSG,QG1,QSAT,
     &             IPTOP,PRECLS,TT_LSC,QT_LSC)

      DO K=2,NLEV
       DO J=1,NGP
        TTEND(J,K)=TTEND(J,K)+TT_CNV(J,K)+TT_LSC(J,K)
        QTEND(J,K)=QTEND(J,K)+QT_CNV(J,K)+QT_LSC(J,K)
       ENDDO
      ENDDO


C--   3. Radiation (shortwave and longwave) and surface fluxes

C     3.1 Compute shortwave tendencies and initialize lw transmissivity

      if (iitest.eq.1) print *, ' 3.1 in PHYPAR'

C     The sw radiation may be called at selected time steps

      IF (LRADSW) THEN

        CALL CLOUD (QG1,RH,PRECNV,PRECLS,IPTOP,
     &              ICLTOP,CLOUDC)

        DO J=1,NGP
          CLTOP(J)=SIGH(ICLTOP(J)-1)*PSG(J)
        ENDDO

        CALL RADSW (PSG,QG1,ALB1,ICLTOP,CLOUDC,
     &              TSR,SSR,TT_RSW)

        DO K=1,NLEV
         DO J=1,NGP
          TT_RSW(J,K)=TT_RSW(J,K)*RPS(J)*GRDSCP(K)
         ENDDO
        ENDDO

      ENDIF

C     3.2 Compute downward longwave fluxes 

      CALL RADLW (-1,TG1,TS,ST4S,
     &            OLR,SLR,TT_RLW)

C     3.3. Compute surface fluxes and land skin temperature
 
      CALL SUFLUX (PSG,UG1,VG1,TG1,QG1,RH,PHIG1,
     &             PHIS0,FMASK1,STL1,SST1,SOILW1,SSR,SLR,
     &             USTR,VSTR,SHF,EVAP,ST4S,HFLUXN,
     &             TS,TSKIN,U0,V0,T0,Q0) 
 
       CALL ADDFLX (HFLUXN)     

C     3.4 Compute upward longwave fluxes, convert them to tendencies 
C         and add shortwave tendencies

      if (iitest.eq.1) print *, ' 3.4 in PHYPAR'

      CALL RADLW (1,TG1,TS,ST4S,
     &            OLR,SLR,TT_RLW)

      DO K=1,NLEV
       DO J=1,NGP
        TT_RLW(J,K)=TT_RLW(J,K)*RPS(J)*GRDSCP(K)
        TTEND (J,K)=TTEND(J,K)+TT_RSW(J,K)+TT_RLW(J,K)
       ENDDO
      ENDDO

C--   4. PBL interactions with lower troposphere

C     4.1 Vertical diffusion and shallow convection

      CALL VDIFSC (UG1,VG1,SE,RH,QG1,QSAT,PHIG1,ICNV,
     &             UT_PBL,VT_PBL,TT_PBL,QT_PBL)

C     4.2 Add tendencies due to surface fluxes 

      DO J=1,NGP
       UT_PBL(J,NLEV)=UT_PBL(J,NLEV)+USTR(J,3)*RPS(J)*GRDSIG(NLEV)
       VT_PBL(J,NLEV)=VT_PBL(J,NLEV)+VSTR(J,3)*RPS(J)*GRDSIG(NLEV)
       TT_PBL(J,NLEV)=TT_PBL(J,NLEV)+ SHF(J,3)*RPS(J)*GRDSCP(NLEV)
       QT_PBL(J,NLEV)=QT_PBL(J,NLEV)+EVAP(J,3)*RPS(J)*GRDSIG(NLEV)
      ENDDO

      DO K=1,NLEV
       DO J=1,NGP
        UTEND(J,K)=UTEND(J,K)+UT_PBL(J,K)
        VTEND(J,K)=VTEND(J,K)+VT_PBL(J,K)
        TTEND(J,K)=TTEND(J,K)+TT_PBL(J,K)
        QTEND(J,K)=QTEND(J,K)+QT_PBL(J,K)
       ENDDO
      ENDDO

C--   5. Random diabatic forcing 

      IF (LRANDF) THEN

        DO K=1,NLEV
         DO J=1,NGP
          TT_PBL(J,K)=RANDFH(J,1)*RANDFV(K,1)
     &               +RANDFH(J,2)*RANDFV(K,2)    
          TTEND(J,K)=TTEND(J,K)+RANDFH(J,1)*RANDFV(K,1)
     &                         +RANDFH(J,2)*RANDFV(K,2)
         ENDDO
        ENDDO

      ENDIF      

C--
      RETURN
      END


