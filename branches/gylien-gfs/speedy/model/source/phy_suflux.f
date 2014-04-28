
      SUBROUTINE SUFLUX (PSA,UA,VA,TA,QA,RH,PHI,
     &                   PHI0,FMASK,TLAND,TSEA,SWAV,SSR,SLRD,
     &                   USTR,VSTR,SHF,EVAP,SLRU,HFLUXN,
     &                   TSFC,TSKIN,U0,V0,T0,Q0)
C--
C--   SUBROUTINE SUFLUX (PSA,UA,VA,TA,QA,RH,PHI,
C--  &                   PHI0,FMASK,TLAND,TSEA,SWAV,SSR,SLRD,
C--  &                   USTR,VSTR,SHF,EVAP,SLRU,HFLUXN,
C--  &                   TSFC,TSKIN,U0,V0,T0,Q0)
C--
C--   Purpose: Compute surface fluxes of momentum, energy and moisture,
C--            and define surface skin temperature from energy balance
C--   Input:   PSA    = norm. surface pressure [p/p0]   (2-dim)
C--            UA     = u-wind                          (3-dim)
C--            VA     = v-wind                          (3-dim)
C--            TA     = temperature                     (3-dim)
C--            QA     = specific humidity [g/kg]        (3-dim)
C--            RH     = relative humidity [0-1]         (3-dim)
C--            PHI    = geopotential                    (3-dim)
C--            PHI0   = surface geopotential            (2-dim)
C--            FMASK  = fractional land-sea mask        (2-dim)
C--            TLAND  = land-surface temperature        (2-dim)
C--            TSEA   =  sea-surface temperature        (2-dim)
C--            SWAV   = soil wetness availability [0-1] (2-dim)
C--            SSR    = sfc sw radiation (net flux)     (2-dim)
C--            SLRD   = sfc lw radiation (downward flux)(2-dim)
C--   Output:  USTR   = u stress                        (2-dim)
C--            VSTR   = v stress                        (2-dim)
C--            SHF    = sensible heat flux              (2-dim)
C--            EVAP   = evaporation [g/(m^2 s)]         (2-dim)
C--            SLRU   = sfc lw radiation (upward flux)  (2-dim)
C--            HFLUXN = net heat flux into land/sea     (2-dim)           
C--            TSFC   = surface temperature (clim.)     (2-dim)
C--            TSKIN  = skin surface temperature        (2-dim)
C--            U0     = near-surface u-wind             (2-dim)
C--            V0     = near-surface v-wind             (2-dim)
C--            T0     = near-surface air temperature    (2-dim)
C--            Q0     = near-surface sp. humidity [g/kg](2-dim)
C--
C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Physical constants + functions of sigma and latitude

      include "com_physcon.h"

C     Surface flux constants

      include "com_sflcon.h"      

      include "com_radcon.h"

      REAL PSA(NGP), UA(NGP,NLEV), VA(NGP,NLEV), TA(NGP,NLEV),
     &     QA(NGP,NLEV), RH(NGP,NLEV), PHI(NGP,NLEV),
     &     PHI0(NGP), FMASK(NGP), TLAND(NGP), TSEA(NGP), SWAV(NGP),
     &     SSR(NGP), SLRD(NGP)

      REAL USTR(NGP,3), VSTR(NGP,3), SHF(NGP,3), EVAP(NGP,3),
     &     SLRU(NGP), HFLUXN(NGP,2), TSFC(NGP), TSKIN(NGP),
     &     U0(NGP), V0(NGP), T0(NGP), Q0(NGP)
									
      REAL T1(NGP), QSAT0(NGP,2), DENVV(NGP), CDENVV(NGP,2)


C--   1. Extrapolation of wind, temp, hum. and density to the surface

C     1.1 Wind components
   
      DO J=1,NGP
        U0(J) = FWIND0*UA(J,NLEV)
        V0(J) = FWIND0*VA(J,NLEV)
      ENDDO

C     1.2 Temperature

      GTEMP0 = 1.-FTEMP0
      RCP = 1./CP
      NL1=NLEV-1
C
      DO J=1,NGP
        T0(J) = TA(J,NLEV)+WVI(NLEV,2)*(TA(J,NLEV)-TA(J,NL1))
        T1(J) = TA(J,NLEV)+RCP*(PHI(J,NLEV)-PHI0(J))
      ENDDO

      DO J=1,NGP
        T0(J) = FTEMP0*T0(J)+GTEMP0*T1(J)
        IF (TA(J,NLEV).LE.TA(J,NL1)) T0(J) = TA(J,NLEV)
      ENDDO

C     1.3 Spec. humidity

      GHUM0 = 1.-FHUM0

      CALL SHTORH (-1,NGP,T0,PSA,1.,Q0,RH(1,NLEV),QSAT0)

      DO J=1,NGP
        Q0(J)=FHUM0*Q0(J)+GHUM0*QA(J,NLEV)
      ENDDO

C     1.4 Density * wind speed (including gustiness factor)

      PRD = P0/RD
      VG2 = VGUST*VGUST

      DO J=1,NGP
        DENVV(J)=(PRD*PSA(J)/T0(J))*
     &           SQRT(U0(J)*U0(J)+V0(J)*V0(J)+VG2)
      ENDDO

C     1.5 Define effective skin temperature to compensate for
C         non-linearity of heat/moisture fluxes during the daily cycle

      DO JLAT=1,NLAT
	J0=NLON*(JLAT-1)
        SQCLAT=SQRT(CLAT(JLAT))
        DO J=J0+1,J0+NLON
          TSKIN(J)=TLAND(J)+CTDAY*SQCLAT*SSR(J)*PSA(J)
        ENDDO
      ENDDO


C--   2. Computation of fluxes over land and sea

C     2.1 Wind stress

C     Orographic correction

      DO J=1,NGP
        CDENVV(J,1)=CDL*DENVV(J)*FOROG(J)
        CDENVV(J,2)=CDS*DENVV(J)
      ENDDO

      DO J=1,NGP
        USTR(J,1) = -CDENVV(J,1)*UA(J,NLEV)
        VSTR(J,1) = -CDENVV(J,1)*VA(J,NLEV)
        USTR(J,2) = -CDENVV(J,2)*UA(J,NLEV)
        VSTR(J,2) = -CDENVV(J,2)*VA(J,NLEV)
      ENDDO

C     2.2 Sensible heat flux 

C     Stability correction

      RDTH = FSTAB/DTHETA

      DO J=1,NGP
        FSLAND=1.+MIN(DTHETA,MAX(-DTHETA,TSKIN(J)-T1(J)))*RDTH
        FSSEA =1.+MIN(DTHETA,MAX(-DTHETA, TSEA(J)-T1(J)))*RDTH
        CDENVV(J,1)=CHL*DENVV(J)*FSLAND
        CDENVV(J,2)=CHS*DENVV(J)*FSSEA
      ENDDO

      DO J=1,NGP
        SHF(J,1) = CDENVV(J,1)*CP*(TSKIN(J)-T0(J))
        SHF(J,2) = CDENVV(J,2)*CP*(TSEA(J) -T0(J))
      ENDDO

C     2.3 Evaporation

      CALL SHTORH (0,NGP,TSKIN,PSA,1.,QDUMMY,RDUMMY,QSAT0(1,1))
      CALL SHTORH (0,NGP,TSEA ,PSA,1.,QDUMMY,RDUMMY,QSAT0(1,2))

      DO J=1,NGP
C       EVAP(J,1) = CDENVV(J,1)*SWAV(J)*MAX(0.,QSAT0(J,1)-Q0(J))
        EVAP(J,1) = CDENVV(J,1)*MAX(0.,SWAV(J)*QSAT0(J,1)-Q0(J))
        EVAP(J,2) = CDENVV(J,2)*              (QSAT0(J,2)-Q0(J))
      ENDDO


C--   3. Weighted average of surface fluxes and temperatures 
C--      according to land-sea mask

      DO J=1,NGP
        USTR(J,3) = USTR(J,2)+FMASK(J)*(USTR(J,1)-USTR(J,2))
        VSTR(J,3) = VSTR(J,2)+FMASK(J)*(VSTR(J,1)-VSTR(J,2))
         SHF(J,3) =  SHF(J,2)+FMASK(J)*( SHF(J,1)- SHF(J,2))
        EVAP(J,3) = EVAP(J,2)+FMASK(J)*(EVAP(J,1)-EVAP(J,2))
      ENDDO

      DO J=1,NGP
        TSFC(J)  = TSEA(J)+FMASK(J)*(TLAND(J)-TSEA(J))
        TSKIN(J) = TSEA(J)+FMASK(J)*(TSKIN(J)-TSEA(J))
      ENDDO


C--   4. Emission of lw radiation from the surface
C--      and net heat flux into land and sea surface

      ESBC=EMISFC*SBC

      DO J=1,NGP

        SLRL    = ESBC*TSKIN(J)**4
        SLRS    = ESBC*TSEA(J) **4
        SLRU(J) = SLRS+FMASK(J)*(SLRL-SLRS)
        
        HFLUXN(J,1) = SSR(J)+SLRD(J)-(SLRL+SHF(J,1)+ALHC*EVAP(J,1))
        HFLUXN(J,2) = SSR(J)+SLRD(J)-(SLRS+SHF(J,2)+ALHC*EVAP(J,2))

      ENDDO

      RETURN
      END

      SUBROUTINE SFLSET (PHI0)
C--
C--   SUBROUTINE SFLSET (PHI0)
C--
C--   Purpose: compute orographic factor for land surface drag
C--   Input:   PHI0   = surface geopotential            (2-dim)
C--            Initialized common blocks: SFLFIX

C     Resolution parameters
C
      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Physical constants + functions of sigma and latitude
      include "com_physcon.h"

C     Surface flux constants
      include "com_sflcon.h"

      REAL PHI0(NGP)

      RHDRAG = 1./(GG*HDRAG)

      DO J=1,NGP
        FOROG(J)=1.+FHDRAG*(1.-EXP(-MAX(PHI0(J),0.)*RHDRAG))
      ENDDO

C--
      RETURN
      END
