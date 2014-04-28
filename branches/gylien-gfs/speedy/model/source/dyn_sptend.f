      SUBROUTINE SPTEND (DIVDT,TDT,PSDT,J4)
C--
C--   SUBROUTINE SPTEND (DIVDT,TDT,PSDT,J4)
C--
C--   Purpose : compute spectral tendencies of divergence, temperature
C--             and log_surf.pressure)
C--   Input/output : DIVDT = divergence tendency (spec.)
C--                  TDT   = temperature tendency (spec.)
C--                  PSDT  = tendency of log_surf.pressure (spec.)
C--                  J4    = time level index (1 or 2)
C--
      include "atparam.h"
      include "atparam1.h"

      include "com_dyncon1.h"
      include "com_dyncon2.h"
      include "com_dynvar.h"

      COMPLEX DIVDT(MX,NX,KX), TDT(MX,NX,KX), PSDT(MX,NX)

      COMPLEX DUMK(MX,NX,KXP),DMEANC(MX,NX),SIGDTC(MX,NX,KXP)
      COMPLEX TEMPC(MX,NX,3)
      COMPLEX DUMC(MX,NX,2),ZERO

      ZERO=(0.,0.)

c   vertical mean div and pressure tendency

      DO 512 M=1,MXNX
        DMEANC(M,1)=ZERO
  512 CONTINUE

      DO 53 K=1,KX
        DO 513 M=1,MXNX
          DMEANC(M,1)=DMEANC(M,1)+DIV(M,1,K,J4)*DHS(K)
  513   CONTINUE
   53 CONTINUE

      DO 114 M=1,MXNX
        PSDT(M,1)=PSDT(M,1)-DMEANC(M,1)
  114 CONTINUE
      PSDT(1,1)=ZERO

c  sigma-dot "velocity" and temperature tendency

      DO 116 M=1,MXNX
        SIGDTC(M,1,1)=ZERO
        SIGDTC(M,1,KXP)=ZERO
  116 CONTINUE

      DO 217 K=1,KXM
        DO 218 M=1,MXNX
          SIGDTC(M,1,K+1)=SIGDTC(M,1,K)
     *     -DHS(K)*(DIV(M,1,K,J4)-DMEANC(M,1))
  218   CONTINUE
  217 CONTINUE

      DO 5701 M=1,MXNX
        DUMK(M,1,1)=ZERO
        DUMK(M,1,KXP)=ZERO
 5701 CONTINUE

      DO 5702 K=2,KX
      DO 5702 M=1,MXNX
        DUMK(M,1,K)=SIGDTC(M,1,K)*(TREF(K)-TREF(K-1))
 5702 CONTINUE

      DO 5502 K=1,KX
      DO 5502 M=1,MXNX
        TDT(M,1,K)= TDT(M,1,K)-(DUMK(M,1,K+1)+DUMK(M,1,K))*DHSR(K)
     *     +TREF3(K)*(SIGDTC(M,1,K+1)+SIGDTC(M,1,K))
     *     -TREF2(K)*DMEANC(M,1)
 5502 CONTINUE

c   geopotential and divergence tendency

      CALL GEOP(J4)

      DO 18 K=1,KX
        DO 221 M=1,MXNX
          DUMC(M,1,1)=PHI(M,1,K)+RGAS*TREF(K)*PS(M,1,J4)
  221 CONTINUE
        CALL LAP(DUMC(1,1,1),DUMC(1,1,2))
        DO 222 M=1,MXNX
          DIVDT(M,1,K)=DIVDT(M,1,K)-DUMC(M,1,2)
  222   CONTINUE
   18 CONTINUE

      RETURN
      END   
