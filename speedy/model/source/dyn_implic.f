      SUBROUTINE IMPLIC (DIVDT,TDT,PSDT)
C--
C--   SUBROUTINE IMPLIC (DIVDT,TDT,PSDT)
C--
C--   Purpose : Correct tendencies for implicit gravity wave model
C--   Input/output : DIVDT = divergence tendency
C--                  TDT   = temperature tendency
C--                  PSDT  = tendency of log(surf.pressure)
C--
      include "atparam.h"
      include "atparam1.h"
      PARAMETER (MXNXKX=MX*NX*KX)

      include "com_dyncon1.h"
      include "com_dyncon2.h"

      COMPLEX DIVDT(MX,NX,KX),TDT(MX,NX,KX),PSDT(MX,NX)
      COMPLEX YE(MX,NX,KX),YF(MX,NX,KX),ZERO

      ZERO=(0.,0.)

      DO 1 K=1,KX
      DO 1 N=1,NX
      DO 1 M=1,MX
        YE(M,N,K)=ZERO
    1 CONTINUE

      DO 2 K1=1,KX
      DO 2 K=1,KX
        DO 11 M=1,MXNX
          YE(M,1,K)=YE(M,1,K)+XD(K,K1)*TDT(M,1,K1)
   11   CONTINUE
    2 CONTINUE

      DO 21 K=1,KX
      DO 22 M=1,MXNX
        YE(M,1,K)=YE(M,1,K)+TREF1(K)*PSDT(M,1)
   22 CONTINUE
   21 CONTINUE

      DO 4 K=1,KX
        DO 44 M=1,MXNX
          YF(M,1,K)=DIVDT(M,1,K)+ELZ(M,1)*YE(M,1,K)
   44   CONTINUE
    4 CONTINUE

      DO 5 M=1,MXNXKX
        DIVDT(M,1,1)=ZERO
    5 CONTINUE

      DO 6 N=1,NX
      DO 6 M=1,MX
        MM=ISC*(M-1)+1
        LL=MM+N-2
        IF(LL.NE.0) THEN
          DO 66 K1=1,KX
          DO 66 K=1,KX
            DIVDT(M,N,K)=DIVDT(M,N,K)+XJ(K,K1,LL)*YF(M,N,K1)
   66     CONTINUE
        ENDIF
    6 CONTINUE

      DO 46 K=1,KX
        DO 47 M=1,MXNX
          PSDT(M,1)=PSDT(M,1)-DIVDT(M,1,K)*DHSX(K)
   47   CONTINUE
   46 CONTINUE

      DO 7 K=1,KX
      DO 7 K1=1,KX
        DO 48 M=1,MXNX
          TDT(M,1,K)=TDT(M,1,K)+XC(K,K1)*DIVDT(M,1,K1)
   48   CONTINUE
    7 CONTINUE

      RETURN
      END
