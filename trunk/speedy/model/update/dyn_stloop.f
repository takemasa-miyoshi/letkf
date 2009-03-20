      SUBROUTINE STLOOP (ISTEP)
C--
C--   SUBROUTINE STLOOP (ISTEP)
C--
C--   Purpose: Perform a series of time steps calling 
C--            post-processing/output routines at selected steps
C--   Input/output : ISTEP = time step index
C--   Updated common block : LFLAG2
C-- 

      include "com_tsteps.h"
      include "com_lflags.h"
      include "com_date.h"                                              !TM
 
      iitest=0

!     DO J=1,NSTEPS                                                     !TM
      DO JJ=1,4                                                         !TM
      DO J=1,NSTEPS/4                                                   !TM

        if (iitest.eq.1) print*, 'STLOOP: calling step ', istep

C       Set logical flags

        LRADSW = (MOD(ISTEP,NSTRAD).EQ.1)
        LRANDF = ((ISTEP.LE.NSTRDF).OR.(NSTRDF.LT.0))

C       Perform one leapfrog time step

        CALL STEP (2,2,DELT2,ALPH,ROB)   

C       Do diagnostic, post-processing and I/O tasks 
 
        CALL DIAGNS (2,ISTEP)
        IF (IHOUT.NE.1) THEN !time-mean are omitted                     !TM
        CALL TMINC_DAILY
        IF (MOD(ISTEP,NSTPPR).EQ.0) CALL TMINC

        IF (MOD(ISTEP,NSTOUT).EQ.0) CALL TMOUT (1)
        IF (MOD(ISTEP,NSTEPS).EQ.0) CALL TMOUT_DAILY (1)
        END IF                                                          !TM

        ISTEP=ISTEP+1

      END DO                                                            !TM
        IHOUR=IHOUR+6                                                   !TM
        IF (IHOUR.GT.23) THEN                                           !TM
          IHOUR=IHOUR-24                                                !TM
          IDATE=IDATE+1                                                 !TM
          IF (IMONTH==2) THEN                                           !TM
            IF (MOD(IYEAR,4)==0) THEN                                   !TM
              IF(IDATE==30) THEN                                        !TM
                IDATE=1                                                 !TM
                IMONTH=IMONTH+1                                         !TM
              END IF                                                    !TM
            ELSE                                                        !TM
              IF(IDATE==29) THEN                                        !TM
                IDATE=1                                                 !TM
                IMONTH=IMONTH+1                                         !TM
              END IF                                                    !TM
            END IF                                                      !TM
          ELSE IF (IMONTH==4.OR.IMONTH==6.OR.IMONTH==9.OR.IMONTH==11)   !TM
     &      THEN                                                        !TM
            IF (IDATE==31) THEN                                         !TM
              IDATE=1                                                   !TM
              IMONTH=IMONTH+1                                           !TM
            END IF                                                      !TM
          ELSE                                                          !TM
            IF (IDATE==32) THEN                                         !TM
              IDATE=1                                                   !TM
              IMONTH=IMONTH+1                                           !TM
            END IF                                                      !TM
          END IF                                                        !TM
          IF (IMONTH==13) THEN                                          !TM
            IMONTH=1                                                    !TM
            IYEAR=IYEAR+1                                               !TM
          END IF                                                        !TM
        END IF                                                          !TM
        IF (IHOUT.EQ.1) THEN                                            !TM
          IF (IPOUT.EQ.1) CALL IOGRID (2) !output for every 6 hours     !TM
          CALL IOGRID (4) !Gridded data output for every 6 hours        !TM
        END IF                                                          !TM

        IF (SIXHRRUN.EQ.1) THEN                                         !TM
          CALL RESTART (2)                                              !TM
          PRINT *,'NORMAL END WITH 6-HR FCST (YEAHHHHHHH!!!!)'          !TM
          STOP 1111                                                     !TM
        END IF                                                          !TM

      ENDDO

      RETURN
      END
  
