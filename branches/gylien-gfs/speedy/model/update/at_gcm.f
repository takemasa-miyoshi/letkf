      PROGRAM AT_GCM
C--
C--   Main program : AT_GCM
C--
C--   Modified common blocks: DATE1
C--
      include "atparam.h"
      include "atparam1.h"

      include "com_tsteps.h"
      include "com_date.h"

      PRINT *, ' Hallo from AT_GCM'

C--   1. Initialization of all model constants and variables 

      CALL INIALL

C     set up the forcing fields for the first time step

      IDAY = 0
      CALL FORDATE 

C     do the initial (2nd-order) time step, initialize the semi-impl. scheme

      CALL STEPONE

C     Re-initialize the surface-flux anomalies  

      CALL SFC_AN (0,0,0,NSTEPS)

C--   2. Time integration

!     For less than one-day integration                                 !TM
      SIXHRRUN=0                                                        !TM
      IF (NDAYSL.EQ.0) THEN                                             !TM
        NDAYSL=1                                                        !TM
        SIXHRRUN=1                                                      !TM
        IHOUT=1                                                         !TM
      END IF                                                            !TM
!     Salinity check                                                    !TM
      IF (IHOUT.EQ.1.AND.NMONTS.GE.4) THEN                              !TM
        PRINT *,'You are going to make 6-hourly output for more than',  !TM
     &   ' 4 months! Check the NMONTS parameter! (AT_GCM)'              !TM
        STOP                                                            !TM
      END IF                                                            !TM

      ISTEP = 1
!     NDAYS = 30                                                        !TM
      NDAYS = 31                                                        !TM

C     Open output files
!     CALL SETGRD (0)                                                   !TM
      IF (IHOUT.NE.1) CALL SETGRD (0)                                   !TM
! Write initial data                                                    !TM
      IF (IHOUT.EQ.1.AND.IPOUT.EQ.1) CALL IOGRID (2)                    !TM
      IF (IHOUT.EQ.1) CALL IOGRID (4)                                   !TM

C     Loop over months

      DO JMONTH=1,NMONTS

!       PRINT *, ' Start of year/month = ', IYEAR, IMONTH               !TM
        PRINT '(A,I4.4,A,I2.2,A,I2.2,A,I2.2)',                          !TM
     &    ' Start of year/month/date/hour = ',                          !TM
     &    IYEAR,'/',IMONTH,'/',IDATE,'/',IHOUR                          !TM
        IF (JMONTH.EQ.NMONTS) NDAYS = NDAYSL

C       Loop over days

        DO JDAY=1,NDAYS

C         Modify the forcing fields according to the date

!         IDAY = JDAY                                                   !TM
          IDAY = IDATE                                                  !TM
          IF (IDAY.GT.30) IDAY=30                                       !TM
          CALL FORDATE

C         Integrate the atmospheric model for 1 day

          CALL STLOOP (ISTEP)

C         Integrate the surface anomaly model (if requested)

          CALL SFC_AN (IALST,IASST,IAICE,NSTEPS)

          IF (NDAYS.EQ.31.AND.IDATE.EQ.1) EXIT                          !TM

        ENDDO

C       Update the year and month indices, and reset output files
!       Time-index updates have been done in the STLOOP                 !TM
!                                                                       !TM
!       IMONTH=IMONTH+1                                                 !TM
!       IF (IMONTH.GT.12) THEN                                          !TM
!                                                                       !TM
!         IMONTH = 1                                                    !TM
!         IYEAR = IYEAR+1                                               !TM
!                                                                       !TM
!         IF (JMONTH.NE.NMONTS) CALL SETGRD (1)                         !TM
        IF (IHOUT.NE.1.AND.IMONTH.EQ.1.AND.JMONTH.NE.NMONTS.AND.        !TM
     &    IYEAR.NE.IYEAR0) CALL SETGRD (1)                              !TM
!                                                                       !TM
!       ENDIF                                                           !TM

C       Write a restart dataset

!       IF (MOD(JMONTH,NMONRS).EQ.0) CALL RESTART (2)                   !TM

      ENDDO

!     Restart dataset is written only at last                           !TM
      CALL RESTART (2)                                                  !TM
        PRINT '(A,I4.4,A,I2.2,A,I2.2,A,I2.2)',                          !TM
     &    ' End of year/month/date/hour = ',                            !TM
     &    IYEAR,'/',IMONTH,'/',IDATE,'/',IHOUR                          !TM

      STOP
      END
