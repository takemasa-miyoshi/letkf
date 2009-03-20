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

      ISTEP = 1
      NDAYS = 30

C     Open output files
      CALL SETGRD (0)

C     Loop over months

      DO JMONTH=1,NMONTS

        PRINT *, ' Start of year/month = ', IYEAR, IMONTH
        IF (JMONTH.EQ.NMONTS) NDAYS = NDAYSL

C       Loop over days

        DO JDAY=1,NDAYS

C         Modify the forcing fields according to the date

          IDAY = JDAY
          CALL FORDATE

C         Integrate the atmospheric model for 1 day

          CALL STLOOP (ISTEP)

C         Integrate the surface anomaly model (if requested)

          CALL SFC_AN (IALST,IASST,IAICE,NSTEPS)

        ENDDO

C       Update the year and month indices, and reset output files

        IMONTH=IMONTH+1
        IF (IMONTH.GT.12) THEN

          IMONTH = 1
          IYEAR = IYEAR+1

          IF (JMONTH.NE.NMONTS) CALL SETGRD (1)

        ENDIF

C       Write a restart dataset

        IF (MOD(JMONTH,NMONRS).EQ.0) CALL RESTART (2)

      ENDDO


      STOP
      END
