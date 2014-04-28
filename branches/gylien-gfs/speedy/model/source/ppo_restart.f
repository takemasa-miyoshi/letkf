
      SUBROUTINE RESTART (IMODE)
C--
C--   SUBROUTINE RESTART (IMODE)
C--
C--   Purpose : read or write a restart file
C--   Input :   IMODE = 1 : read model variables from a restart file
C--                   = 2 : write model variables  to a restart file
C--                   = 3 : write end-of-data flag to a restart file
C--   Initialized common blocks (if IMODE = 1) : DYNSP1, SFCANOM 
C--
      include "atparam.h"
      include "atparam1.h"

      include "com_dynvar.h"
      include "com_anomvar.h"

      include "com_date.h"
      include "com_tsteps.h"


      IF (IMODE.EQ.1) THEN

C--   1. Read the restart dataset corresponding to the specified initial date

  100    CONTINUE
         READ (3,END=200) IYEAR, IMONTH

         IF (IYEAR.EQ.IYEAR0.AND.IMONTH.EQ.IMONT0) THEN

           print*, 'Read restart dataset for year/month: ', IYEAR,IMONTH
           
           READ (3) VOR
           READ (3) DIV
           READ (3) T
           READ (3) PS
           READ (3) TR

           READ (3) STANOM

         ELSE
									
           print*, 'Skip restart dataset for year/month: ', IYEAR,IMONTH
           
           DO JREC=1,5
             READ (3) ADUMMY
           ENDDO

           READ (3) SDUMMY

           GO TO 100

         ENDIF

      ELSE IF (IMODE.EQ.2) THEN

C--   2. Write date and model variables to the restart file

         print*, 'Write restart dataset for year/month: ', IYEAR,IMONTH

         WRITE (10) IYEAR, IMONTH

         WRITE (10) VOR
         WRITE (10) DIV
         WRITE (10) T
         WRITE (10) PS
         WRITE (10) TR

         WRITE (10) STANOM

      ELSE

C--   3. Write end-of-data flag to the restart file

C         IZERO=0
C         WRITE (10) IZERO, IZERO

      ENDIF
C--
      RETURN

C--   4. Stop integration if restart file is not found

  200 CONTINUE

      print*, ' No restart dataset for the specified initial date'

      STOP 'invalid restart'

C--
      END

