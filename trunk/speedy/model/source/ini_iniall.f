      SUBROUTINE INIALL
C--
C--   SUBROUTINE INIALL
C--
C--   Purpose : Call initializion routines for all model common blocks 
C--
      include "atparam.h"
      include "atparam1.h"
      include "par_tmean.h"

      include "com_tsteps.h"
      include "com_dyncon1.h"

      include "com_outfiles.h"

      REAL PPL(KX)

      iitest=1

C--   1. Initialize ffts

      if (iitest.eq.1) print*, 'calling INIFFT'
      CALL INIFFT

C--   2. Initialize constants for time-stepping and dynamics

C     2.1 Check if the run is a restart run:  ISTART=1 restart
C                                             ISTART=0 no restart
      READ (2,*) ISTART
      IF (ISTART.NE.0) ISTART = 1

C     2.2 Initialize constants and operators
 
      if (iitest.eq.1) print*, 'calling INDYNS'
      CALL INDYNS

C     2.3 Set post-processing levels

      DO K=1,KX
        PPL(K)=PRLEV(FSG(K))
      ENDDO

C--   3. Initialize constants for physical parametrization

      if (iitest.eq.1) print*, 'calling INPHYS'
      CALL INPHYS (HSG,PPL,RADANG)

C--   4. Initialize forcing fields (boundary cond. + random forcing)

      if (iitest.eq.1) print*, 'calling INFORC'
      CALL INFORC (GRAV)

      if (iitest.eq.1) print*, 'calling INIRDF'
      CALL INIRDF (INDRDF)

C--   5. Initialize model variables (atmosphere + surface)

      if (iitest.eq.1) print*, 'calling INVARS'
      CALL INVARS

      if (iitest.eq.1) print*, 'calling SFC_IN'
      CALL SFC_IN

C--   6. Initialize time-mean arrays

      if (iitest.eq.1) print*, 'calling TMOUT'

      CALL TMOUT (0)
      CALL TMOUT_DAILY (0)
 
C--   7. Set up the time-mean and daily-mean output (grads format)

      READ(2,'(A3)') NORUN

      NTM = ((NMONTS-1)*30+NDAYSL)*NSTEPS/NSTOUT
      NDM = ((NMONTS-1)*30+NDAYSL)
      NDAYTM = NSTOUT/NSTEPS
      NDAYDM = 1

      if (iitest.eq.1) print*, 'calling SETCTL'

      IS3D=1
      CALL SETCTL (12,IX,IL,KX,NTM,NDAYTM,IS3D,NS3D1,NS2D,NS2D_D-NS2D2,
     *             RADANG,PPL,'attm',NORUN,IYEAR0,IMONT0)

C--  Daily-mean
      if(IDOUT .gt. 0 .and. IDOUT .le. 3 ) then
         CALL SETCTL_DAILY (18,IX,IL,KX,NDM,NDAYDM,IS3D,NS3D_D,NS2D_D,
     *                      RADANG,PPL,'daytm',NORUN,IYEAR0,IMONT0)
      endif
C--
      IS3D=IS3D+NS3D1
      CALL SETCTL (14,IX,IL,KX,NTM,NDAYTM,IS3D,NS3D2,0,0,
     *             RADANG,PPL,'atva',NORUN,IYEAR0,IMONT0)

      IS3D=IS3D+NS3D2
      CALL SETCTL (16,IX,IL,KX,NTM,NDAYTM,IS3D,NS3D3,0,0,
     *             RADANG,PPL,'atdf',NORUN,IYEAR0,IMONT0)

C--
      RETURN
      END
  

      FUNCTION PRLEV (SIGLEV)
C--									
C--   FUNCTION PRLEV (SIGLEV)
C--   Purpose : select the closest standard pressure level for post-proc.
C--   Input :   SIGLEV = sigma level

c      REAL PLEV(16)

c      DATA PLEV/ 0.925, 0.85, 0.775, 0.7, 0.6, 0.5, 0.4, 0.3,  
c     *           0.25, 0.2, 0.15, 0.10, 0.07, 0.05, 0.03, 0.01/

      REAL PLEV(12)

      DATA PLEV/ 0.925, 0.850, 0.775, 0.700, 0.600, 0.500, 0.400,  
     *           0.300, 0.250, 0.200, 0.150, 0.100/

      PRLEV = 1.
      DIF = 1.-SIGLEV

      DO K=1,16
        ADIF = ABS(PLEV(K)-SIGLEV)
        IF (ADIF.LE.DIF) THEN
          DIF = ADIF
          PRLEV = PLEV(K)
        ENDIF
      ENDDO

      RETURN
      END



      

