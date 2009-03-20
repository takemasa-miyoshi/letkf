      SUBROUTINE FORDATE
C--
C--   SUBROUTINE FORDATE 
C--   
C--   Purpose :	Compute forcing fields for the current date
C--             and correction terms for horiz. diffusion
C--   Modified common blocks: DATE1, FORDAY, HDIFC4, ANOM
C--
      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      include "com_tsteps.h"
      include "com_date.h"

      include "com_forcon.h"
      include "com_dyncon0.h"
      include "com_physcon.h"

      include "com_hdifcon.h"
      include "com_forcing.h"
      include "com_anomfor.h"


      REAL GAMLAT(NLAT),
     &     CORH(NLON,NLAT), TSFC(NLON,NLAT), TREF(NLON,NLAT),
     &     PSFC(NLON,NLAT), QSFC(NLON,NLAT), QREF(NLON,NLAT)

      iitest = 0

C--   1. Define the interpolation parameters
 
C     1.1 Define the fraction of month for interpolation 
C         (middle of the month for fixed-month integrations) 

      NDAYS = 30
      IDAYH = 1+NDAYS/2

      IF (IDAY.EQ.0) THEN
        PRINT *, ' Start of FORDATE routine (init.)'
        IF (ISEASC.EQ.1) THEN
          FDAY = 0.5
        ELSE
          FDAY = 0.5*NDAYS
        ENDIF
      ELSE
        IF (ISEASC.EQ.1) THEN
          FDAY = IDAY-0.5
        ELSE
          GO TO 900
        ENDIF
      ENDIF

      FMON = FDAY/NDAYS
      TYEAR = (IMONTH-1+FMON)/12.

C--   2. Interpolate monthly-mean surface fields

C     Non-linear, mean-conserving interpolation of land/sea temperature
      CALL FORIN5 (NGP,IMONTH,FMON,SST12,SST01)
      CALL FORIN5 (NGP,IMONTH,FMON,STL12,STL01)

C     Linear interpolation of water/snow/ice fields
      CALL FORINT (NGP,IMONTH,FMON,OICE12,OICE1)
      CALL FORINT (NGP,IMONTH,FMON,SNOW12,SNOW1)
      CALL FORINT (NGP,IMONTH,FMON,SOILW12,SOILW1)

C--   3. Interpolate SST anomaly and heat fluxes 

      IF (ISSTAN.GT.0) THEN

        IF (IDAY.EQ.IDAYH) CALL NEWSST

        IF (IDAY.LT.IDAYH) THEN
          CALL FORINT (NGP,2,FMON,SSTAN2,SSTAN1)
        ELSE
          CALL FORINT (NGP,1,FMON,SSTAN2,SSTAN1)
        ENDIF

      ENDIF

      IF (IALST.GT.0.OR.IASST.GT.0.OR.IAICE.GT.0) THEN

        CALL FORIN5 (NGP,IMONTH,FMON,HFLXL12,HFLXL1)
        CALL FORIN5 (NGP,IMONTH,FMON,HFLXS12,HFLXS1)

      ENDIF

C--   4. Set maximum soil water availability to 1

      DO J=1,NLAT
        DO I=1,NLON
           SOILW1(I,J)=MIN(1.,SOILW1(I,J))
        ENDDO
      ENDDO

C--   5. Surface albedo:
C         defined as a weighed average of land and ocean albedos, where
C         land albedo depends linearly on snow depth (up to the SDALB
C         threshold) and sea albedo depends linearly on sea-ice fraction. 

      DALB=ALBICE-ALBSEA
      RSD=1./SDALB

      DO J=1,NLAT
        DO I=1,NLON
          SNOWC1(I,J)=MIN(1.,RSD*SNOW1(I,J))
          ALBL=ALB0(I,J)+MAX(ALBSN-ALB0(I,J),0.0)*SNOWC1(I,J)
          ALBS=ALBSEA+DALB*OICE1(I,J)
          ALB1(I,J)=FMASK1(I,J)*ALBL+FMASK0(I,J)*ALBS
        ENDDO
      ENDDO

C--   6. Call flow-independent parts of physical parametrizations

      IF (IDAY.EQ.0) THEN

        CALL RADSET

        CALL SFLSET (PHIS0)

      ENDIF

      CALL SOL_OZ (SOLC,TYEAR)

C--   7. Temperature correction term for horizontal diffusion

      CALL SETGAM (TYEAR,GAMLAT)

      DO J=1,NLAT
        DO I=1,NLON
          CORH(I,J)=GAMLAT(J)*PHIS0(I,J)
        ENDDO
      ENDDO

      if (iitest.gt.1.and.iday.eq.0) then
         call outest (19,PHIS0)
         call outest (19,CORH)
      endif

      CALL SPEC (CORH,TCORH)

C--   8. Humidity correction term for horizontal diffusion

      DO J=1,NLAT
        PEXP=1./(RD*GAMLAT(J))
        DO I=1,NLON
          TSFC(I,J)=FMASK1(I,J)*STL01(I,J)+FMASK0(I,J)*SST01(I,J)
          TREF(I,J)=TSFC(I,J)+CORH(I,J)
          PSFC(I,J)=(TSFC(I,J)/TREF(I,J))**PEXP
        ENDDO
      ENDDO

      CALL SHTORH (0,NGP,TREF,  1.,-1.,DUMMY,DUMMY,QREF)
      CALL SHTORH (0,NGP,TSFC,PSFC, 1.,DUMMY,DUMMY,QSFC)

      DO J=1,NLAT
        DO I=1,NLON
          CORH(I,J)=REFRH1*(QREF(I,J)-QSFC(I,J))
        ENDDO
      ENDDO

      if (iitest.gt.1.and.iday.eq.0) call outest (19,CORH)

      CALL SPEC (CORH,QCORH)

C--   9. Add sfc temp. anomaly to climatological sfc temp.

 900  CONTINUE

      CALL ADDSTA (IALST,IASST,IAICE)

C--
      RETURN
      END

      SUBROUTINE FORINT (NGP,IMON,FMON,FOR12,FOR1)  

C--   Aux. routine FORINT : linear interpolation of monthly-mean forcing

      REAL FOR12(NGP,*), FOR1(NGP)

      IF (FMON.LE.0.5) THEN
        IMON2 = IMON-1
        IF (IMON.EQ.1) IMON2 = 12
        WMON = 0.5-FMON
      ELSE
        IMON2 = IMON+1
        IF (IMON.EQ.12) IMON2 = 1
        WMON = FMON-0.5
      ENDIF

      DO J=1,NGP
        FOR1(J) = FOR12(J,IMON)+WMON*(FOR12(J,IMON2)-FOR12(J,IMON))
      ENDDO
C--
      RETURN
      END

      subroutine FORIN5 (ngp,imon,fmon,for12,for1)

C--   Aux. routine FORIN5 : non-linear, mean-conserving interpolation 
C--                         of monthly-mean forcing fields

      real for12(ngp,12), for1(ngp)

      im2 = imon-2
      im1 = imon-1
      ip1 = imon+1
      ip2 = imon+2

      if (im2.lt.1)  im2 = im2+12
      if (im1.lt.1)  im1 = im1+12
      if (ip1.gt.12) ip1 = ip1-12
      if (ip2.gt.12) ip2 = ip2-12
 
      c0 = 1./12.
      t0 = c0*fmon
      t1 = c0*(1.-fmon)
      t2 = 0.25*fmon*(1-fmon)

      wm2 =        -t1   +t2
      wm1 =  -c0 +8*t1 -6*t2
      w0  = 7*c0      +10*t2     
      wp1 =  -c0 +8*t0 -6*t2
      wp2 =        -t0   +t2 

      do j=1,ngp
        for1(j) = wm2*for12(j,im2)+wm1*for12(j,im1)
     &           + w0*for12(j,imon)
     &           +wp1*for12(j,ip1)+wp2*for12(j,ip2)
      enddo

      return
      end

      SUBROUTINE SETGAM (TYEAR,GAMLAT)  

C--   Aux. routine GAMLAT : compute reference lapse rate 
C--                         as a function of latitude and date

      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      include "com_dyncon0.h"
      include "com_physcon.h"

      REAL GAMLAT(NLAT)

      GAMLAT(1) = GAMMA/(1000.*GG)
      DO J=2,NLAT
        GAMLAT(J) = GAMLAT(1)
      ENDDO
C--
      RETURN
      END

      SUBROUTINE NEWSST  

C--   Aux. routine NEWSST : update SST anomaly field 

      include "atparam.h"

      include "com_forcing.h"
      include "com_anomfor.h"

      real*4 r4inp(ix,il)

      do j = 1,il
        do i = 1,ix
          sstan2(i,j,1) = sstan2(i,j,2)
        enddo
      enddo

      read(30,end=100) ((r4inp(i,j),i=1,ix),j=il,1,-1)

      do j = 1,il
        do i = 1,ix
          sstan2(i,j,2) = r4inp(i,j)
        enddo
      enddo

      CALL FORCHK (fmasko1,sstan2(1,1,2),ix*il,1,-50.,50.,273.)

      RETURN

 100  continue

      print *, ' WARNING: end-of-file reached on SST anomaly file'
      print *, ' SST anomaly will be kept constant'

C--
      RETURN
      END

      SUBROUTINE OUTEST (iunit,fout)

C--   Aux. routine OUTEST : write one field on a test output file 

      include "atparam.h"

      real*4 r4out(ix,il)

      do j=1,il
        do i=1,ix
          r4out(i,j)=fout(i,j)
        enddo
      enddo

      write (iunit) r4out

C--
      RETURN
      END
