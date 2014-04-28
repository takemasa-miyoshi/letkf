      SUBROUTINE INIRDF (INDRDF)
C--
C--   SUBROUTINE INIRDF (INDRDF)
C--
C--   Purpose : Initialize random diabatic forcing 
C--   Input :   INIRDF = index of forcing perturbation
C--   Initialized common blocks: RANDF
C--
      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

      include "com_physcon.h"

      include "com_randfor.h"

      real randf1(ngp), randf2(nlon,nlat)
      equivalence (randf1, randf2)

c     Test forcing: Upper tropospheric heating at lon=180, lat=0

      do j=1,ngp
        randf1(j) = 0.
      enddo

      randf2(nlon/2  ,nlat/2  ) = 1.
      randf2(nlon/2+1,nlat/2  ) = 1.
      randf2(nlon/2  ,nlat/2+1) = 1.
      randf2(nlon/2+1,nlat/2+1) = 1.

      do j=1,ngp
        randfh(j,1) = randf1(j)
        randfh(j,2) = 0.
      enddo

      if (indrdf.eq.0) then
         ampl = 0.
      else if (indrdf.gt.0) then
         ampl = 4.
      else if (indrdf.lt.0) then
         ampl =-4.
      endif

      ampl = ampl/86400.

      randfv(1,1) = 0.
      randfv(2,1) = ampl*0.3
      randfv(3,1) = ampl
      randfv(4,1) = ampl*0.9
      randfv(5,1) = ampl*0.6
      randfv(6,1) = ampl*0.3
      randfv(7,1) = 0.

      do k=1,nlev
        randfv(k,2) = 0.
      enddo

c     if (indrdf.eq.0) return

      return
      end

