      SUBROUTINE TMINC_DAILY
C--
C--   SUBROUTINE TMINC_DAILY
C--
C--   Purpose : perform increment daily-mean arrays
C--   Modified common blocks : TMSAVE_D
C--

C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Parameters for post-processing arrays
      include "par_tmean.h"

C     Post-processing arrays (time means)
      include "com_tmean_daily.h"

C     Constants and conversion factors
      include "com_physcon.h"

C     Model variables, tendencies and fluxes on gaussian grid
      include "com_physvar.h"

      data FACT2D_D / 3*86.400, 5*1.0, 1000., 1. /

      iitest=0
      if (iitest.eq.1) print *, ' inside TMINC_DAILY'
      
C--   1. Increment 2-d daily-mean fields

      if (iitest.eq.1) print*,' store 2d fields'

      n0=1

c     nb:  prec. and evap. are in (g/m^2)/s

      call ADD1F_D (SAVE2D_D,PRECLS,   ngp,n0)

      n0=1      
 
c     nb: prec. is saved daily considering PRECLS and PRECNV together

      call ADD1F_D (SAVE2D_D,PRECNV,   ngp,n0)
      call ADD1F_D (SAVE2D_D,EVAP(1,3),ngp,n0)

      call ADD1F_D (SAVE2D_D,OLR,      ngp,n0)

      call ADD1F_D (SAVE2D_D,USTR(1,3),  ngp,n0)
      call ADD1F_D (SAVE2D_D,VSTR(1,3),  ngp,n0)

      call ADD1F_D (SAVE2D_D,HFLUXN(1,1),ngp,n0)
      call ADD1F_D (SAVE2D_D,HFLUXN(1,2),ngp,n0)

      n0=0

      call ADD1F_D (SAVE2D_L,PRECLS,   ngp,n0)
      call ADD1F_D (SAVE2D_L,PRECNV,   ngp,n0)
      call ADD1F_D (SAVE2D_L,EVAP(1,3),ngp,n0)

      call ADD1F_D (SAVE2D_L,OLR,      ngp,n0)

      call ADD1F_D (SAVE2D_L,USTR(1,3),  ngp,n0)
      call ADD1F_D (SAVE2D_L,VSTR(1,3),  ngp,n0)

      call ADD1F_D (SAVE2D_L,HFLUXN(1,1),ngp,n0)
      call ADD1F_D (SAVE2D_L,HFLUXN(1,2),ngp,n0)

      if (iitest.eq.1) print *, 'end of TMINC_DAILY'

      RETURN
      END

      SUBROUTINE ADD1F_D (FSAVE,FADD,NGP,NF)

C *** Add one field to storage array 

      REAL FSAVE(NGP,4), FADD(NGP)

      NF=NF+1

      DO J=1,NGP
        FSAVE(J,NF)=FSAVE(J,NF)+FADD(J)
      ENDDO
C
      RETURN
      END

