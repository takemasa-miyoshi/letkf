      SUBROUTINE TMINC
C--
C--   SUBROUTINE TMINC
C--
C--   Purpose : perform post-processing on pressure levels
C--             and increment time-mean arrays
C--   Modified common blocks : TMSAVE
C--

C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"
C
      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Parameters for post-processing arrays
      include "par_tmean.h"

C     Post-processing arrays (time means)
      include "com_tmean.h"
      include "com_tmean_daily.h"

C     Constants and conversion factors
      include "com_physcon.h"

C     Model variables, tendencies and fluxes on gaussian grid
      include "com_physvar.h"

C     Surface anomaly variables on gaussian grid
      include "com_anomvar.h"

C     Forcing fields
      include "com_forcing.h"

C     Logical flags
      include "com_lflags.h"
  
      real ADSAVE(ngp,6), FPQ(ngp), PHISG(ngp)
      equivalence (PHISG,PHIS0)

C     fields for vertical interpolation
      integer K0(ngp)
      real    W0(ngp), ZOUT(ngp), ZINP(nlev), RDZINP(nlev)

      data FACT2D / 1000., 1000., 2*1., 100., 3*1., 100., 100., 1000.,
     &              6*1.0 /

      iitest=0
      if (iitest.eq.1) print *, ' inside TMINC'

      rg    = 1./gg
      rdr2  = 0.5*rd
      gam0  = 0.006*rg
      rgam  = rd*gam0
      rrgam = 1./rgam

C--   1. Increment 2-d time-mean fields

      if (iitest.eq.1) print*,' store 2d fields'

c     mean-sea-level pressure
      do j=1,ngp
        tsg=0.5*(T0(j)+max(255.,min(295.,T0(j))))
        FPQ(j)=PSG(j)*(1.+gam0*PHISG(j)/tsg)**rrgam
      enddo

      n0=0
      call ADD1F (SAVE2D,PSG,       ngp,n0)
      call ADD1F (SAVE2D,FPQ,       ngp,n0)
      call ADD1F (SAVE2D,TS,        ngp,n0)
      call ADD1F (SAVE2D,TSKIN,     ngp,n0)
      call ADD1F (SAVE2D,SOILW1,    ngp,n0)
      call ADD1F (SAVE2D,U0,        ngp,n0)
      call ADD1F (SAVE2D,V0,        ngp,n0)
      call ADD1F (SAVE2D,T0,        ngp,n0)
      call ADD1F (SAVE2D,RH(1,NLEV),ngp,n0)
      call ADD1F (SAVE2D,CLOUDC,    ngp,n0)
      call ADD1F (SAVE2D,CLTOP,     ngp,n0)

      call ADD1F (SAVE2D,SHF(1,3), ngp,n0)
      call ADD1F (SAVE2D,TSR,      ngp,n0)
      call ADD1F (SAVE2D,SSR,      ngp,n0)
      call ADD1F (SAVE2D,SLR,      ngp,n0)

      call ADD1F (SAVE2D,STANOM(1,1,1),ngp,n0)
      call ADD1F (SAVE2D,SSTAN,        ngp,n0)

      n0=NS2D_D-NS2D2
      call ADD1F (SAVE2D_D,FPQ,    ngp,n0)
      call ADD1F (SAVE2D_D,T0,     ngp,n0)


C--   2. Perform vertical interpolation from sigma to pressure levels
C--      and increment 3-d time-mean fields

      if (iitest.eq.1) print*, ' store 3d fields'

      ZINP(1)  =-SIGL(1)
      do k=2,nlev
        ZINP(k)  =-SIGL(k)
        RDZINP(k)= 1./(ZINP(k-1)-ZINP(k))
C       RDZINP(k)= WVI(k-1,1)
      enddo

      zmin  = ZINP(nlev)

      do k=1,kx

C       2.1 Set coefficients for vertical interpolation
C           using coordinate Z = log (p_s/p) 

        if (lppres) then
          plog=log(POUT(k))
          do j=1,ngp
            ZOUT(j)=PSLG1(j)-plog
          enddo
        else
c         Set zout=zinp(k) to do post-proc. on sigma levels
          do j=1,ngp
            ZOUT(j)=ZINP(k)
          enddo
        endif

        call SETVIN (ZINP,RDZINP,ZOUT,ngp,kx,K0,W0)

C       2.2 Interpolate 3-d fields

c       Temperature (extrapolated below the lowest level when W0(j)<0)

        call VERINT (ADSAVE(1,2),TG1,ngp,kx,K0,W0) 

c       Remove extrapolation of temperature inversions 
c       and correct extrap. values using a reference lapse rate

        wref = 0.7

        do j=1,ngp
          if (ZOUT(j).lt.zmin) then
            textr = max(ADSAVE(j,2),TG1(j,nlev))
            aref = rgam*(zmin-ZOUT(j))
            tref = TG1(j,nlev)*(1.+aref+0.5*aref*aref)
            ADSAVE(j,2) = textr+wref*(tref-textr)
          endif
        enddo

c       Geopotential (computed from the closest levels 
c                     using the hydrostatic equation)

        do j=1,ngp
          W0(j)=max(W0(j),0.)
        enddo

        do j=1,ngp
          kj=K0(j)
          kj1=kj-1
          phi1=PHIG1(j,kj)
     &         +rdr2*(ADSAVE(j,2)+TG1(j,kj ))*(ZOUT(j)-ZINP(kj ))
          phi2=PHIG1(j,kj1)
     &         +rdr2*(ADSAVE(j,2)+TG1(j,kj1))*(ZOUT(j)-ZINP(kj1))
          ADSAVE(j,1)=phi1+W0(j)*(phi2-phi1)
        enddo

c       Wind and relative humidity 

c       a) Interpolate above the lowest level

C       call VERINT (ADSAVE(1,1),PHIG1,ngp,kx,K0,W0) 
        call VERINT (ADSAVE(1,3),UG1,  ngp,kx,K0,W0) 
        call VERINT (ADSAVE(1,4),VG1,  ngp,kx,K0,W0) 
        call VERINT (ADSAVE(1,6),RH,   ngp,kx,K0,W0) 

c       b) Decrease wind speed below the lowest level

        do j=1,ngp
          if (ZOUT(j).lt.zmin) then
            fwind=ADSAVE(j,1)/PHIG1(j,nlev)
            ADSAVE(j,3)=ADSAVE(j,3)*fwind
            ADSAVE(j,4)=ADSAVE(j,4)*fwind  
          endif
        enddo

c       Estimate specific humidity using interpolated rel.hum. and
c       sat. spec.hum. at interpolated temperature

        call SHTORH (-1,ngp,ADSAVE(1,2),POUT(k),-1.,
     *               ADSAVE(1,5),ADSAVE(1,6),FPQ)

c       Below the surface, set spec.hum. = near-surface value 

        do j=1,ngp
          if (ZOUT(j).lt.0.0) then
            ADSAVE(j,5)=Q0(j)
            ADSAVE(j,6)=Q0(j)/FPQ(j)
          endif
        enddo

c        rescale geopotential and rel. humidity

        do j=1,ngp
          ADSAVE(j,1)=ADSAVE(j,1)*rg
          ADSAVE(j,6)=ADSAVE(j,6)*100.
        enddo

C       2.3 Add 3-d fields to time-mean and daily-means arrays

        do n=1,6
         do j=1,ngp
           SAVE3D(j,k,n)=SAVE3D(j,k,n)+ADSAVE(j,n)
         enddo
        enddo

        if(k.eq.6)then
         do j=1,ngp
           SAVE3D_D(j,1)=SAVE3D_D(j,1)+ADSAVE(j,3)
           SAVE3D_D(j,2)=SAVE3D_D(j,2)+ADSAVE(j,4)
           SAVE3D_D(j,3)=SAVE3D_D(j,3)+ADSAVE(j,5)
         enddo
        elseif(k.eq.4)then
          do j=1,ngp
           SAVE3D_D(j,4)=SAVE3D_D(j,4)+ADSAVE(j,1)
          enddo
        elseif(k.eq.2)then
          do j=1,ngp
           SAVE3D_D(j,5)=SAVE3D_D(j,5)+ADSAVE(j,3)
           SAVE3D_D(j,6)=SAVE3D_D(j,6)+ADSAVE(j,4)
          enddo
        endif

C       2.4 Increment variances on pressure levels

        if (ns3d2.gt.0) then

         do n=1,4
          nv=n+ns3d1
          do j=1,ngp
            SAVE3D(j,k,nv)=SAVE3D(j,k,nv)+ADSAVE(j,n)*ADSAVE(j,n)
          enddo
         enddo

         nuv=ns3d1+5
         nvt=ns3d1+6
         do j=1,ngp
           SAVE3D(j,k,nuv)=SAVE3D(j,k,nuv)+ADSAVE(j,3)*ADSAVE(j,4)
           SAVE3D(j,k,nvt)=SAVE3D(j,k,nvt)+ADSAVE(j,2)*ADSAVE(j,4)
         enddo

        endif

C       end-of-loop over pressure levels

      enddo

C--   3. Save diabatic forcing terms on model levels

      if (ns3d3.gt.0) then

        n0=ns3d1+ns3d2
        call ADD1F (SAVE3D,TT_LSC,ngp*nlev,n0)
        call ADD1F (SAVE3D,TT_CNV,ngp*nlev,n0)
        call ADD1F (SAVE3D,TT_RSW,ngp*nlev,n0) 
        call ADD1F (SAVE3D,TT_RLW,ngp*nlev,n0) 
        call ADD1F (SAVE3D,TT_PBL,ngp*nlev,n0) 

      endif

      if (iitest.eq.1) print *, 'end of TMINC'

      RETURN
      END

      SUBROUTINE SETVIN (ZINP,RDZINP,ZOUT,NGP,NLEV,K0,W0)
C
      INTEGER K0(NGP)
      REAL    ZINP(NLEV), RDZINP(NLEV), ZOUT(NGP), W0(NGP)
C
C *** 1. Select closest vertical levels
C
      DO J=1,NGP
        K0(J)=2
      ENDDO
C
      DO K=2,NLEV-1
       DO J=1,NGP
         IF (ZOUT(J).LT.ZINP(K)) K0(J)=K+1
       ENDDO
      ENDDO
C
C *** 2. Compute interpolation weight
C
      DO J=1,NGP
        W0(J)=(ZOUT(J)-ZINP(K0(J)))*RDZINP(K0(J))
      ENDDO
C
      RETURN
      END

      SUBROUTINE VERINT (F2D,F3D,NGP,NLEV,K0,W0)

C *** 1. Perform vertical interpolation 

      INTEGER K0(NGP)
      REAL    F2D(NGP), F3D(NGP,NLEV), W0(NGP)

      DO J=1,NGP
        F2D(J)=F3D(J,K0(J))+W0(J)*(F3D(J,K0(J)-1)-F3D(J,K0(J)))
      ENDDO
C
      RETURN
      END

      SUBROUTINE ADD1F (FSAVE,FADD,NGP,NF)

C *** Add one field to storage array 

      REAL FSAVE(NGP,*), FADD(NGP)

      NF=NF+1

      DO J=1,NGP
        FSAVE(J,NF)=FSAVE(J,NF)+FADD(J)
      ENDDO
C
      RETURN
      END


