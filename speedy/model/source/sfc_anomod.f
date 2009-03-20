
      subroutine addflx (hfluxn)

      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL )

      include "com_forcing.h"
      include "com_anomvar.h"

      real hfluxn (nlon,nlat,2)

      do j=1,nlat
        do i=1,nlon
          hfluxn(i,j,1)=hfluxn(i,j,1)*fmask1(i,j)
          hfluxn(i,j,2)=hfluxn(i,j,2)*fmask0(i,j)
        enddo
      enddo

      do ks=1,2
        do j=1,nlat
          do i=1,nlon
            hfint(i,j,ks)=hfint(i,j,ks)+hfluxn(i,j,ks)
          enddo
        enddo
      enddo

      return
      end


      subroutine addsta (ialst,iasst,iaice)

      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL )

      include "com_forcing.h"
      include "com_anomfor.h"
      include "com_anomvar.h"

C--   1. Define SST anomaly from prescribed and/or mixed-layer values

      if (iasst.eq.1.or.iasst.eq.3) then

        sstice=273.2
        do j=1,nlat
          do i=1,nlon
            if (sst01(i,j).le.sstice) sstan1(i,j)=0.
          enddo
        enddo

      endif

      if (iasst.eq.1) then

        do j=1,nlat
          do i=1,nlon
            sstan(i,j)=sstan1(i,j)
          enddo
        enddo

      else if (iasst.eq.2) then

        do j=1,nlat
          do i=1,nlon
            sstan(i,j)=stanom(i,j,2)
          enddo
        enddo

      else if (iasst.eq.3) then

        do j=1,nlat
          do i=1,nlon
            sstan(i,j)=stanom(i,j,2)+
     &                 wobsst(i,j)*(sstan1(i,j)-stanom(i,j,2))
          enddo
        enddo

      endif

C--   2. Correct SST using sea-ice temperature anomaly if requested

      if (iaice.gt.0) then

        do j=1,nlat
          do i=1,nlon
            sstan(i,j)=sstan(i,j)+
     &                 oice1(i,j)*(stanom(i,j,3)-sstan(i,j))
          enddo
        enddo

      endif

C--   3. Superimpose anomalies to land and sea temperature climatology

      do j=1,nlat
        do i=1,nlon
          stl1(i,j)=stl01(i,j)+stanom(i,j,1)
          sst1(i,j)=sst01(i,j)+sstan(i,j)
        enddo
      enddo

      return
      end


      subroutine sfc_an (ialst,iasst,iaice,nstep)

      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL )

      include "com_forcing.h"
      include "com_anomfor.h"
      include "com_anomvar.h"

      real hfanom(nlon,nlat)

C--   1. Compute sfc temp. anomalies from flux anomalies

      rnstep = 1./nstep

C--   1.1 Land

      if (ialst.gt.0) then

        do j=1,nlat
          do i=1,nlon
            if (fmaskl1(i,j).gt.0.) then
              fsnow=1.-0.5*snowc1(i,j)
              hfanom(i,j)=fsnow*(hfint(i,j,1)*rnstep-hflxl1(i,j))
              stanom(i,j,1)=stdis(1)*
     &                      (stanom(i,j,1)+hfanom(i,j)*rhcap2(i,j,1))
            endif
          enddo
        enddo

      endif

C--   1.2 Sea-ice

      do j=1,nlat
        do i=1,nlon
          hfanom(i,j)=hfint(i,j,2)*rnstep-hflxs1(i,j)
        enddo
      enddo

      if (iaice.gt.0) then

        do j=1,nlat
          do i=1,nlon
            if (oice1(i,j).gt.0.) then
              stanom(i,j,3)=stdis(3)*
     &                      (stanom(i,j,3)+hfanom(i,j)*rhcap(3))
              hfice=flxice*stanom(i,j,3)-hfanom(i,j)
              hfanom(i,j)=hfanom(i,j)+oice1(i,j)*hfice
            endif
          enddo
        enddo

      else

        do j=1,nlat
          do i=1,nlon
            hfanom(i,j)=hfanom(i,j)*(1.-oice1(i,j))
          enddo
        enddo

      endif

C--   1.3 Ocean mixed layer

      if (iasst.gt.1) then

        do j=1,nlat
          do i=1,nlon
            if (rhcap2(i,j,2).gt.0.) then
              stanom(i,j,2)=stdis(2)*
     &                      (stanom(i,j,2)+hfanom(i,j)*rhcap2(i,j,2))
            endif
          enddo
        enddo

      endif

C--   2. Set flux integral to zero for next step 

      do ks=1,2
        do j=1,nlat
          do i=1,nlon
            hfint(i,j,ks)=0.
          enddo
        enddo
      enddo

      return
      end


      subroutine sfc_in 

      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL )

      include "com_tsteps.h"
      include "com_dyncon1.h"

      include "com_forcing.h"
      include "com_anomvar.h"

C--   1. Set heat capacities and dissipation times
C--      for land (root-layer), sea (mixed layer) and sea-ice 

C     organic soil layer (depth*heat_cap/m)
C      hcapl =  0.6*2.50e+6
      hcapl =  1.*2.50e+6
C     oceanic mixed layer (depth*heat_cap/m)
      hcaps = 50.0*4.18e+6
C     sea-ice layer (depth*heat_cap/m)
      hcapi =  1.8*1.93e+6

C     Dissipation time (days) for land-surface temperature anomalies
C      tdlsta  = 20.
      tdlsta  = 40.
C     Dissipation time (days) for sea-surface temperature anomalies
      tdssta  = 90.
C     Dissipation time (days) for sea-ice temperature anomalies
C      tdiceta = 10.
      tdiceta = 20.

      rhcap(1) = 86400./hcapl
      rhcap(2) = 86400./hcaps
      rhcap(3) = 86400./hcapi

      stdis(1) = tdlsta/(tdlsta+1.)
      stdis(2) = tdssta/(tdssta+1.)
      stdis(3) = tdiceta/(tdiceta+1.)

      flxice = hcapi/(86400.*tdiceta)

C--   2. Set variable mixed-layer depth over selected regions

      nlat2=nlat/2

c     Continents/land-ice
      do j=1,nlat
        do i=1,nlon
           if (alb0(i,j).lt.0.4) then
             rhcap2(i,j,1) = rhcap(1)
           else
C             rhcap2(i,j,1) = 0.5*rhcap(3)
             rhcap2(i,j,1) = rhcap(3)
           endif
        enddo
      enddo

c     Oceans
      do j=1,nlat
        do i=1,nlon
          rhcap2(i,j,2) = -1.
        enddo
      enddo

c     Atlantic Ocean
      dlon=360./nlon
      reast=42.
      rwest=300.-dlon
      do j=nlat2+1,nlat
        do i=1,nlon
          rlon=(i-1)*dlon
          if (rlon.lt.reast.or.rlon.gt.rwest) rhcap2(i,j,2) = rhcap(2)
        enddo
        rwest=max(260.,rwest-2*dlon)
      enddo

      do i=1,nlon
        rhcap2(i,nlat2+1,2)=rhcap2(i,nlat2+1,2)*0.25
        rhcap2(i,nlat2+2,2)=rhcap2(i,nlat2+2,2)*0.5
      enddo

C--   3. Set weight mask for observed SST anomalies

      do j=1,nlat
        do i=1,nlon
          wobsst(i,j) = 0.
        enddo
      enddo

c     Tropics ( weight = sqrt(cos(3*lat)) for abs(lat) < 30 )
c     Indian + Pacific only

      nwest=1+nint(nlon/16.)
      neast=1+nint(nlon*290./360.)
      rad30=asin(0.5)

      do j=1,nlat
        if (abs(radang(j)).lt.rad30) then
          wob1 = sqrt(cos(3.*radang(j)))
          do i=nwest,neast
            if (rhcap2(i,j,2).lt.0.) wobsst(i,j) = wob1
          enddo
        endif
      enddo

C--   4. Initialize anomaly model variables

      if (ialst.eq.0.or.istart.eq.0) then
        do j=1,nlat
          do i=1,nlon
            stanom(i,j,1)=0.
          enddo
        enddo
      endif

      if (iasst.eq.0.or.istart.eq.0) then
        do j=1,nlat
          do i=1,nlon
            stanom(i,j,2)=0.
          enddo
        enddo
      endif

      if (iaice.eq.0.or.istart.eq.0) then
        do j=1,nlat
          do i=1,nlon
            stanom(i,j,3)=0.
          enddo
        enddo
      endif

      do j=1,nlat
        do i=1,nlon
          sstan(i,j)=0.
        enddo
      enddo

      do ks=1,2
        do j=1,nlat
          do i=1,nlon
            hfint(i,j,ks)=0.
          enddo
        enddo
      enddo

      return
      end
