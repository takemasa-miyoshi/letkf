 
      SUBROUTINE INFORC (GRAV)
C--
C--   SUBROUTINE INFORC (GRAV)
C--
C--   Purpose : Read forcing (boundary condition) fields 
C--   Input :   GRAV = gravity accel.
C--             Initialized common blocks: LSMASK, FORFIX, FORMON, ANOM
C--																		

      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL)

      include "com_tsteps.h" 
      include "com_forcon.h" 
      include "com_forcing.h"    
      include "com_anomfor.h"

      real*4 r4inp(ix,il), dummy4
      real*4 veg(ix,il), swl1(ix,il), swl2(ix,il)

      iitest=1

c     set threshold for land-sea mask definition
c     (ie minimum fraction of either land or sea)

      thrsh = 0.1

C--   1. Read time-invariant fields (orography, land-sea mask, sfc albedo)

      if (iitest.ge.1) print*,' read orography' 

      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)

      do j = 1,il
        do i = 1,ix
          phi0(i,j) = grav*r4inp(i,j)
        enddo
      enddo

      call truncg (ntrun,phi0,phis0)
 
      if (iitest.ge.1) print*,' read fractional land-sea mask'  

      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)

      do j = 1,il
        do i = 1,ix
          fmask1(i,j) = r4inp(i,j)
        enddo
      enddo

      if (iitest.ge.1) print*,' read surface albedo'  

      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)

      do j = 1,il
        do i = 1,ix
          alb0(i,j) = 0.01*r4inp(i,j)
        enddo
      enddo

C--   2. Compute additional land-sea masks
c
c     fmask0  = 1 - fmask1 = sea fraction
c     fmaskl1 = 1 wherever there is any land, 0 elsewhere
c     fmasko1 = 1 wherever there is any sea,  0 elsewhere
c
c     Note that because some points have both land and ocean
c     it is *not* true that: fmasko1+fmaskl1=1 everywhere

      do i=1,ix
       do j=1,il

         fmask0(i,j)=1.0-fmask1(i,j)

         if (fmask1(i,j).ge.thrsh) then
           fmaskl1(i,j)=1.0
         else
           fmaskl1(i,j)=0.0
           fmask1 (i,j)=0.0
           fmask0 (i,j)=1.0
         endif

         if (fmask0(i,j).ge.thrsh) then
           fmasko1(i,j)=1.0
         else
           fmasko1(i,j)=0.0
           fmask0 (i,j)=0.0
           fmask1 (i,j)=1.0
         endif

       enddo
      enddo

C--   3. Read monthly-mean climatologies of surface fields

c     3.1 Read SST 

      if (iitest.ge.1) print*,' reading sst' 

      do it = 1,12
        read (21) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
            sst12(i,j,it) = r4inp(i,j)
          enddo
        enddo
      enddo

      if (iitest.ge.1) print*,' checking sst'

      CALL FORCHK (fmasko1,sst12,ix*il,12,0.,400.,273.)

c     3.2 Read sea ice

      if (iitest.ge.1) print*,' reading sea ice'  

      do it = 1,12
        read (22) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
            oice12(i,j,it) = r4inp(i,j)
          enddo
        enddo
      enddo

      if (iitest.ge.1) print*,' checking sea ice'

      CALL FORCHK (fmasko1,oice12,ix*il,12,0.,1.,0.)


c     3.3 Read land-surface temp

      if (iitest.ge.1) print*,' reading land-surface temp.'
  
      do it = 1,12
        read (23) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
            stl12(i,j,it) = r4inp(i,j)
          enddo
        enddo
      enddo

      if (iitest.eq.1) print*,' checking land-surface temp.'

      CALL FORCHK (fmaskl1,stl12,ix*il,12,0.,400.,273.)

      do it = 1,12

        call ftland (stl12(1,1,it),phi0,phis0,fmaskl1)

        if (iitest.gt.1) then
          do j = 1,il
            do i = 1,ix
              r4inp(i,j) = stl12(i,j,it)
            enddo
          enddo
          write (18) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        endif

      enddo

c     3.4 Read snow depth

      if (iitest.ge.1) print*,' reading snow depth'  

      do it = 1,12
        read (24) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
            snow12(i,j,it) = r4inp(i,j)
          enddo
        enddo
      enddo

      if (iitest.ge.1) print*,' checking snow depth'

      CALL FORCHK (fmaskl1,snow12,ix*il,12,0.,20000.,0.)

c    3.5 Read soil moisture and compute soil water availability 
c        using vegetation fraction

      if (iitest.ge.1) print*,' reading soil moisture'  

c     read vegetation fraction (in %)
      read (25) ((veg(i,j),i=1,ix),j=il,1,-1)

      do j = 1,il
        do i = 1,ix
          veg(i,j)=max(0.,0.01*veg(i,j))
        enddo
      enddo

      sdep1 = 70.
      idep2 = 3
      sdep2 = idep2*sdep1

      swwil2= sdep2*swwil
      rsw   = 1./(sdep1*swcap+sdep2*(swcap-swwil))

      do it = 1,12
        read (26) ((swl1(i,j),i=1,ix),j=il,1,-1)
        read (26) ((swl2(i,j),i=1,ix),j=il,1,-1)
        read (26) dummy4
        do j = 1,il
          do i = 1,ix
            swroot = idep2*swl2(i,j)
            soilw12(i,j,it) = rsw*(swl1(i,j)+
     &                        veg(i,j)*max(0.,swroot-swwil2))									
          enddo
        enddo
      enddo

      if (iitest.ge.1) print*,' checking soil moisture'

      CALL FORCHK (fmaskl1,soilw12,ix*il,12,0.,10.,0.)

C--   4. Read SST anomalies for initial and preceeding month

      if (isstan.gt.0) then

        if (iitest.ge.1) print*,' reading sst anomalies' 

        do jrec=1,isst0-2
          read (30) dummy4
        enddo

        read (30) ((r4inp(i,j),i=1,ix),j=il,1,-1)

        do j = 1,il
          do i = 1,ix
            sstan2(i,j,1) = r4inp(i,j)
          enddo
        enddo

        if (isst0.gt.1) read (30) ((r4inp(i,j),i=1,ix),j=il,1,-1)

        do j = 1,il
          do i = 1,ix
            sstan2(i,j,2) = r4inp(i,j)
          enddo
        enddo

        if (iitest.ge.1) print*,' checking sst anomalies'

        CALL FORCHK (fmasko1,sstan2,ix*il,2,-50.,50.,273.)

      endif

C--   4. Read heat fluxes for surface anomaly model

      if (ialst.gt.0.or.iasst.gt.1.or.iaice.gt.0) then

        if (iitest.ge.1) print*,' reading sfc heat fluxes' 

        irecl=4*ix*il
        irec =0

        open ( unit=31, file='fort.31', status='old', 
     &         form='unformatted', access='direct', recl=irecl )

        do it = 1,12

          irec=irec+1
          read (31,rec=irec) r4inp
          do j = 1,il
            do i = 1,ix
              hflxl12(i,j,it) = r4inp(i,j)
            enddo
          enddo

          irec=irec+1
          read (31,rec=irec) r4inp
          do j = 1,il
            do i = 1,ix
              hflxs12(i,j,it) = r4inp(i,j)
            enddo
          enddo  

        enddo

        if (iitest.ge.1) print*,' checking sfc heat fluxes'

        CALL FORCHK (fmaskl1,hflxl12,ix*il,12,-1000.,1000.,0.)
        CALL FORCHK (fmasko1,hflxs12,ix*il,12,-1000.,1000.,0.)

      endif
C--
      RETURN
      END

      SUBROUTINE FORCHK (FMASK,FIELD,NGP,NF,FMIN,FMAX,FSET)
										
C--   Aux. routine FORCHK: Check consistency of sfc fields with land-sea mask 
C--   and set undefined values to a constant (to avoid over/underflow)

      real fmask(ngp), field(ngp,nf)

      do jf = 1,nf

        nfault=0

        do jgp = 1,ngp
          if (fmask(jgp).gt.0.0) then
            if (field(jgp,jf).lt.fmin.or.field(jgp,jf).gt.fmax)
     *          nfault = nfault+1
          else
            field(jgp,jf) = fset
          endif
        enddo

        print *, ' field: ', jf, '   no. of faulty points:', nfault

      enddo

      print *, ' undefined values set to', fset

      RETURN
      END

      SUBROUTINE FTLAND (STL,PHI0,PHIS0,FMASKL)

      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL)

      include "com_dyncon0.h" 
      include "com_dyncon1.h" 

      REAL STL(NLON,NLAT), PHI0(NLON,NLAT), PHIS0(NLON,NLAT),
     &     FMASKL(NLON,NLAT)

      REAL STL2(NLON,NLAT)

      NL8 = NLAT/8
      GAM = 0.001*GAMMA/GRAV

      NLAT1 = 1
      NLAT2 = NL8

      DO JBAND=1,8

        SUMT=0.
        SUMW=0.

        DO J=NLAT1,NLAT2
         DO I=1,NLON
           STL(I,J)=STL(I,J)+GAM*PHI0(I,J)
           SUMT=SUMT+GCOS(J)*FMASKL(I,J)*STL(I,J)
           SUMW=SUMW+GCOS(J)*FMASKL(I,J)
         ENDDO
        ENDDO

        SUMT=SUMT/SUMW

        DO J=NLAT1,NLAT2
         DO I=1,NLON
           IF (FMASKL(I,J).EQ.0.) STL(I,J)=SUMT
         ENDDO
        ENDDO
  
        NLAT1=NLAT1+NL8
        NLAT2=NLAT2+NL8

      ENDDO

      ITR=7
      IDTR=(NTRUN-6)/3

      DO JFIL=1,4

        CALL TRUNCG (ITR,STL,STL2)

        DO J=1,NLAT
         DO I=1,NLON
           IF (FMASKL(I,J).EQ.0.) STL(I,J)=STL2(I,J)
         ENDDO
        ENDDO

        ITR=MIN(ITR+IDTR,NTRUN)

      ENDDO

      CALL TRUNCG (ITR,STL,STL2)

      DO J=1,NLAT
       DO I=1,NLON
         STL(I,J)=STL2(I,J)-GAM*PHIS0(I,J)
       ENDDO
      ENDDO       

      RETURN
      END

      SUBROUTINE TRUNCG (ITR,FG1,FG2)

C--   SUBROUTINE TRUNCG (ITR,FG1,FG2)
C--   Purpose : compute a spectrally-filtered grid-point field
C--   Input   : ITR : spectral truncation (triangular)
C--           : FG1 : original grid-point field
C--   Output  : FG2 : filtered grid-point field

      include "atparam.h"

      REAL FG1 (IX,IL), FG2(IX,IL)
      COMPLEX FSP(MX,NX), ZERO 

      ZERO = (0.,0.)

      CALL SPEC (FG1,FSP)

      DO N=1,NX
        DO M=1,MX
          ITWN=ISC*(M-1)+N-1
          IF (ITWN.GT.ITR) FSP(M,N)=ZERO
        ENDDO
      ENDDO

      CALL GRID (FSP,FG2,1)

      RETURN
      END
