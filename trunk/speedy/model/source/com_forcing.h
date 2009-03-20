C--
C--   /LSMASK/ land-sea masks (initial. in INFORC)
      common /LSMASK/ fmask1(ix,il),  fmask0(ix,il),
     .                fmasko1(ix,il), fmaskl1(ix,il)
C--									
C--   /FORFIX/ Time invariant forcing fields 
C--            (initial. in INFORC, except for phis0 initial. in INVARS)
      common /FORFIX/ phi0(ix,il),
     .                phis0(ix,il),
     .                alb0(ix,il)
C--
C--   /FORMON/ Monthly-mean forcing fields (initial. in INFORC)
      common /FORMON/ sst12(ix,il,12),
     .                oice12(ix,il,12),
     .                stl12(ix,il,12),
     .                snow12(ix,il,12),
     .                soilw12(ix,il,12)
C--
C--   /FORDAY/ Daily forcing fields (updated in FORDATE)
      common /FORDAY/ sst1(ix,il), sst01(ix,il),
     .                stl1(ix,il), stl01(ix,il),
     .                oice1(ix,il),
     .                snow1(ix,il), snowc1(ix,il),
     .                soilw1(ix,il),
     .                alb1(ix,il)
