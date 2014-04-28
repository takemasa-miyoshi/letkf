C--
C--:  /RADCON/: Radiation and cloud constants (initial. in INPHYS)
C--:   RHCL1  = relative hum. threshold corr. to cloud cover = 0
C--:   RHCL2  = relative hum. corr. to cloud cover = 1
C--:   QACL   = specific hum. threshold for cloud cover
C--:   WPCL   = cloud c. weight for the sq. root of precip. (for p = 1 mm/day)
C--:   PMAXCL = max. value of precip. (mm/day) contributing to cloud cover 	 
C--:   ALBCL  = cloud albedo (for cloud cover = 1)
C--:   ALBMIN = min. albedo for stratocumulus clouds correction over sea
C--:   EPSSW  = fraction of incoming solar radiation absorbed by ozone
C--:   EPSLW  = fraction of blackbody spectrum absorbed/emitted by PBL only
C--:   EMISFC = longwave surface emissivity
C--:            shortwave absorptivities (for dp = 10^5 Pa) :
C--:   ABSDRY = abs. of dry air      (visible band)
C--:   ABSAER = abs. of aerosols     (visible band)
C--:   ABSWV1 = abs. of water vapour (visible band, for dq = 1 g/kg)
C--:   ABSWV2 = abs. of water vapour (near IR band, for dq = 1 g/kg)
C--:   ABSCL2 = abs. of clouds       (visible band, for dq_base = 1 g/kg)
C--:   ABSCL1 = abs. of clouds       (visible band, maximum value)
C--:            longwave absorptivities (per dp = 10^5 Pa) :
C--:   ABLWIN = abs. of air in "window" band
C--:   ABLCO2 = abs. of air in CO2 band
C--:   ABLWV1 = abs. of water vapour in H2O band 1 (weak),   for dq = 1 g/kg
C--:   ABLWV2 = abs. of water vapour in H2O band 2 (strong), for dq = 1 g/kg
C--:   ABLCL1 = abs. of "thick" clouds in window band (below cloud top) 
C--:   ABLCL2 = abs. of "thin" upper clouds in window and H2O bands

      COMMON /RADCON/ RHCL1,  RHCL2,  QACL,   WPCL,   PMAXCL, 
     &                ALBCL, ALBMIN, EPSSW,  EPSLW,  EMISFC, 
     &                ABSDRY, ABSAER, ABSWV1, ABSWV2, ABSCL1, ABSCL2, 
     &                ABLWIN, ABLCO2, ABLWV1, ABLWV2, ABLCL1, ABLCL2

C--								
C--   /RADFIX/: Time-invariant fields (initial. in RADSET)
C--    FBAND  = energy fraction emitted in each LW band = f(T)

      COMMON /RADFIX/ FBAND(100:400,0:4)
C--								
C--   /RADZON/: Zonally-averaged fields for SW scheme
C--             (updated in SOL_OZ)
C--    FSOL   = flux of incoming solar radiation
C--    OZONE  = flux absorbed by ozone (lower stratos.)
C--    OZUPP  = flux absorbed by ozone (upper stratos.)
C--    ZENIT  = optical depth ratio (function of solar zenith angle) 

      COMMON /RADZON/ FSOL(NGP), OZONE(NGP), OZUPP(NGP), ZENIT(NGP), 
     &                STRATZ(NGP)
C--								
C--   /RADFLD/: Transmissivity and blackbody rad. (updated in RADSW/RADLW)
C--    TAU2   = transmissivity of atmospheric layers
C--    ST4A   = blackbody emission from full and half atmospheric levels
C--    STRATC = stratospheric correction term 
C--    FLUX   = radiative flux in different spectral bands

      COMMON /RADFLD/ TAU2(NGP,NLEV,4), ST4A(NGP,NLEV,2),
     &                STRATC(NGP), FLUX(NGP,4)
C--								
C--   /RADCLD/: Radiative properties of clouds (updated in CLOUD)
C--    QCLOUD = Equivalent specific humidity of clouds 

      COMMON /RADCLD/ QCLOUD(NGP)
