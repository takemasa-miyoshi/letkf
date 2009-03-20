C--
C--   Constants for physical parametrization routines:

C--   Constants for boundary forcing (common FORCON):

      SOLC   = 342.

      ALBSEA = 0.07
      ALBICE = 0.60
      ALBSN  = 0.60
      SDALB  = 60.

      SWCAP = 0.30
      SWWIL = 0.17

C--   Constants for convection (common CNVCON):

      PSMIN  =  0.8
      TRCNV  =  6.
      RHBL   =  0.9
      RHIL   =  0.7
      ENTMAX =  0.5
      SMF    =  0.8

C--   Constants for large-scale condensation (common LSCCON):

      TRLSC  = 4.
      RHLSC  = 0.9
      DRHLSC = 0.1
      RHBLSC = 0.95

C--   Constants for radiation (common RADCON):

      RHCL1  =  0.30
      RHCL2  =  1.00
      QACL   =  0.20
      WPCL   =  0.2
      PMAXCL = 10.0

      ALBCL  =  0.35
      ALBMIN = ALBCL*0.3

      EPSSW  =  0.015
      EPSLW  =  0.05
      EMISFC =  0.98

      ABSDRY =  0.033
      ABSAER =  0.033
      ABSWV1 =  0.022
      ABSWV2 = 15.000

      ABSCL1 =  0.015
      ABSCL2 =  0.15

      ABLWIN =  0.3
      ABLCO2 =  4.0
      ABLWV1 =  0.7
      ABLWV2 = 50.0

      ABLCL1 = 12.0
      ABLCL2 =  0.6

C--   Constants for surface fluxes (common SFLCON):

      FWIND0 = 0.6
      FTEMP0 = 1.
      FHUM0  = 0.

      CDL = 2.4e-3
      CDS = 0.8e-3
      CHL = 1.2e-3
      CHS = 0.8e-3

      VGUST  = 5.
      CTDAY  = 1.0e-2
      DTHETA = 3.
      FSTAB  = 0.67
      HDRAG  = 2000.
      FHDRAG = 0.5

C--   Constants for vertical diffusion and sh. conv. (common VDICON):

      TRSHC  =  6.
      TRVDI  = 24.
      TRVDS  =  6.
      REDSHC =  0.5
      RHGRAD =  0.5
      SEGRAD =  0.1
