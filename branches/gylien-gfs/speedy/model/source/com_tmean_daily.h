C--
C--   /TMSAVE_D/ : Arrays for the computation of daily-means calculated
C--               at every time-step
C--              (initial./updated in TMINC_DAILY and TMOUT_DAILY)
      COMMON /TMSAVE_D/ SAVE2D_D(NGP,NS2D_D), SAVE3D_D(NGP,NS3D_D),
     &		      SAVE2D_L(NGP,NS2D_D-NS2D2), FACT2D_D(NS2D_D)
     
