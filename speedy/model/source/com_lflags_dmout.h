C--:
C--:  /LFLAG3/: Logical flags to include daily output 
C--:   LDOUT =  Flag to output specific daily field 
C--:            (.false. do not output)
      COMMON /LFLAG3/ LDOUT(NS3D_D+NS2D_D)
   
      LOGICAL LDOUT 
