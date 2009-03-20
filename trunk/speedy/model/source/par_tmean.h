C--
C--   Post-processing parameters (time-mean and variance fields) 
C--   NS2D   : no. of 2-d fields saved every t-mean
C--   NS2D_D : no. of 2-d fields saved daily 
C--   NS2D2: no. of 2-d fields saved daily with means incremented every NSTPPR
C--   NS3D1  : no. of 3-d means of model variables
C--   NS3D_D : no. of 3-d daily means of model variables
C--   NS3D2  : no. of 3-d variances and covariances
C--   NS3D3  : no. of 3-d means of diabatic forcing fields

      PARAMETER ( NS2D=17, NS2D_D=10, NS2D2=2,
     &            NS3D1=8, NS3D_D= 6, NS3D2=6, NS3D3=5 )
      PARAMETER ( NS3D=NS3D1+NS3D2+NS3D3 )
