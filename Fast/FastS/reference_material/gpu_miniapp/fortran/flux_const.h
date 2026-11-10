c     Constantes du flux de Roe (identiques a fluroe_euler_o3_3dfull)
      gam    = param_real( GAMMA )
      rgp    = param_real( CVINF )*(gam-1.)
      gam1   = gam/(gam-1.)
      psiroe = param_real( PSIROE )
      c6     = 1./6.
      c4     = 5.*c6
      c5     = 2.*c6
      c6     =-1.*c6
      v1 = 0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)
      v5 = 4*param_int(NDIMDX)
      v1mtr =   0
      v2mtr =   param_int(NDIMDX_MTR)
      v3mtr = 2*param_int(NDIMDX_MTR)
