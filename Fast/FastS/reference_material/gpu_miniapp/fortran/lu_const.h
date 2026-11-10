c     Constantes LU-SGS (identiques a invlu_l / invlu_u)
      gam1    = param_real(GAMMA)
      gamm1   = gam1 - 1.
      cp      = param_real(GAMMA)*param_real(CVINF)
      gam2    = gamm1-1.
      cvinv   = 1./param_real(CVINF)
      cvinv2  = 0.5*cvinv
