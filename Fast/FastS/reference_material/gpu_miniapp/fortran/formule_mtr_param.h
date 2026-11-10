C     adresse point courant pour tableau metrique
C     (copie de FastC/formule_mtr_param.h)
      INTEGER_E indmtr, i_3,j_3,k_3
      indmtr(i_3,j_3,k_3) =  1
     &            + (i_3+param_int(NIJK_MTR+3)-1)*param_int(NIJK_MTR)
     &            + (j_3+param_int(NIJK_MTR+3)-1)*param_int(NIJK_MTR+1)
     &            + (k_3+param_int(NIJK_MTR+4)-1)*param_int(NIJK_MTR+2)
