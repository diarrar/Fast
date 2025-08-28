C    GPU-compatible device functions for statement functions
C    These replace the statement functions inddm, indmtr, indven for GPU compatibility
C    
C    Device function for inddm - adresse point courant pour tableau de la taille d'un domaine
!$OMP DECLARE TARGET(inddm_gpu)
      INTEGER_E FUNCTION inddm_gpu(i_1, j_1, k_1, param_int)
      implicit none
      INTEGER_E i_1, j_1, k_1
      INTEGER_E param_int(0:136)
      
      inddm_gpu = 1 
     &     + (i_1 + param_int(3) - 1)
     &     + (j_1 + param_int(3) - 1) * param_int(0)
     &     + (k_1 + param_int(4) - 1) * param_int(0) * param_int(1)
     
      END FUNCTION inddm_gpu

C    Device function for indmtr - adresse interface pour tableau metric  
!$OMP DECLARE TARGET(indmtr_gpu)
      INTEGER_E FUNCTION indmtr_gpu(i_3, j_3, k_3, param_int)
      implicit none
      INTEGER_E i_3, j_3, k_3  
      INTEGER_E param_int(0:136)
      
      indmtr_gpu = 1
     &            + (i_3 + param_int(8) - 1) * param_int(5)
     &            + (j_3 + param_int(8) - 1) * param_int(6)  
     &            + (k_3 + param_int(9) - 1) * param_int(7)
     
      END FUNCTION indmtr_gpu

C    Device function for indven - adresse interface pour tableau vitesse entrainement
!$OMP DECLARE TARGET(indven_gpu)
      INTEGER_E FUNCTION indven_gpu(i_4, j_4, k_4, param_int)
      implicit none
      INTEGER_E i_4, j_4, k_4
      INTEGER_E param_int(0:136)
      
      indven_gpu = 1
     &           + (i_4 + param_int(18) - 1) * param_int(15)
     &           + (j_4 + param_int(18) - 1) * param_int(16)
     &           + (k_4 + param_int(19) - 1) * param_int(17)
     
      END FUNCTION indven_gpu