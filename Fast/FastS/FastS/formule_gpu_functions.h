C    GPU-compatible macros for statement functions with AMD GPU memory access fix
C    These replace the statement functions inddm, indmtr, indven for GPU compatibility
C    Using pre-computed parameter values to avoid array access within GPU kernels
C    This eliminates AMD GPU memory access faults by removing parameter array dependencies
C    
C    Macro for inddm - adresse point courant pour tableau de la taille d'un domaine
#define inddm_gpu(i,j,k) 1 + (i+nijk3_val-1) + (j+nijk3_val-1)*nijk_val + (k+nijk4_val-1)*nijk_val*nijk1_val

C    Additional inddm macros for i+1, j+1, k+1 cases
#define inddm_gpu_i_p1(i,j,k) 1 + ((i+1)+nijk3_val-1) + (j+nijk3_val-1)*nijk_val + (k+nijk4_val-1)*nijk_val*nijk1_val
#define inddm_gpu_j_p1(i,j,k) 1 + (i+nijk3_val-1) + ((j+1)+nijk3_val-1)*nijk_val + (k+nijk4_val-1)*nijk_val*nijk1_val
#define inddm_gpu_k_p1(i,j,k) 1 + (i+nijk3_val-1) + (j+nijk3_val-1)*nijk_val + ((k+1)+nijk4_val-1)*nijk_val*nijk1_val

C    Macro for indmtr - adresse interface pour tableau metric
#define indmtr_gpu(i,j,k) 1 + (i+nijk_mtr3_val-1)*nijk_mtr_val + (j+nijk_mtr3_val-1)*nijk_mtr1_val + (k+nijk_mtr4_val-1)*nijk_mtr2_val

C    Additional indmtr macros for i+1, j+1, k+1 cases
#define indmtr_gpu_i_p1(i,j,k) 1 + ((i+1)+nijk_mtr3_val-1)*nijk_mtr_val + (j+nijk_mtr3_val-1)*nijk_mtr1_val + (k+nijk_mtr4_val-1)*nijk_mtr2_val
#define indmtr_gpu_j_p1(i,j,k) 1 + (i+nijk_mtr3_val-1)*nijk_mtr_val + ((j+1)+nijk_mtr3_val-1)*nijk_mtr1_val + (k+nijk_mtr4_val-1)*nijk_mtr2_val
#define indmtr_gpu_k_p1(i,j,k) 1 + (i+nijk_mtr3_val-1)*nijk_mtr_val + (j+nijk_mtr3_val-1)*nijk_mtr1_val + ((k+1)+nijk_mtr4_val-1)*nijk_mtr2_val

C    Macro for indven - adresse interface pour tableau vitesse entrainement  
#define indven_gpu(i,j,k) 1 + (i+nijk_vent3_val-1)*nijk_vent_val + (j+nijk_vent3_val-1)*nijk_vent1_val + (k+nijk_vent4_val-1)*nijk_vent2_val

C    Additional indven macros for i+1, j+1, k+1 cases
#define indven_gpu_i_p1(i,j,k) 1 + ((i+1)+nijk_vent3_val-1)*nijk_vent_val + (j+nijk_vent3_val-1)*nijk_vent1_val + (k+nijk_vent4_val-1)*nijk_vent2_val
#define indven_gpu_j_p1(i,j,k) 1 + (i+nijk_vent3_val-1)*nijk_vent_val + ((j+1)+nijk_vent3_val-1)*nijk_vent1_val + (k+nijk_vent4_val-1)*nijk_vent2_val
#define indven_gpu_k_p1(i,j,k) 1 + (i+nijk_vent3_val-1)*nijk_vent_val + (j+nijk_vent3_val-1)*nijk_vent1_val + ((k+1)+nijk_vent4_val-1)*nijk_vent2_val