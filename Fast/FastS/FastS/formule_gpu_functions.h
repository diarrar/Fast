C    GPU-compatible macros for statement functions
C    These replace the statement functions inddm, indmtr, indven for GPU compatibility
C    Using C preprocessor macros to avoid AMD Fortran function declaration issues
C    
C    Macro for inddm - adresse point courant pour tableau de la taille d'un domaine
#define inddm_gpu(i,j,k) (1 + ((i) + param_int(NIJK+3) - 1) + ((j) + param_int(NIJK+3) - 1)*param_int(NIJK) + ((k) + param_int(NIJK+4) - 1)*param_int(NIJK)*param_int(NIJK+1))

C    Additional inddm macros for i+1, j+1, k+1 cases
#define inddm_gpu_i_p1(i,j,k) (1 + (((i)+1) + param_int(NIJK+3) - 1) + ((j) + param_int(NIJK+3) - 1)*param_int(NIJK) + ((k) + param_int(NIJK+4) - 1)*param_int(NIJK)*param_int(NIJK+1))
#define inddm_gpu_j_p1(i,j,k) (1 + ((i) + param_int(NIJK+3) - 1) + (((j)+1) + param_int(NIJK+3) - 1)*param_int(NIJK) + ((k) + param_int(NIJK+4) - 1)*param_int(NIJK)*param_int(NIJK+1))
#define inddm_gpu_k_p1(i,j,k) (1 + ((i) + param_int(NIJK+3) - 1) + ((j) + param_int(NIJK+3) - 1)*param_int(NIJK) + (((k)+1) + param_int(NIJK+4) - 1)*param_int(NIJK)*param_int(NIJK+1))

C    Macro for indmtr - adresse interface pour tableau metric
#define indmtr_gpu(i,j,k) (1 + ((i) + param_int(NIJK_MTR+3) - 1)*param_int(NIJK_MTR) + ((j) + param_int(NIJK_MTR+3) - 1)*param_int(NIJK_MTR+1) + ((k) + param_int(NIJK_MTR+4) - 1)*param_int(NIJK_MTR+2))

C    Additional indmtr macros for i+1, j+1, k+1 cases
#define indmtr_gpu_i_p1(i,j,k) (1 + (((i)+1) + param_int(NIJK_MTR+3) - 1)*param_int(NIJK_MTR) + ((j) + param_int(NIJK_MTR+3) - 1)*param_int(NIJK_MTR+1) + ((k) + param_int(NIJK_MTR+4) - 1)*param_int(NIJK_MTR+2))
#define indmtr_gpu_j_p1(i,j,k) (1 + ((i) + param_int(NIJK_MTR+3) - 1)*param_int(NIJK_MTR) + (((j)+1) + param_int(NIJK_MTR+3) - 1)*param_int(NIJK_MTR+1) + ((k) + param_int(NIJK_MTR+4) - 1)*param_int(NIJK_MTR+2))
#define indmtr_gpu_k_p1(i,j,k) (1 + ((i) + param_int(NIJK_MTR+3) - 1)*param_int(NIJK_MTR) + ((j) + param_int(NIJK_MTR+3) - 1)*param_int(NIJK_MTR+1) + (((k)+1) + param_int(NIJK_MTR+4) - 1)*param_int(NIJK_MTR+2))

C    Macro for indven - adresse interface pour tableau vitesse entrainement  
#define indven_gpu(i,j,k) (1 + ((i) + param_int(NIJK_VENT+3) - 1)*param_int(NIJK_VENT) + ((j) + param_int(NIJK_VENT+3) - 1)*param_int(NIJK_VENT+1) + ((k) + param_int(NIJK_VENT+4) - 1)*param_int(NIJK_VENT+2))

C    Additional indven macros for i+1, j+1, k+1 cases
#define indven_gpu_i_p1(i,j,k) (1 + (((i)+1) + param_int(NIJK_VENT+3) - 1)*param_int(NIJK_VENT) + ((j) + param_int(NIJK_VENT+3) - 1)*param_int(NIJK_VENT+1) + ((k) + param_int(NIJK_VENT+4) - 1)*param_int(NIJK_VENT+2))
#define indven_gpu_j_p1(i,j,k) (1 + ((i) + param_int(NIJK_VENT+3) - 1)*param_int(NIJK_VENT) + (((j)+1) + param_int(NIJK_VENT+3) - 1)*param_int(NIJK_VENT+1) + ((k) + param_int(NIJK_VENT+4) - 1)*param_int(NIJK_VENT+2))
#define indven_gpu_k_p1(i,j,k) (1 + ((i) + param_int(NIJK_VENT+3) - 1)*param_int(NIJK_VENT) + ((j) + param_int(NIJK_VENT+3) - 1)*param_int(NIJK_VENT+1) + (((k)+1) + param_int(NIJK_VENT+4) - 1)*param_int(NIJK_VENT+2))