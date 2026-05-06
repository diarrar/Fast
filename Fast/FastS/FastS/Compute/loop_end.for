#ifndef E_SCALAR_COMPUTER
             enddo
             enddo
             enddo
#ifdef _GPU_OFFLOAD
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#endif
#else
             enddo
             enddo
             enddo
#endif
