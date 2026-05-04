       do k = ind_loop(5), ind_loop(6)
       do j = ind_loop(3), ind_loop(4)
#ifdef _OPENMP_GPU_OFFLOAD
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#endif
       do i = ind_loop(1), ind_loop(2)

        l   = inddm(i,j,k)
        lt  = indmtr(i,j,k)
        lvo = lt
