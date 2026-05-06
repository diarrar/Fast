#ifdef _GPU_OFFLOAD
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
!$OMP& DEFAULT(PRESENT) FIRSTPRIVATE(ind_loop)
#endif
      do k = ind_loop(5), ind_loop(6)
      do j = ind_loop(3), ind_loop(4)
      do i = ind_loop(1), ind_loop(2)

        l   = inddm(i,j,k)
