/*
    Copyright 2013-2025 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * gpu_data_map.h  —  Phase 1: OpenMP target data residency macros
 *
 * PURPOSE
 * -------
 * Provides C++ macros that establish a persistent OpenMP target data region
 * (device-side resident arrays) around the time-stepping loop in rhs.cpp /
 * lhs.cpp.  When _GPU_OFFLOAD is not defined every macro expands to nothing,
 * keeping the CPU build exactly unchanged.
 *
 * USAGE (in rhs.cpp / lhs.cpp, around the time-step loop)
 * --------------------------------------------------------
 *   E_Int ndimdx = param_int[nd][NDIMDX];
 *   E_Int neq    = param_int[nd][NEQ];
 *   E_Int ndimdx_mtr = param_int[nd][NDIMDX_MTR];
 *   E_Int neq_ij     = param_int[nd][NEQ_IJ];
 *   E_Int neq_k      = param_int[nd][NEQ_K];
 *
 *   FAST_GPU_ENTER_DATA(nd, ndimdx, neq, ndimdx_mtr, neq_ij, neq_k,
 *                       iptro[nd], iptrodm[nd],
 *                       ipti[nd], iptj[nd], iptk[nd], iptvol[nd],
 *                       iptmut[nd], param_int[nd], param_real[nd])
 *   {
 *       // ... time loop ...
 *       // For MPI halo exchange / I/O, bring rop back to host:
 *       FAST_GPU_UPDATE_FROM(iptro[nd], ndimdx*neq)
 *       // ... MPI calls ...
 *       FAST_GPU_UPDATE_TO(iptro[nd], ndimdx*neq)
 *   }
 *   FAST_GPU_EXIT_DATA(iptro[nd], ndimdx*neq)
 *
 * DATA MAPPING POLICY
 * -------------------
 *   Metrics (ti/tj/tk/vol)     : map(to:...)    — read-only, fixed for steady grids
 *   Solution (rop)             : map(tofrom:...) — in at start, out for MPI/I/O
 *   RHS increment (drodm)      : map(alloc:...)  — scratch only, never leaves GPU
 *   Viscosity (xmut)           : map(tofrom:...) — updated every RHS call
 *   param_int / param_real     : map(to:...)     — read-only metadata
 *
 * NOTES
 * -----
 *  - DEFAULT(PRESENT) in the Fortran loop macros (loopGpu_begin.for) requires
 *    that ALL arrays touched inside the loop are already mapped here before
 *    the first kernel call.
 *  - ALE / moving-grid cases: re-map ti/tj/tk/vol as tofrom if the grid moves.
 *  - Multi-GPU / MPI: call FAST_GPU_UPDATE_FROM / FAST_GPU_UPDATE_TO around
 *    every MPI_Isend / MPI_Irecv halo exchange.
 */

#ifndef _FASTS_GPU_DATA_MAP_H_
#define _FASTS_GPU_DATA_MAP_H_

#ifdef _GPU_OFFLOAD

/* Enter persistent device data region for one zone (nd).
 * Must be paired with FAST_GPU_EXIT_DATA. */
#define FAST_GPU_ENTER_DATA(nd,                                            \
        ndimdx_, neq_, ndimdx_mtr_, neq_ij_, neq_k_,                     \
        rop_,  drodm_,                                                     \
        ti_,   tj_,   tk_,  vol_,                                         \
        xmut_, param_int_, param_real_)                                   \
    _Pragma("omp target enter data map(to:                               " \
            "  param_int_[0:200],                                        " \
            "  param_real_[0:100],                                       " \
            "  ti_  [0:(ndimdx_mtr_)*(neq_ij_)],                        " \
            "  tj_  [0:(ndimdx_mtr_)*(neq_ij_)],                        " \
            "  tk_  [0:(ndimdx_mtr_)*(neq_k_) ],                        " \
            "  vol_ [0:(ndimdx_mtr_)           ])")                       \
    _Pragma("omp target enter data map(tofrom:                           " \
            "  rop_  [0:(ndimdx_)*(neq_)],                               " \
            "  xmut_ [0:(ndimdx_)       ])")                              \
    _Pragma("omp target enter data map(alloc:                            " \
            "  drodm_[0:(ndimdx_)*(neq_)])")

/* Exit persistent device data region and copy solution back to host. */
#define FAST_GPU_EXIT_DATA(rop_, ndimdx_, neq_, xmut_)                   \
    _Pragma("omp target exit data map(from:                              " \
            "  rop_  [0:(ndimdx_)*(neq_)],                               " \
            "  xmut_ [0:(ndimdx_)       ])")

/* Bring rop from device to host (before MPI exchange or file output). */
#define FAST_GPU_UPDATE_FROM(ptr_, size_)                                 \
    _Pragma("omp target update from(ptr_[0:(size_)])")

/* Send updated rop from host to device (after MPI halo fill). */
#define FAST_GPU_UPDATE_TO(ptr_, size_)                                   \
    _Pragma("omp target update to(ptr_[0:(size_)])")

/* Re-map metrics when the grid has moved (ALE). */
#define FAST_GPU_UPDATE_METRICS(ti_, tj_, tk_, vol_,                      \
                                ndimdx_mtr_, neq_ij_, neq_k_)            \
    _Pragma("omp target update to(                                       " \
            "  ti_  [0:(ndimdx_mtr_)*(neq_ij_)],                        " \
            "  tj_  [0:(ndimdx_mtr_)*(neq_ij_)],                        " \
            "  tk_  [0:(ndimdx_mtr_)*(neq_k_) ],                        " \
            "  vol_ [0:(ndimdx_mtr_)           ])")

#else /* CPU-only build: macros expand to nothing */

#define FAST_GPU_ENTER_DATA(nd,                                           \
        ndimdx_, neq_, ndimdx_mtr_, neq_ij_, neq_k_,                     \
        rop_, drodm_, ti_, tj_, tk_, vol_, xmut_,                         \
        param_int_, param_real_)
#define FAST_GPU_EXIT_DATA(rop_, ndimdx_, neq_, xmut_)
#define FAST_GPU_UPDATE_FROM(ptr_, size_)
#define FAST_GPU_UPDATE_TO(ptr_, size_)
#define FAST_GPU_UPDATE_METRICS(ti_, tj_, tk_, vol_,                      \
                                ndimdx_mtr_, neq_ij_, neq_k_)

#endif /* _GPU_OFFLOAD */

#endif /* _FASTS_GPU_DATA_MAP_H_ */
