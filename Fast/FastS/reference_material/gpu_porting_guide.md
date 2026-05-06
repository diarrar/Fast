# GPU Porting Guide — Fast/FastS (OpenMP Target Offload)

> **Technology choice:** OpenMP 4.5/5.x target offload  
> **Supported platforms:** NVIDIA (NVHPC SDK `nvfortran`) · AMD (ROCm `amdflang`)  
> **OpenACC:** excluded (insufficient compiler support on NVIDIA HPC SDK ≥ 24.x)

---

## 1. Architecture overview

Fast/FastS uses ~2 200 Fortran 77 files (`.for`, preprocessed by CPP) for all
numerical kernels.  The loop structure of every kernel is abstracted behind
CPP `#include` macros so that a single change to those files propagates to the
entire codebase:

```
loop_begin.for          ──►  loopGpu_begin.for   (GPU path, E_SCALAR_COMPUTER unset)
   ... kernel body ...
loop_end.for            ──►  loopGpu_end.for     (GPU path)
```

The GPU path is activated by **not** defining `E_SCALAR_COMPUTER`
(which is the default in SCons GPU builds, see §3).

When `_GPU_OFFLOAD` is also defined, `loopGpu_begin.for` emits:

```fortran
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
!$OMP& DEFAULT(PRESENT) FIRSTPRIVATE(ind_loop)
do k = ind_loop(5), ind_loop(6)
do j = ind_loop(3), ind_loop(4)
do i = ind_loop(1), ind_loop(2)
   l   = inddm(i,j,k)
   lt  = indmtr(i,j,k)
   ...
```

This single macro covers **all explicit RHS kernels** (Roe flux, SA model,
DTLOC, viscosity, etc.) — roughly 1 500 Fortran files.

---

## 2. CPP flags

| Flag | GPU build | CPU build | Purpose |
|---|---|---|---|
| `E_SCALAR_COMPUTER` | **absent** | defined | CPU vectorised flat-loop path |
| `E_OMP_SOUS_DOMAIN` | **absent** | defined | CPU OpenMP sub-domain splitting |
| `_GPU_OFFLOAD` | defined | absent | activates `!$OMP TARGET` directives |
| `_GPU_NVIDIA` | set when `GPU_VENDOR=nvidia` | absent | vendor-specific tuning |
| `_GPU_AMD` | set when `GPU_VENDOR=amd` | absent | vendor-specific tuning |

---

## 3. Build instructions

### 3.1 SCons (recommended)

```bash
# NVIDIA A100/H100 (NVIDIA HPC SDK ≥ 24.1, nvfortran)
export GPU_VENDOR=nvidia
export ELSAPROD=gpu_nvidia_r8
scons -j8

# AMD MI250X/MI300X (ROCm ≥ 6.0, amdflang)
export GPU_VENDOR=amd
export ELSAPROD=gpu_amd_r8
scons -j8

# CPU-only (unchanged, gfortran/ifort)
unset GPU_VENDOR
export ELSAPROD=x86_r8
scons -j8
```

### 3.2 Makefile (legacy)

```bash
# NVIDIA
make GPU_VENDOR=nvidia FC=nvfortran

# AMD
make GPU_VENDOR=amd FC=amdflang

# CPU (unchanged)
make
```

### 3.3 Hardware target overrides

Edit `setup.scons` (or `Make.inc`) to change the GPU architecture:

```python
# NVIDIA — change cc80 (A100) to cc90 (H100) or cc70 (V100)
addFFlags += ['-mp=gpu', '-gpu=cc80,cc90', '-Minfo=accel']

# AMD — change gfx90a (MI250X) to gfx942 (MI300X)
addFFlags += ['-fopenmp-targets=amdgcn-amd-amdhsa',
              '-Xopenmp-target=amdgcn-amd-amdhsa', '-march=gfx942']
```

---

## 4. Porting phases

### Phase 0 — Infrastructure ✅ (implemented in this PR)

- [x] `loopGpu_begin.for` / `loop3dcartGpu_begin.for` — `!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)` directive
- [x] `loopGpu_end.for` — matching `!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD`
- [x] SCons GPU build detection (`GPU_VENDOR` environment variable)
- [x] Makefile GPU build support (`GPU_VENDOR` make variable)
- [x] `gpu_data_map.h` — C++ macros for Phase 1 data residency

### Phase 1 — Data residency (pending)

**Goal:** keep all zone arrays permanently on the GPU between time steps.
Transfers happen only at MPI halo exchanges and file output.

**How to use `gpu_data_map.h`** (in `FastS/FastS/Compute/rhs.cpp` and `lhs.cpp`):

```cpp
#include "FastS/gpu_data_map.h"

// Before the time-step loop for zone nd:
E_Int ndimdx     = param_int[nd][NDIMDX];
E_Int neq        = param_int[nd][NEQ];
E_Int ndimdx_mtr = param_int[nd][NDIMDX_MTR];
E_Int neq_ij     = param_int[nd][NEQ_IJ];
E_Int neq_k      = param_int[nd][NEQ_K];

FAST_GPU_ENTER_DATA(nd,
    ndimdx, neq, ndimdx_mtr, neq_ij, neq_k,
    iptro[nd], iptrodm[nd],
    ipti[nd], iptj[nd], iptk[nd], iptvol[nd],
    iptmut[nd], param_int[nd], param_real[nd])

for (E_Int it = 0; it < niter; it++)
{
    // ... navier_stokes_struct_() call (kernels run on GPU) ...

    // Before MPI halo exchange:
    FAST_GPU_UPDATE_FROM(iptro[nd], ndimdx * neq)
    // ... MPI_Isend / MPI_Irecv ...
    // After halo fill:
    FAST_GPU_UPDATE_TO(iptro[nd], ndimdx * neq)
}

FAST_GPU_EXIT_DATA(iptro[nd], ndimdx, neq, iptmut[nd])
```

**Data mapping policy:**

| Array | Direction | Rationale |
|---|---|---|
| `ti`, `tj`, `tk`, `vol` | `map(to:...)` | Read-only metrics; re-map only if ALE |
| `param_int`, `param_real` | `map(to:...)` | Read-only zone metadata |
| `rop` (solution) | `map(tofrom:...)` | Enters with IC; leaves for MPI/output |
| `drodm` (RHS scratch) | `map(alloc:...)` | Never needs to reach the CPU |
| `xmut` (turbulent viscosity) | `map(tofrom:...)` | Updated every RHS iteration |

### Phase 2 — Explicit RHS offload (pending)

Once Phase 1 data residency is in place, the `!$OMP TARGET` directives added
to `loopGpu_begin.for` automatically offload **all** explicit kernels:

- Roe flux: `fluroe_[euler/lamin/SA]_[o3/minmod]_3dfull.for`
- Laminar viscosity: `invist.for`
- Spalart-Allmaras: `Compute/SA/*.for`
- Local time step: `Compute/DTLOC/*.for`
- Boundary conditions: `BC/bvbs_*.for` (the `_GPU_gpu.for` variants already exist)

**Benchmarking target:** ×10–30 speedup for Roe flux vs. 16 CPU cores.

### Phase 3 — Implicit LU-SSOR offload (pending)

The LU-SSOR sweeps (`lu_[i/j/k]_3dfull.for`, `invlu_[l/u].for`) have
wavefront dependencies that prevent trivial GPU offload.

**Recommended approach (option B):**  
Keep LU-SSOR on CPU initially.  Offload only the RHS.  Profile whether the
remaining CPU cost is significant before investing in wavefront parallelism.

**Option A (wavefront parallelism):**  
For sweep in direction `i`, planes `(j,k)` are independent.  Use:
```fortran
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
!$OMP& DEFAULT(PRESENT)
do k = ind_loop(5), ind_loop(6)
do j = ind_loop(3), ind_loop(4)
  ! sequential sweep in i inside
  do i = ind_loop(1), ind_loop(2)
    ...
```

### Phase 4 — Multi-GPU MPI (pending)

- Map 1 MPI rank per GPU (current 1 rank/socket → 1 rank/GPU)
- Enable CUDA-aware or ROCm-aware MPI for peer-to-peer transfers
- Use `!$OMP TARGET NOWAIT` + `!$OMP TASKWAIT` to overlap halo exchange
  with interior-domain computation

---

## 5. Debugging and verification

### Detect missing data maps

With `DEFAULT(PRESENT)`, any array not yet in the device data environment
causes a **runtime error** rather than silent wrong results.  Use this to
discover missing maps systematically.

### Unified memory fallback (NVIDIA only)

For initial kernel validation before Phase 1 is complete, add `-gpu=managed`
to the nvfortran flags.  This enables CUDA Unified Memory so all allocations
are automatically accessible on the GPU.  **Not for production** (slower due
to page-fault overhead).

```bash
GPU_VENDOR=nvidia scons  # edit setup.scons: '-gpu=cc80,managed'
```

### Numerical regression test

A reference L2 norm on a small test case should be compared between CPU and
GPU results.  Acceptable tolerance: ε ≈ 1e-12 for double precision (round-off
differences are expected due to SIMD/FMA reordering).

---

## 6. Key files changed

| File | Change |
|---|---|
| `FastS/FastS/Compute/loopGpu_begin.for` | Added `!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)` |
| `FastS/FastS/Compute/loopGpu_end.for` | Added `!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD` |
| `FastS/FastS/Compute/loop3dcartGpu_begin.for` | Same (Cartesian 3D variant) |
| `FastC/FastC/HPC_LAYER/loopGpu_begin.for` | Same (FastC HPC layer) |
| `FastC/FastC/HPC_LAYER/loopGpu_end.for` | Same |
| `FastS/setup.scons` | Added `GPU_VENDOR` detection, nvfortran/amdflang flags |
| `FastC/setup.scons` | Same |
| `FastS/FastS/Make.inc` | Added `GPU_VENDOR` detection for Makefile builds |
| `FastC/FastC/Make.inc` | Same |
| `FastS/FastS/gpu_data_map.h` | New — Phase 1 data residency C++ macros |

---

## 7. Compiler compatibility notes

| Compiler | OpenMP target support | Min version | Notes |
|---|---|---|---|
| `nvfortran` (NVIDIA HPC SDK) | Full OpenMP 5.0 target | 22.x | Use `-mp=gpu -gpu=ccXX` |
| `amdflang` (ROCm) | OpenMP 4.5 stable, 5.0 partial | ROCm 6.0 | Use `-fopenmp-targets=amdgcn-amd-amdhsa` |
| Cray Fortran (`ftn`) | OpenMP 5.1 target | PE 23.x | Supports both NVIDIA and AMD backends |
| `gfortran` | No GPU target | — | CPU fallback only |
| `ifort` / `ifx` | Limited GPU target | oneAPI 2024.x | Not recommended for production |
