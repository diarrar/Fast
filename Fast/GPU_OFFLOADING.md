# Fast GPU Offloading with OpenMP

This document describes the GPU offloading support that has been added to Fast using OpenMP target directives.

## Overview

Fast now includes OpenMP GPU offloading support for the main computational loops and some boundary conditions. The implementation uses conditional compilation to maintain compatibility with existing builds.

## Build Configuration

### Compiler Support
- **Intel Fortran (ifort)**: Uses `-qopenmp -qopenmp-offload=generic`
- **GNU Fortran (gfortran)**: Uses `-fopenmp -foffload=<target>`

### Environment Variables

#### `OMP_OFFLOAD_TARGET`
Controls the offloading target for GCC builds:
- `disable`: Disable GPU offloading (CPU fallback) - **Default for compatibility**
- `nvptx-none`: NVIDIA GPU offloading
- `amdgcn-amdhsa`: AMD GPU offloading

#### Usage Examples:
```bash
# Build with GPU offloading disabled (default)
./install

# Build with NVIDIA GPU offloading
export OMP_OFFLOAD_TARGET=nvptx-none
./install

# Build with AMD GPU offloading  
export OMP_OFFLOAD_TARGET=amdgcn-amdhsa
./install
```

## Code Changes

### Affected Files
- **Loop templates**: `loopGpu_begin.for`, `loopGpu_end.for`, `loop3dcartGpu_begin.for`
- **Main computation**: `template_FluxAndBalance_gpu.for`
- **Boundary conditions**: `bvbs_wall_viscous_GPU_gpu.for`

### Template Architecture

Fast now maintains two master templates for flux computations:

- **`template_FluxAndBalance.for`**: Base template with conditional GPU directives
  - Contains `#ifdef _OPENMP_GPU_OFFLOAD` blocks for unified codebase
  - Used when building CPU-only versions
  
- **`template_FluxAndBalance_gpu.for`**: GPU master template
  - Optimized structure for GPU execution
  - Used when building GPU-enabled versions

### Code Generation

The flux generation process uses the appropriate template based on the command:

```bash
cd Fast/FastS/FastS/Compute

# Generate CPU flux routines (uses template_FluxAndBalance.for)
python generate_flu.py ROE

# Generate GPU flux routines (uses template_FluxAndBalance_gpu.for)  
python generate_flu.py ROE GPU
```

Where `<FLUX_TYPE>` can be:
- `ROE` - Roe flux scheme 
- `AUSM` - AUSM flux scheme
- `SENSOR` - Sensor flux scheme
- `SENSOR_INIT` - Sensor initialization flux scheme

Example to regenerate GPU-enabled ROE flux routines:
```bash
python generate_flu.py ROE GPU
```

This generates GPU-accelerated versions of all flux routines in the corresponding subdirectories (`ROE/3dfull/`, `ROE/3dcart/`, etc.).

### Implementation Details

#### Target Data Regions
GPU data mapping is handled with `TARGET DATA` directives:
```fortran
#ifdef _OPENMP_GPU_OFFLOAD
!$OMP TARGET DATA MAP(to: input_arrays) MAP(tofrom: output_arrays)
#endif
```

#### Parallel Execution
Main computational loops use:
```fortran
#ifdef _OPENMP_GPU_OFFLOAD
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#endif
```

#### Inner Loop Vectorization
Innermost loops maintain SIMD directives for optimization:
```fortran
!$OMP SIMD
```

### Mapped Arrays
Key arrays transferred to GPU:
- `rop`, `drodm` (flow variables)
- `wig` (grid metrics)
- `venti`, `ventj`, `ventk` (velocities)
- `ti`, `tj`, `tk` (metrics)
- `vol` (volumes)
- `xmut` (viscosity)

## Runtime Usage

### Environment Variables
- `OMP_TARGET_OFFLOAD=MANDATORY`: Force GPU execution
- `OMP_TARGET_OFFLOAD=DISABLED`: Force CPU execution
- `OMP_NUM_TEAMS`: Control number of GPU teams

### Example Runtime:
```bash
# Force GPU execution
export OMP_TARGET_OFFLOAD=MANDATORY
export OMP_NUM_TEAMS=1024
./your_fast_application

# Profile GPU usage (NVIDIA)
nvprof ./your_fast_application
```

## Testing GPU Offloading

### Recommended Test Cases

For verifying GPU functionality, use these Fast test cases:

**Simple 2D Tests:**
- `lambPT.py` - Lamb vortex (good for initial validation)
- `cylindrePT.py` - 2D cylinder flow

**3D Tests:**  
- `flatPlatePT.py` - 3D flat plate boundary layer
- `PplaneNasaSAMesh136x96_t1.py` - NASA turbulent flow test

### Test Procedure

1. **Generate GPU flux routines**:
   ```bash
   cd Fast/FastS/FastS/Compute
   python generate_flu.py ROE GPU
   ```

2. **Build with GPU support**:
   ```bash
   # For NVIDIA GPUs
   export OMP_OFFLOAD_TARGET=nvptx-none
   ./install
   
   # For AMD GPUs  
   export OMP_OFFLOAD_TARGET=amdgcn-amdhsa
   ./install
   ```

3. **Run test with GPU enforcement**:
   ```bash
   export OMP_TARGET_OFFLOAD=MANDATORY
   cd Fast/FastS/test
   python lambPT.py
   ```

### Verification

- Check OpenMP runtime output: `export LIBOMPTARGET_DEBUG=1`
- Compare results with CPU-only runs for accuracy
- Monitor GPU utilization with `nvidia-smi` or `rocm-smi`

## Compatibility

The implementation maintains full backward compatibility:
- Existing builds continue to work unchanged
- GPU features only activate when explicitly enabled
- CPU fallback automatic when GPU unavailable

## Performance Notes

- GPU offloading is most beneficial for large computational domains
- Memory transfer overhead may impact performance for small problems
- Consider using `OMP_TARGET_OFFLOAD=MANDATORY` for debugging GPU code paths

## Troubleshooting

### Common Issues
1. **Compilation fails**: Check if compiler supports OpenMP 4.5+
2. **Runtime errors**: Ensure GPU drivers and OpenMP runtime are installed
3. **No performance gain**: Profile to check if data transfers dominate

### Debug Information
Enable verbose OpenMP output:
```bash
export LIBOMPTARGET_DEBUG=1
```