# Mini-app FastS GPU, phase 0 : version Fortran (forme fixe + cpp)
# Comme FastS (Cassiopee Dist): reels par defaut en double precision (-r8)
#   make COMPILER=gnu                  : gfortran, regions target executees sur le host
#   make COMPILER=nvhpc GPU_ARCH=cc90  : nvfortran, offload NVIDIA
#   make COMPILER=amd   GPU_ARCH=gfx942: amdflang (ROCm), offload AMD
#   make COMPILER=cray                 : Cray ftn (module craype-accel-* charge)
#   make COMPILER=gnu-offload GPU_ARCH=sm_90 : gfortran avec offload nvptx (controle)
COMPILER ?= gnu
GPU_ARCH ?= cc90

SRCS = init.for flux_cpu.for flux_gpu.for lusgs_cpu.for lusgs_gpu.for main.for
EXE  = miniapp_fortran

ifeq ($(COMPILER),gnu)
  FC     = gfortran
  FFLAGS = -O3 -cpp -ffixed-form -ffixed-line-length-132 -fopenmp -march=native -fdefault-real-8 -fdefault-double-8
  LDFLAGS= -fopenmp
endif
ifeq ($(COMPILER),gnu-offload)
  FC     = gfortran
  FFLAGS = -O3 -cpp -ffixed-form -ffixed-line-length-132 -fdefault-real-8 -fdefault-double-8 -fopenmp -foffload=nvptx-none -foffload-options=-misa=$(GPU_ARCH)
  LDFLAGS= -fopenmp -foffload=nvptx-none
endif
ifeq ($(COMPILER),nvhpc)
  FC     = nvfortran
  FFLAGS = -O3 -Mpreprocess -Mfixed -Mextend -r8 -mp=gpu -gpu=$(GPU_ARCH),lineinfo -Minfo=mp
  LDFLAGS= -mp=gpu -gpu=$(GPU_ARCH)
endif
ifeq ($(COMPILER),amd)
  FC     = amdflang
  FFLAGS = -O3 -cpp -ffixed-form -ffixed-line-length=132 -fdefault-real-8 -fopenmp --offload-arch=$(GPU_ARCH)
  LDFLAGS= -fopenmp --offload-arch=$(GPU_ARCH)
endif
ifeq ($(COMPILER),cray)
  FC     = ftn
  FFLAGS = -O3 -eZ -f fixed -N 132 -s real64 -homp
  LDFLAGS= -homp
endif

OBJS = $(SRCS:.for=.o)

all: $(EXE)

$(EXE): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)

%.o: %.for mini_param.h formule_param.h formule_mtr_param.h flux_face_roe_o3.h flux_locvar.h flux_const.h lu_mat.h lu_locvar.h lu_const.h mjr_newton.h
	$(FC) $(FFLAGS) -c $< -o $@

run: $(EXE)
	./$(EXE) 64 64 64 10

clean:
	rm -f $(OBJS) $(EXE) *.mod
