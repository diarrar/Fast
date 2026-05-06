#
# -- Compilation des fortrans --
#
# Definis dans le Makefile :
# SOURCEF77 : liste des source f77
# LIBNAME : nom librairie
# Repertoire de production des .o fortrans
PRODDIR   = ../../build/$(ELSAPROD)/

# -- Ne pas modifier ci-dessous --
SRC77     = $(shell ../MakeObjs.mk "$(SOURCEF77)")
OBJS      = $(SRC77:%.for=%.o)
OBJS2     = $(shell ../MakeObjsPath.mk "$(OBJS)") 
LIB2      = $(shell ../MakeObjsPath.mk "lib$(LIBNAME).a")
all: link $(LIB2)

.SUFFIXES: .for .f

# GPU offload support via OpenMP target offload (OpenMP 4.5/5.x)
# Usage: make GPU_VENDOR=nvidia  (NVIDIA HPC SDK, nvfortran)
#        make GPU_VENDOR=amd     (ROCm, amdflang)
# When GPU_VENDOR is set:
#   - E_SCALAR_COMPUTER and E_OMP_SOUS_DOMAIN are omitted  -> 3D i,j,k loop path
#   - _GPU_OFFLOAD is defined                              -> !$OMP TARGET directives active
#   - _GPU_NVIDIA or _GPU_AMD is defined for vendor tuning
ifeq ($(GPU_VENDOR),nvidia)
  CPP_GPU_DEFS = -D_GPU_OFFLOAD -D_GPU_NVIDIA
  GPU_FFLAGS   = -mp=gpu -gpu=cc80,cc90 -Minfo=accel
else ifeq ($(GPU_VENDOR),amd)
  CPP_GPU_DEFS = -D_GPU_OFFLOAD -D_GPU_AMD
  GPU_FFLAGS   = -fopenmp-targets=amdgcn-amd-amdhsa \
                 -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a
else
  CPP_GPU_DEFS =
  GPU_FFLAGS   =
endif

# CPU-only flags: dropped in GPU builds
ifeq ($(GPU_VENDOR),)
  CPP_CPU_ONLY = -DE_SCALAR_COMPUTER -DE_OMP_SOUS_DOMAIN
else
  CPP_CPU_ONLY =
endif

.for.o:
	echo "Compiling fortran 77: " $<
	@(cpp -traditional -I../../../KCore/KCore \
	  -DINTEGER_E="INTEGER*4" -DREAL_E="REAL*8" \
	  $(CPP_CPU_ONLY) $(CPP_GPU_DEFS) \
	  $< > $(PRODDIR)/temp.f; \
	$(FC) $(FOPT0) $(FOPT1) $(FOPT2) $(GPU_FFLAGS) \
	  -c $(PRODDIR)/temp.f -o $@; rm $(PRODDIR)/temp.f)

link:	$(SRC77)
	@if [ -d ../../build/ ]; then :; else mkdir ../../build; fi
	@if [ -d $(PRODDIR) ]; then :; else mkdir $(PRODDIR); fi
	@for i in $(SRC77); do \
	../MakeLink.mk "`pwd`" $$i $(PRODDIR);\
	done

$(LIB2) : $(OBJS2)
	@echo Making $(DIRECTORY) library
	@ar clr  $@ $(OBJS2) ; ranlib $(LIB2)
	@echo done
