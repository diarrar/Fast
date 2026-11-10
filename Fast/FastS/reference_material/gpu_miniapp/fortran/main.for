c***********************************************************************
c     Mini-application FastS GPU, phase 0 (version Fortran)
c     usage : miniapp_fortran [ni nj nk nrep]
c     Mesure et verifie :
c       flux_cpu   : flux Roe o3, formulation scatter FastS (reference)
c       flux_gpu   : flux Roe o3, formulation gather, OpenMP target
c       lusgs_cpu  : LU-SGS lexicographique (reference)
c       lusgs_gpu  : LU-SGS par hyperplans, OpenMP target
c***********************************************************************
      program miniapp_fortran
      implicit none
#include "mini_param.h"
      INTEGER_E param_int(0:NPARAM-1)
      REAL_E    param_real(0:NPARAM-1)
      REAL_E, allocatable :: rop(:), rop_1(:), drodm(:), drodm_ref(:),
     &  drodm_out(:), coe(:), ti(:), tj(:), tk(:), vol(:), rop_ref(:)
      integer, allocatable :: plane_ptr(:), plane_cell(:)
      integer ni,nj,nk,nrep,ndim,nx,ny,nz,it,nplane,narg
      integer*8 t0,t1,rate
      real*8 tcpu,tgpu,ncell
      character*32 arg

      ni = 64
      nj = 64
      nk = 64
      nrep = 10
      narg = command_argument_count()
      if (narg.ge.3) then
        call get_command_argument(1,arg)
        read(arg,*) ni
        call get_command_argument(2,arg)
        read(arg,*) nj
        call get_command_argument(3,arg)
        read(arg,*) nk
      endif
      if (narg.ge.4) then
        call get_command_argument(4,arg)
        read(arg,*) nrep
      endif

c     parametres entiers (cles FastS)
      nx = ni+4
      ny = nj+4
      nz = nk+4
      ndim = nx*ny*nz
      param_int(:) = 0
      param_int(NIJK)   = nx
      param_int(NIJK+1) = ny
      param_int(NIJK+2) = nz
      param_int(NIJK+3) = 2
      param_int(NIJK+4) = 2
      param_int(NIJK_MTR)   = 1
      param_int(NIJK_MTR+1) = nx
      param_int(NIJK_MTR+2) = nx*ny
      param_int(NIJK_MTR+3) = 2
      param_int(NIJK_MTR+4) = 2
      param_int(IJKV)   = ni
      param_int(IJKV+1) = nj
      param_int(IJKV+2) = nk
      param_int(NEQ)     = 5
      param_int(NEQ_IJ)  = 3
      param_int(NEQ_K)   = 3
      param_int(NEQ_COE) = 5
      param_int(NDIMDX)     = ndim
      param_int(NDIMDX_MTR) = ndim
      param_real(:) = 0.d0
      param_real(GAMMA)  = 1.4d0
      param_real(CVINF)  = 2.5d0
      param_real(PSIROE) = 0.1d0
      param_real(DTC)    = 2.d0/ni

      allocate(rop(ndim*5), rop_1(ndim*5), drodm(ndim*5))
      allocate(drodm_ref(ndim*5), drodm_out(ndim*5), coe(ndim*5))
      allocate(rop_ref(ndim*5))
      allocate(ti(ndim*3), tj(ndim*3), tk(ndim*3), vol(ndim))
      allocate(plane_ptr(ni+nj+nk), plane_cell(ni*nj*nk))

      call init_fields(param_int, param_real, rop, ti,tj,tk, vol, coe)
      call build_planes(ni,nj,nk,nx,ny, nplane, plane_ptr, plane_cell)
      ncell = dble(ni)*dble(nj)*dble(nk)
      rop_1(:) = 0.d0
      drodm(:) = 0.d0
      drodm_out(:) = 0.d0
      write(*,'(a,4(1x,a,i0))') 'INFO impl=fortran', 'ni=',ni,
     &  'nj=',nj, 'nk=',nk, 'nrep=',nrep
      call system_clock(count_rate=rate)

c     residence des donnees sur le device (equivalent warmup FastS)
!$OMP TARGET ENTER DATA MAP(TO: rop, ti, tj, tk, coe, param_int,
!$OMP&  param_real, plane_ptr, plane_cell)
!$OMP& MAP(ALLOC: drodm, drodm_ref, drodm_out, rop_1)

c----- flux CPU (reference)
      call flux_roe_cpu(param_int, param_real, rop, ti,tj,tk, drodm)
      call system_clock(t0)
      do it = 1, nrep
        call flux_roe_cpu(param_int, param_real, rop, ti,tj,tk, drodm)
      enddo
      call system_clock(t1)
      tcpu = dble(t1-t0)/dble(rate)/nrep
      write(*,'(a,1x,a,es12.4)') 'TIME flux_cpu', 'ns_per_cell=',
     &  tcpu/ncell*1.d9
      call checksum(param_int, 'flux_cpu', drodm)
      drodm_ref(:) = drodm(:)
!$OMP TARGET UPDATE TO(drodm_ref)

c----- flux GPU (gather, target)
      call flux_roe_gpu(param_int, param_real, rop, ti,tj,tk, drodm)
      call system_clock(t0)
      do it = 1, nrep
        call flux_roe_gpu(param_int, param_real, rop, ti,tj,tk, drodm)
      enddo
      call system_clock(t1)
      tgpu = dble(t1-t0)/dble(rate)/nrep
      write(*,'(a,1x,a,es12.4)') 'TIME flux_gpu', 'ns_per_cell=',
     &  tgpu/ncell*1.d9
!$OMP TARGET UPDATE FROM(drodm)
      call checksum(param_int, 'flux_gpu', drodm)
      call maxdiff(param_int, 'flux_gpu_vs_cpu', drodm_ref, drodm)

c----- LU-SGS CPU (reference)
      call lusgs_cpu(param_int, param_real, rop, rop_1, coe, ti,tj,tk,
     &               drodm_ref, drodm_out)
      call system_clock(t0)
      do it = 1, nrep
        call lusgs_cpu(param_int, param_real, rop, rop_1, coe,
     &                 ti,tj,tk, drodm_ref, drodm_out)
      enddo
      call system_clock(t1)
      tcpu = dble(t1-t0)/dble(rate)/nrep
      write(*,'(a,1x,a,es12.4)') 'TIME lusgs_cpu', 'ns_per_cell=',
     &  tcpu/ncell*1.d9
      call checksum(param_int, 'lusgs_cpu_dw', drodm_out)
      call checksum(param_int, 'lusgs_cpu_rop1', rop_1)
      drodm(:)   = drodm_out(:)
      rop_ref(:) = rop_1(:)

c----- LU-SGS GPU (hyperplans, target)
      call lusgs_gpu(param_int, param_real, nplane, ni*nj*nk,
     &               plane_ptr,
     &               plane_cell, rop, rop_1, coe, ti,tj,tk,
     &               drodm_ref, drodm_out)
      call system_clock(t0)
      do it = 1, nrep
        call lusgs_gpu(param_int, param_real, nplane, ni*nj*nk,
     &               plane_ptr,
     &                 plane_cell, rop, rop_1, coe, ti,tj,tk,
     &                 drodm_ref, drodm_out)
      enddo
      call system_clock(t1)
      tgpu = dble(t1-t0)/dble(rate)/nrep
      write(*,'(a,1x,a,es12.4)') 'TIME lusgs_gpu', 'ns_per_cell=',
     &  tgpu/ncell*1.d9
!$OMP TARGET UPDATE FROM(drodm_out, rop_1)
      call checksum(param_int, 'lusgs_gpu_dw', drodm_out)
      call checksum(param_int, 'lusgs_gpu_rop1', rop_1)
      call maxdiff(param_int, 'lusgs_gpu_vs_cpu_dw', drodm, drodm_out)
      call maxdiff(param_int, 'lusgs_gpu_vs_cpu_rop1', rop_ref,
     &             rop_1)

!$OMP TARGET EXIT DATA MAP(DELETE: rop, ti, tj, tk, coe, param_int,
!$OMP&  param_real, plane_ptr, plane_cell, drodm, drodm_ref,
!$OMP&  drodm_out, rop_1)
      end
