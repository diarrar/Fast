c***********************************************************************
c     LU-SGS 5 equations, maillage 3dfull, formulation GPU par
c     hyperplans i+j+k = cte : boucle sequentielle sur les plans
c     (host), noyau OpenMP target parallele sur les cellules du plan.
c     Les listes de cellules par plan sont construites sur le host
c     (build_planes) et residentes sur le device.
c***********************************************************************
      subroutine build_planes(ni, nj, nk, nx, ny, nplane, plane_ptr,
     &                        plane_cell)
      implicit none
      integer ni,nj,nk,nx,ny,nplane
      integer plane_ptr(*), plane_cell(*)
      integer i,j,k,p,m,l
      nplane = ni+nj+nk-2
      do p = 1, nplane+1
        plane_ptr(p) = 0
      enddo
c     comptage
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
        p = i+j+k-2
        plane_ptr(p+1) = plane_ptr(p+1) + 1
      enddo
      enddo
      enddo
      plane_ptr(1) = 1
      do p = 1, nplane
        plane_ptr(p+1) = plane_ptr(p+1) + plane_ptr(p)
      enddo
c     remplissage (ordre k,j,i dans chaque plan)
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
        p = i+j+k-2
        m = plane_ptr(p)
        l = 1 + (i+1) + (j+1)*nx + (k+1)*nx*ny
        plane_cell(m) = l
        plane_ptr(p)  = m + 1
      enddo
      enddo
      enddo
      do p = nplane, 1, -1
        plane_ptr(p+1) = plane_ptr(p)
      enddo
      plane_ptr(1) = 1
      end

c***********************************************************************
      subroutine lusgs_gpu(param_int, param_real, nplane, ncell,
     &                     plane_ptr,
     &                     plane_cell, rop, rop_1, coe,
     &                     ti,tj,tk, drodm_in, drodm_out)
      implicit none
#include "mini_param.h"
      INTEGER_E param_int(0:NPARAM-1)
      REAL_E    param_real(0:NPARAM-1)
      integer   nplane, ncell, plane_ptr(nplane+1), plane_cell(ncell)
      REAL_E rop      (param_int(NDIMDX),param_int(NEQ))
      REAL_E rop_1    (param_int(NDIMDX),param_int(NEQ))
      REAL_E coe      (param_int(NDIMDX),param_int(NEQ_COE))
      REAL_E drodm_in (param_int(NDIMDX),param_int(NEQ))
      REAL_E drodm_out(param_int(NDIMDX),param_int(NEQ))
      REAL_E ti(param_int(NDIMDX_MTR),3), tj(param_int(NDIMDX_MTR),3)
      REAL_E tk(param_int(NDIMDX_MTR),3)
c     variables locales
      INTEGER_E i,j,k,l,lt,m,p,lm,ii,jj,kk,rem,nx,ny,
     & inci,incj,inck,inci_mtr,incj_mtr,inck_mtr,ni,nj,nk
#include "lu_locvar.h"
#include "formule_mtr_param.h"

      ni   = param_int(IJKV)
      nj   = param_int(IJKV+1)
      nk   = param_int(IJKV+2)
      nx   = param_int(NIJK)
      ny   = param_int(NIJK+1)
      inci = 1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)
      inci_mtr = param_int(NIJK_MTR)
      incj_mtr = param_int(NIJK_MTR+1)
      inck_mtr = param_int(NIJK_MTR+2)
#include "lu_const.h"

c-----------------------------------------------------------------------
c     balayage inferieur par hyperplans croissants
c-----------------------------------------------------------------------
      signe = 0.5
      DO p = 1, nplane
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
!$OMP& PRIVATE(l,lt,lm,ii,jj,kk,rem,i,j,k,l1,lt1,lt2,idir,xal,diag,rad,
!$OMP& b11,b12,b13,b14,b15,b21,b22,b23,b24,b25,b31,b32,b33,b34,b35,
!$OMP& b41,b42,b43,b44,b45,b51,b52,b53,b54,b55,b1,b2,b3,b4,b5,
!$OMP& x1,x2,x3,x4,x5,u,v,w,t,q2,h,ph2,qn,tcx,tcy,tcz)
       do m = plane_ptr(p), plane_ptr(p+1)-1
        l  = plane_cell(m)
        lm = l - 1
        kk = lm/(nx*ny)
        rem= lm - kk*nx*ny
        jj = rem/nx
        ii = rem - jj*nx
        i  = ii - 1
        j  = jj - 1
        k  = kk - 1
        lt = indmtr(i,j,k)
        xal = coe(l,1)*signe
        x1 = drodm_in(l,1)
        x2 = drodm_in(l,2)
        x3 = drodm_in(l,3)
        x4 = drodm_in(l,4)
        x5 = drodm_in(l,5)
        if (k.gt.1) then
          l1   = l  - inck
          lt1  = lt - inck_mtr
          lt2  = lt
          idir = 4
#define TN tk
#include "lu_mat.h"
#undef TN
          x1 = x1 + b1*xal
          x2 = x2 + b2*xal
          x3 = x3 + b3*xal
          x4 = x4 + b4*xal
          x5 = x5 + b5*xal
        endif
        if (j.gt.1) then
          l1   = l  - incj
          lt1  = lt - incj_mtr
          lt2  = lt
          idir = 3
#define TN tj
#include "lu_mat.h"
#undef TN
          x1 = x1 + b1*xal
          x2 = x2 + b2*xal
          x3 = x3 + b3*xal
          x4 = x4 + b4*xal
          x5 = x5 + b5*xal
        endif
        if (i.gt.1) then
          l1   = l  - inci
          lt1  = lt - inci_mtr
          lt2  = lt
          idir = 2
#define TN ti
#include "lu_mat.h"
#undef TN
          x1 = x1 + b1*xal
          x2 = x2 + b2*xal
          x3 = x3 + b3*xal
          x4 = x4 + b4*xal
          x5 = x5 + b5*xal
        endif
        diag = 1./coe(l,5)
        drodm_out(l,1) = x1*diag
        drodm_out(l,2) = x2*diag
        drodm_out(l,3) = x3*diag
        drodm_out(l,4) = x4*diag
        drodm_out(l,5) = x5*diag
       enddo
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
      ENDDO

c-----------------------------------------------------------------------
c     balayage superieur par hyperplans decroissants + Newton
c-----------------------------------------------------------------------
      signe = -0.5
      DO p = nplane, 1, -1
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
!$OMP& PRIVATE(l,lt,lm,ii,jj,kk,rem,i,j,k,l1,lt1,lt2,idir,xal,diag,rad,
!$OMP& b11,b12,b13,b14,b15,b21,b22,b23,b24,b25,b31,b32,b33,b34,b35,
!$OMP& b41,b42,b43,b44,b45,b51,b52,b53,b54,b55,b1,b2,b3,b4,b5,
!$OMP& x1,x2,x3,x4,x5,u,v,w,t,q2,h,ph2,qn,tcx,tcy,tcz,
!$OMP& ro_old,u_old,v_old,w_old,t_old,roe_old,r_1)
       do m = plane_ptr(p), plane_ptr(p+1)-1
        l  = plane_cell(m)
        lm = l - 1
        kk = lm/(nx*ny)
        rem= lm - kk*nx*ny
        jj = rem/nx
        ii = rem - jj*nx
        i  = ii - 1
        j  = jj - 1
        k  = kk - 1
        lt = indmtr(i,j,k)
        xal  = coe(l,1)*signe
        diag = coe(l,5)
        x1 = drodm_out(l,1)*diag
        x2 = drodm_out(l,2)*diag
        x3 = drodm_out(l,3)*diag
        x4 = drodm_out(l,4)*diag
        x5 = drodm_out(l,5)*diag
        if (k.lt.nk) then
          l1   = l  + inck
          lt1  = lt + inck_mtr
          lt2  = lt + 2*inck_mtr
          idir = 4
#define TN tk
#include "lu_mat.h"
#undef TN
          x1 = x1 + b1*xal
          x2 = x2 + b2*xal
          x3 = x3 + b3*xal
          x4 = x4 + b4*xal
          x5 = x5 + b5*xal
        endif
        if (j.lt.nj) then
          l1   = l  + incj
          lt1  = lt + incj_mtr
          lt2  = lt + 2*incj_mtr
          idir = 3
#define TN tj
#include "lu_mat.h"
#undef TN
          x1 = x1 + b1*xal
          x2 = x2 + b2*xal
          x3 = x3 + b3*xal
          x4 = x4 + b4*xal
          x5 = x5 + b5*xal
        endif
        if (i.lt.ni) then
          l1   = l  + inci
          lt1  = lt + inci_mtr
          lt2  = lt + 2*inci_mtr
          idir = 2
#define TN ti
#include "lu_mat.h"
#undef TN
          x1 = x1 + b1*xal
          x2 = x2 + b2*xal
          x3 = x3 + b3*xal
          x4 = x4 + b4*xal
          x5 = x5 + b5*xal
        endif
        diag = 1./coe(l,5)
        drodm_out(l,1) = x1*diag
        drodm_out(l,2) = x2*diag
        drodm_out(l,3) = x3*diag
        drodm_out(l,4) = x4*diag
        drodm_out(l,5) = x5*diag
#include "mjr_newton.h"
       enddo
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
      ENDDO
      end
