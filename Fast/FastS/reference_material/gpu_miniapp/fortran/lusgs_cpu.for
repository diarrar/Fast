c***********************************************************************
c     LU-SGS 5 equations, maillage 3dfull, formulation CPU FastS:
c     balayage lexicographique inferieur puis superieur
c     (transcription de invlu_l.for / invlu_u.for sans sous-domaines).
c     drodm_out = (D+U)^-1 D (D+L)^-1 drodm_in ; rop_1 = Newton(rop)
c***********************************************************************
      subroutine lusgs_cpu(param_int, param_real, rop, rop_1, coe,
     &                     ti,tj,tk, drodm_in, drodm_out)
      implicit none
#include "mini_param.h"
      INTEGER_E param_int(0:*)
      REAL_E    param_real(0:*)
      REAL_E rop      (param_int(NDIMDX),param_int(NEQ))
      REAL_E rop_1    (param_int(NDIMDX),param_int(NEQ))
      REAL_E coe      (param_int(NDIMDX),param_int(NEQ_COE))
      REAL_E drodm_in (param_int(NDIMDX),param_int(NEQ))
      REAL_E drodm_out(param_int(NDIMDX),param_int(NEQ))
      REAL_E ti(param_int(NDIMDX_MTR),3), tj(param_int(NDIMDX_MTR),3)
      REAL_E tk(param_int(NDIMDX_MTR),3)
c     variables locales
      INTEGER_E i,j,k,l,lt,lij,ltij,inci,incj,inck,inci_mtr,incj_mtr,
     & inck_mtr,ni,nj,nk
#include "lu_locvar.h"
#include "formule_param.h"
#include "formule_mtr_param.h"

      ni   = param_int(IJKV)
      nj   = param_int(IJKV+1)
      nk   = param_int(IJKV+2)
      inci = 1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)
      inci_mtr = param_int(NIJK_MTR)
      incj_mtr = param_int(NIJK_MTR+1)
      inck_mtr = param_int(NIJK_MTR+2)
#include "lu_const.h"

c-----------------------------------------------------------------------
c     balayage inferieur : (D+L) x = b
c-----------------------------------------------------------------------
      signe = 0.5
      DO k = 1, nk
       DO j = 1, nj
        lij  =       inddm(1 , j, k) -1
        ltij = lij - indmtr(1 , j, k) +1
c       contributions k-1 et j-1 (lignes/plans deja finalises)
!$OMP SIMD
        do l = lij+1, lij+ni
          lt  = l - ltij
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
          drodm_out(l,1) = x1
          drodm_out(l,2) = x2
          drodm_out(l,3) = x3
          drodm_out(l,4) = x4
          drodm_out(l,5) = x5
        enddo
c       recurrence en i
        do i = 1, ni
          l   = lij + i
          lt  = l - ltij
          xal = coe(l,1)*signe
          if (i.gt.1) then
            l1   = l  - inci
            lt1  = lt - inci_mtr
            lt2  = lt
            idir = 2
#define TN ti
#include "lu_mat.h"
#undef TN
            drodm_out(l,1) = drodm_out(l,1) + b1*xal
            drodm_out(l,2) = drodm_out(l,2) + b2*xal
            drodm_out(l,3) = drodm_out(l,3) + b3*xal
            drodm_out(l,4) = drodm_out(l,4) + b4*xal
            drodm_out(l,5) = drodm_out(l,5) + b5*xal
          endif
          diag = 1./coe(l,5)
          drodm_out(l,1) = drodm_out(l,1)*diag
          drodm_out(l,2) = drodm_out(l,2)*diag
          drodm_out(l,3) = drodm_out(l,3)*diag
          drodm_out(l,4) = drodm_out(l,4)*diag
          drodm_out(l,5) = drodm_out(l,5)*diag
        enddo
       ENDDO
      ENDDO

c-----------------------------------------------------------------------
c     balayage superieur : (D+U) x = D y, puis mise a jour Newton
c-----------------------------------------------------------------------
      signe = -0.5
      DO k = nk, 1, -1
       DO j = nj, 1, -1
        lij  =       inddm(1 , j, k) -1
        ltij = lij - indmtr(1 , j, k) +1
!$OMP SIMD
        do l = lij+1, lij+ni
          lt  = l - ltij
          xal = coe(l,1)*signe
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
          drodm_out(l,1) = x1
          drodm_out(l,2) = x2
          drodm_out(l,3) = x3
          drodm_out(l,4) = x4
          drodm_out(l,5) = x5
        enddo
c       recurrence en i decroissant
        do i = ni, 1, -1
          l   = lij + i
          lt  = l - ltij
          xal = coe(l,1)*signe
          if (i.lt.ni) then
            l1   = l  + inci
            lt1  = lt + inci_mtr
            lt2  = lt + 2*inci_mtr
            idir = 2
#define TN ti
#include "lu_mat.h"
#undef TN
            drodm_out(l,1) = drodm_out(l,1) + b1*xal
            drodm_out(l,2) = drodm_out(l,2) + b2*xal
            drodm_out(l,3) = drodm_out(l,3) + b3*xal
            drodm_out(l,4) = drodm_out(l,4) + b4*xal
            drodm_out(l,5) = drodm_out(l,5) + b5*xal
          endif
          diag = 1./coe(l,5)
          drodm_out(l,1) = drodm_out(l,1)*diag
          drodm_out(l,2) = drodm_out(l,2)*diag
          drodm_out(l,3) = drodm_out(l,3)*diag
          drodm_out(l,4) = drodm_out(l,4)*diag
          drodm_out(l,5) = drodm_out(l,5)*diag
        enddo
c       Newton
!$OMP SIMD
        do l = lij+1, lij+ni
#include "mjr_newton.h"
        enddo
       ENDDO
      ENDDO
      end
