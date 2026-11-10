c***********************************************************************
c     Flux de Roe o3 Euler 3dfull, formulation CPU FastS ("scatter"):
c     chaque cellule calcule le flux de sa face gauche (l-1/2) et le
c     distribue sur les deux cellules voisines, puis faces manquantes
c     en imax+1, jmax+1, kmax+1 (cf. fluroe_euler_o3_3dfull.for).
c***********************************************************************
      subroutine flux_roe_cpu(param_int, param_real, rop, ti,tj,tk,
     &                        drodm)
      implicit none
#include "mini_param.h"
      INTEGER_E param_int(0:*)
      REAL_E    param_real(0:*)
      REAL_E rop  ( param_int(NDIMDX)*param_int(NEQ) )
      REAL_E drodm( param_int(NDIMDX)*param_int(NEQ) )
      REAL_E ti( param_int(NDIMDX_MTR)*3 ), tj( param_int(NDIMDX_MTR)*3 )
      REAL_E tk( param_int(NDIMDX_MTR)*3 )
c     variables locales
      INTEGER_E i,j,k,l,lt,l0,lij,ltij,inci,incj,inck,inci_mtr,
     & incj_mtr,inck_mtr,ni,nj,nk,ndim
#include "flux_locvar.h"
#include "formule_param.h"
#include "formule_mtr_param.h"

      ni   = param_int(IJKV)
      nj   = param_int(IJKV+1)
      nk   = param_int(IJKV+2)
      ndim = param_int(NDIMDX)
      inci = 1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)
      inci_mtr = param_int(NIJK_MTR)
      incj_mtr = param_int(NIJK_MTR+1)
      inck_mtr = param_int(NIJK_MTR+2)
#include "flux_const.h"

c     init_rhs: drodm = 0 (fait par src_term dans FastS)
      do l = 1, ndim*5
        drodm(l) = 0.
      enddo

      DO k = 1, nk
       DO j = 1, nj
        lij  =       inddm(1 , j, k) -1
        ltij = lij - indmtr(1 , j, k) +1
c       faces k-1/2
!$OMP SIMD
        do l = lij+1, lij+ni
          lt  = l - ltij
          lf  = l
          ltf = lt
          inc = inck
          l0  = l - inck
#define TN tk
#include "flux_face_roe_o3.h"
#undef TN
          drodm(l  +v1) = drodm(l  +v1) + flu1
          drodm(l0 +v1) = drodm(l0 +v1) - flu1
          drodm(l  +v2) = drodm(l  +v2) + flu2
          drodm(l0 +v2) = drodm(l0 +v2) - flu2
          drodm(l  +v3) = drodm(l  +v3) + flu3
          drodm(l0 +v3) = drodm(l0 +v3) - flu3
          drodm(l  +v4) = drodm(l  +v4) + flu4
          drodm(l0 +v4) = drodm(l0 +v4) - flu4
          drodm(l  +v5) = drodm(l  +v5) + flu5
          drodm(l0 +v5) = drodm(l0 +v5) - flu5
        enddo
c       faces j-1/2
!$OMP SIMD
        do l = lij+1, lij+ni
          lt  = l - ltij
          lf  = l
          ltf = lt
          inc = incj
          l0  = l - incj
#define TN tj
#include "flux_face_roe_o3.h"
#undef TN
          drodm(l  +v1) = drodm(l  +v1) + flu1
          drodm(l0 +v1) = drodm(l0 +v1) - flu1
          drodm(l  +v2) = drodm(l  +v2) + flu2
          drodm(l0 +v2) = drodm(l0 +v2) - flu2
          drodm(l  +v3) = drodm(l  +v3) + flu3
          drodm(l0 +v3) = drodm(l0 +v3) - flu3
          drodm(l  +v4) = drodm(l  +v4) + flu4
          drodm(l0 +v4) = drodm(l0 +v4) - flu4
          drodm(l  +v5) = drodm(l  +v5) + flu5
          drodm(l0 +v5) = drodm(l0 +v5) - flu5
        enddo
c       faces i-1/2
!$OMP SIMD
        do l = lij+1, lij+ni
          lt  = l - ltij
          lf  = l
          ltf = lt
          inc = inci
          l0  = l - inci
#define TN ti
#include "flux_face_roe_o3.h"
#undef TN
          drodm(l  +v1) = drodm(l  +v1) + flu1
          drodm(l0 +v1) = drodm(l0 +v1) - flu1
          drodm(l  +v2) = drodm(l  +v2) + flu2
          drodm(l0 +v2) = drodm(l0 +v2) - flu2
          drodm(l  +v3) = drodm(l  +v3) + flu3
          drodm(l0 +v3) = drodm(l0 +v3) - flu3
          drodm(l  +v4) = drodm(l  +v4) + flu4
          drodm(l0 +v4) = drodm(l0 +v4) - flu4
          drodm(l  +v5) = drodm(l  +v5) + flu5
          drodm(l0 +v5) = drodm(l0 +v5) - flu5
        enddo
c       flux manquant en imax+1
        lf  = inddm( ni+1, j, k)
        ltf = indmtr(ni+1, j, k)
        inc = inci
        l0  = lf - inci
#define TN ti
#include "flux_face_roe_o3.h"
#undef TN
        drodm(l0 +v1) = drodm(l0 +v1) - flu1
        drodm(lf +v1) = drodm(lf +v1) + flu1
        drodm(l0 +v2) = drodm(l0 +v2) - flu2
        drodm(lf +v2) = drodm(lf +v2) + flu2
        drodm(l0 +v3) = drodm(l0 +v3) - flu3
        drodm(lf +v3) = drodm(lf +v3) + flu3
        drodm(l0 +v4) = drodm(l0 +v4) - flu4
        drodm(lf +v4) = drodm(lf +v4) + flu4
        drodm(l0 +v5) = drodm(l0 +v5) - flu5
        drodm(lf +v5) = drodm(lf +v5) + flu5
       ENDDO
c      complement flux j en jmax+1
       lij  =       inddm(1 , nj+1, k) -1
       ltij = lij - indmtr(1 , nj+1, k) +1
!$OMP SIMD
       do l = lij+1, lij+ni
         lt  = l - ltij
         lf  = l
         ltf = lt
         inc = incj
         l0  = l - incj
#define TN tj
#include "flux_face_roe_o3.h"
#undef TN
         drodm(l0 +v1) = drodm(l0 +v1) - flu1
         drodm(l  +v1) = drodm(l  +v1) + flu1
         drodm(l0 +v2) = drodm(l0 +v2) - flu2
         drodm(l  +v2) = drodm(l  +v2) + flu2
         drodm(l0 +v3) = drodm(l0 +v3) - flu3
         drodm(l  +v3) = drodm(l  +v3) + flu3
         drodm(l0 +v4) = drodm(l0 +v4) - flu4
         drodm(l  +v4) = drodm(l  +v4) + flu4
         drodm(l0 +v5) = drodm(l0 +v5) - flu5
         drodm(l  +v5) = drodm(l  +v5) + flu5
       enddo
      ENDDO
c     complement flux k en kmax+1
      do j = 1, nj
       lij  =       inddm(1 , j, nk+1) -1
       ltij = lij - indmtr(1 , j, nk+1) +1
!$OMP SIMD
       do l = lij+1, lij+ni
         lt  = l - ltij
         lf  = l
         ltf = lt
         inc = inck
         l0  = l - inck
#define TN tk
#include "flux_face_roe_o3.h"
#undef TN
         drodm(l0 +v1) = drodm(l0 +v1) - flu1
         drodm(l  +v1) = drodm(l  +v1) + flu1
         drodm(l0 +v2) = drodm(l0 +v2) - flu2
         drodm(l  +v2) = drodm(l  +v2) + flu2
         drodm(l0 +v3) = drodm(l0 +v3) - flu3
         drodm(l  +v3) = drodm(l  +v3) + flu3
         drodm(l0 +v4) = drodm(l0 +v4) - flu4
         drodm(l  +v4) = drodm(l  +v4) + flu4
         drodm(l0 +v5) = drodm(l0 +v5) - flu5
         drodm(l  +v5) = drodm(l  +v5) + flu5
       enddo
      enddo
      end
