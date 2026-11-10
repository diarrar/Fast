c***********************************************************************
c     Flux de Roe o3 Euler 3dfull, formulation GPU "gather":
c     chaque cellule calcule ses 6 flux de face et ecrit drodm sans
c     dependance entre cellules (pas de course d ecriture).
c     Un seul noyau OpenMP target sur toute la zone, boucles i,j,k
c     collapsees (cf. loopGpu_begin.for de la premiere tentative).
c***********************************************************************
      subroutine flux_roe_gpu(param_int, param_real, rop, ti,tj,tk,
     &                        drodm)
      implicit none
#include "mini_param.h"
      INTEGER_E param_int(0:NPARAM-1)
      REAL_E    param_real(0:NPARAM-1)
      REAL_E rop  ( param_int(NDIMDX)*param_int(NEQ) )
      REAL_E drodm( param_int(NDIMDX)*param_int(NEQ) )
      REAL_E ti( param_int(NDIMDX_MTR)*3 ), tj( param_int(NDIMDX_MTR)*3 )
      REAL_E tk( param_int(NDIMDX_MTR)*3 )
c     variables locales
      INTEGER_E i,j,k,l,lt,inci,incj,inck,inci_mtr,incj_mtr,inck_mtr,
     & ni,nj,nk
      REAL_E d1,d2,d3,d4,d5
#include "flux_locvar.h"
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
#include "flux_const.h"

!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3)
!$OMP& PRIVATE(l,lt,lf,ltf,inc,nm,nm2,np,d1,d2,d3,d4,d5,
!$OMP& si,tcx,tcy,tcz,r1,p1,h1,r2,p2,h2,r,r_1,u,v,w,h,q,c,
!$OMP& qm1,qm2,qm3,qm4,qm5,qp1,qp2,qp3,qp4,qp5,qn1,qn2,du,dv,dw,dp,
!$OMP& dqn,s_1,nx,ny,nz,qn,f1,f2,f3,f4,rou1,rou2,p1p2,
!$OMP& flu1,flu2,flu3,flu4,flu5)
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
        l  = inddm(i,j,k)
        lt = indmtr(i,j,k)
c       direction i: face i-1/2 (+) puis face i+1/2 (-)
        inc = inci
        lf  = l
        ltf = lt
#define TN ti
#include "flux_face_roe_o3.h"
        d1 = flu1
        d2 = flu2
        d3 = flu3
        d4 = flu4
        d5 = flu5
        lf  = l  + inci
        ltf = lt + inci_mtr
#include "flux_face_roe_o3.h"
#undef TN
        d1 = d1 - flu1
        d2 = d2 - flu2
        d3 = d3 - flu3
        d4 = d4 - flu4
        d5 = d5 - flu5
c       direction j
        inc = incj
        lf  = l
        ltf = lt
#define TN tj
#include "flux_face_roe_o3.h"
        d1 = d1 + flu1
        d2 = d2 + flu2
        d3 = d3 + flu3
        d4 = d4 + flu4
        d5 = d5 + flu5
        lf  = l  + incj
        ltf = lt + incj_mtr
#include "flux_face_roe_o3.h"
#undef TN
        d1 = d1 - flu1
        d2 = d2 - flu2
        d3 = d3 - flu3
        d4 = d4 - flu4
        d5 = d5 - flu5
c       direction k
        inc = inck
        lf  = l
        ltf = lt
#define TN tk
#include "flux_face_roe_o3.h"
        d1 = d1 + flu1
        d2 = d2 + flu2
        d3 = d3 + flu3
        d4 = d4 + flu4
        d5 = d5 + flu5
        lf  = l  + inck
        ltf = lt + inck_mtr
#include "flux_face_roe_o3.h"
#undef TN
        drodm(l +v1) = d1 - flu1
        drodm(l +v2) = d2 - flu2
        drodm(l +v3) = d3 - flu3
        drodm(l +v4) = d4 - flu4
        drodm(l +v5) = d5 - flu5
      enddo
      enddo
      enddo
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
      end
