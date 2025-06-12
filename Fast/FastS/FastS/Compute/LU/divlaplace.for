c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 40 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine divlaplace(ndo, it_jacobi, ithread,
     &                     param_int, param_real,
     &                     ind_loop, 
     &                     ti, tj, tk, vol,
     &                     rop, phi1, phi2, rhs, src, residu)
c***********************************************************************
c     _U   USER : PECHIER 
c     _U   USER : DECK
c     ACT
c     _A      Factorisation LU
c     LU+SSOR pour Spalart Allmaras
c     VAL
c     _V      LCI + Jameson-Turkel
c     I/O
c     _/    rop_ssiter
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndo, ithread, it_jacobi, ind_loop(6), param_int(0:*)

      REAL_E phi1(param_int(NDIMDX)),phi2(param_int(NDIMDX))
      REAL_E rhs(param_int(NDIMDX))
      REAL_E ti(param_int(NDIMDX_MTR),2)
      REAL_E tj(param_int(NDIMDX_MTR),2)
      REAL_E tk(param_int(NDIMDX_MTR),1)
      REAL_E vol(NDIMDX_MTR)
      REAL_E rop(param_int(NDIMDX), param_int(NEQ))
      REAL_E src(param_int(NDIMDX), param_int(NEQ)+1)
      REAL_E param_real(0:*), residu(*)

c     Var loc
      INTEGER_E m, k, j, lij, ls, i, l,lt, lvo, ltij, inci,incj,inck
     
      REAL_E dx_inv, dy_inv, dz_inv, b, diag, diff_i, diff_j, diff_k

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

! on blinde si pas assez de travail pour tous les threads

       inci      =1
       incj      =param_int(NIJK)
       inck      =param_int(NIJK)*param_int(NIJK+1)
       !inci_mtr  =param_int(NIJK_MTR)
       !incj_mtr  =param_int(NIJK_MTR+1)
       !inck_mtr  =param_int(NIJK_MTR+2)

      if(ind_loop(2).lt.ind_loop(1)) return

      IF(it_jacobi.eq.1) THEN

#include "FastC/HPC_LAYER/loop_begin.for"
            phi1(l)= 0.
#include "FastC/HPC_LAYER/loop_end.for"

      if (ind_loop(1).eq.1) then
      i =0
#include "FastC/HPC_LAYER/loopPlanI_begin.for"
            phi1(l)= 0.
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif
 
      if (ind_loop(2).eq.param_int(IJKV)) then
      i =param_int(IJKV)+1
#include "FastC/HPC_LAYER/loopPlanI_begin.for"
            phi1(l)= 0.
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif        
 
      if (ind_loop(3).eq.1) then
         j =0
#include "FastC/HPC_LAYER/loopPlanJ_begin.for"
            phi1(l)= 0.
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif     

      if (ind_loop(4).eq.param_int(IJKV+1)) then
         j =param_int(IJKV+1)+1
#include "FastC/HPC_LAYER/loopPlanJ_begin.for"
            phi1(l)= 0.
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif        
 
      if (ind_loop(5).eq.1) then
          k =0
#include "FastC/HPC_LAYER/loopPlanK_begin.for"
            phi1(l)= 0.
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif     

      if (ind_loop(6).eq.param_int(IJKV+2)) then
          k =param_int(IJKV+2)+1
#include "FastC/HPC_LAYER/loopPlanK_begin.for"
            phi1(l)= 0.
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif
  
      ENDIF

#include "FastC/HPC_LAYER/loop_begin.for"
            ! soustraire vitesse!!!
            dx_inv = ti(lt,1)/vol(lt)
            dy_inv = tj(lt,2)/vol(lt)
            dz_inv = tk(lt,1)/vol(lt)
 
            diff_i = ( src( l+inci, 2)-rop(l+inci,2) )*src( l+inci, 6)
     &              -( src( l-inci, 2)-rop(l-inci,2) )*src( l-inci, 6)
            diff_j = ( src( l+incj, 3)-rop(l+incj,3) )*src( l+incj, 6)
     &              -( src( l-incj, 3)-rop(l-incj,3) )*src( l-incj, 6)
            diff_k = ( src( l+inck, 4)-rop(l+inck,4) )*src( l+inck, 6)
     &              -( src( l-inck, 4)-rop(l-inck,4) )*src( l-inck, 6)
           
c            rhs(l) =  0.5* diff_i*ti(lt,1)
c     &              + 0.5* diff_j*tj(lt,2)
c     &              + 0.5* diff_k*tk(lt,1)
c     &              +      (phi1(l+inci) - phi1(l-inci))*ti(lt,1)*dx_inv
c     &              +      (phi1(l+incj) - phi1(l-incj))*tj(lt,2)*dy_inv
c     &              +      (phi1(l+inck) - phi1(l-inck))*tk(lt,1)*dz_inv
c           diag = 0.5*( 1./(dx_inv**2)+ 1./(dy_inv**2) +1./(dz_inv**2) )
c     &          /vol(lt)

            rhs(l) =  0.5* diff_i*dx_inv
     &              + 0.5* diff_j*dy_inv
     &              + 0.5* diff_k*dz_inv
     &       +  (phi1(l+inci) + phi1(l-inci))*dx_inv*dx_inv
     &       +  (phi1(l+incj) + phi1(l-incj))*dy_inv*dy_inv
     &       +  (phi1(l+inck) + phi1(l-inck))*dz_inv*dz_inv


           diag =  2./(dx_inv**2)+ 2./(dy_inv**2) +2./(dz_inv**2) 

           phi2(l) = rhs(l)/diag
           
           residu(ithread) = residu(ithread) + (rhs(l)-phi1(l))**2
c          if(k.eq.20.and.j.le.30.and.l-lij.le.10)
c     & write(*,*) residu(ithread), k,j,l-lij

#include "FastC/HPC_LAYER/loop_end.for"

      if (ind_loop(1).eq.1) then
      i =0
#include "FastC/HPC_LAYER/loopPlanI_begin.for"
            phi2(l)=phi1(l+inci)
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif
 
      if (ind_loop(2).eq.param_int(IJKV)) then
      i =param_int(IJKV)+1
#include "FastC/HPC_LAYER/loopPlanI_begin.for"
            phi2(l)=phi1(l-inci)
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif        
 
      if (ind_loop(3).eq.1) then
         j =0
#include "FastC/HPC_LAYER/loopPlanJ_begin.for"
            phi2(l)=phi1(l+incj)
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif     

      if (ind_loop(4).eq.param_int(IJKV+1)) then
         j =param_int(IJKV+1)+1
#include "FastC/HPC_LAYER/loopPlanJ_begin.for"
            phi2(l)=phi1(l-incj)
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif        
 
      if (ind_loop(5).eq.1) then
          k =0
#include "FastC/HPC_LAYER/loopPlanK_begin.for"
            phi2(l)=phi1(l+inck)
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif     

      if (ind_loop(6).eq.param_int(IJKV+2)) then
          k =param_int(IJKV+2)+1
#include "FastC/HPC_LAYER/loopPlanK_begin.for"
            phi2(l)=phi1(l-inck)
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif

      !if(ithread.eq.1) write(*,*)'RES', residu(1:4)

      end

c***********************************************************************
      subroutine cl_vec(ndo,  ithread,  param_int,
     &                  ind_loop, 
     &                  phi1)
c***********************************************************************
c     _U   USER : PECHIER 
c     _U   USER : DECK
c     ACT
c     _A      Factorisation LU
c     LU+SSOR pour Spalart Allmaras
c     VAL
c     _V      LCI + Jameson-Turkel
c     I/O
c     _/    rop_ssiter
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndo, ithread, ind_loop(6), param_int(0:*)

      REAL_E phi1(param_int(NDIMDX))

c     Var loc
      INTEGER_E m, k, j, lij, ls, i, l,lt, lvo, ltij
     
#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

! on blinde si pas assez de travail pour tous les threads
      if(ind_loop(2).lt.ind_loop(1)) return


c#include "FastC/HPC_LAYER/loop_begin.for"
c            phi1(l)= 0.
c#include "FastC/HPC_LAYER/loop_end.for"

      if (ind_loop(1).eq.1) then
      i =0
#include "FastC/HPC_LAYER/loopPlanI_begin.for"
            phi1(l)= phi1(l+1)
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif
 
      if (ind_loop(2).eq.param_int(IJKV)) then
      i =param_int(IJKV)+1
#include "FastC/HPC_LAYER/loopPlanI_begin.for"
            phi1(l)= phi1(l-1)
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif        
 
      if (ind_loop(3).eq.1) then
         j =0
#include "FastC/HPC_LAYER/loopPlanJ_begin.for"
            phi1(l)=  phi1(l+param_int(NIJK))
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif     

      if (ind_loop(4).eq.param_int(IJKV+1)) then
         j =param_int(IJKV+1)+1
#include "FastC/HPC_LAYER/loopPlanJ_begin.for"
            phi1(l)=  phi1(l-param_int(NIJK))
            !phi1(l)= 0.
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif        
 
      if (ind_loop(5).eq.1) then
          k =0
#include "FastC/HPC_LAYER/loopPlanK_begin.for"
            phi1(l)=  phi1(l+param_int(NIJK)*param_int(NIJK+1))
            !phi1(l)= 0.
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif     

      if (ind_loop(6).eq.param_int(IJKV+2)) then
          k =param_int(IJKV+2)+1
#include "FastC/HPC_LAYER/loopPlanK_begin.for"
            phi1(l)=  phi1(l-param_int(NIJK)*param_int(NIJK+1))
            !phi1(l)= 0.
#include "FastC/HPC_LAYER/loopPlan_end.for"
      endif
  
      end

c***********************************************************************
      subroutine matvec(ndo, ithread,
     &               param_int,  ind_loop, 
     &               ti, tj, tk, vol,
     &               phi1,phi2)
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndo, ithread, ind_loop(6), param_int(0:*)

      REAL_E phi1(param_int(NDIMDX)),phi2(param_int(NDIMDX))
      REAL_E ti(param_int(NDIMDX_MTR),2)
      REAL_E tj(param_int(NDIMDX_MTR),2)
      REAL_E tk(param_int(NDIMDX_MTR),1)
      REAL_E vol(NDIMDX_MTR)

c     Var loc
      INTEGER_E m, k, j, lij, ls, i, l,lt, lvo, ltij, inci,incj,inck
     
      REAL_E dx_inv, dy_inv, dz_inv, b, diag, diff_i, diff_j, diff_k

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

! on blinde si pas assez de travail pour tous les threads

       inci      =1
       incj      =param_int(NIJK)
       inck      =param_int(NIJK)*param_int(NIJK+1)
       !inci_mtr  =param_int(NIJK_MTR)
       !incj_mtr  =param_int(NIJK_MTR+1)
       !inck_mtr  =param_int(NIJK_MTR+2)

      if(ind_loop(2).lt.ind_loop(1)) return

#include "FastC/HPC_LAYER/loop_begin.for"
            ! soustraire vitesse!!!
            dx_inv = ti(lt,1)/vol(lt)
            dy_inv = tj(lt,2)/vol(lt)
            dz_inv = tk(lt,1)/vol(lt)
 
           diag = - 2.*(dx_inv**2)- 2.*(dy_inv**2) -2.*(dz_inv**2) 

            phi2(l) = 
     &         (phi1(l+inci) + phi1(l-inci) -2*phi1(l) )*dx_inv*dx_inv
     &      +  (phi1(l+incj) + phi1(l-incj) -2*phi1(l) )*dy_inv*dy_inv
     &      +  (phi1(l+inck) + phi1(l-inck) -2*phi1(l) )*dz_inv*dz_inv

            phi2(l) = phi2(l)/diag


#include "FastC/HPC_LAYER/loop_end.for"

       end


c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 40 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine rhs(ndo, ithread, param_int, 
     &               ind_loop, 
     &               ti, tj, tk, vol,
     &               rop, rhs_out, src)
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndo, ithread, ind_loop(6), param_int(0:*)

      REAL_E rhs_out(param_int(NDIMDX))
      REAL_E ti(param_int(NDIMDX_MTR),2)
      REAL_E tj(param_int(NDIMDX_MTR),2)
      REAL_E tk(param_int(NDIMDX_MTR),1)
      REAL_E vol(NDIMDX_MTR)
      REAL_E rop(param_int(NDIMDX), param_int(NEQ))
      REAL_E src(param_int(NDIMDX), param_int(NEQ)+1)

c     Var loc
      INTEGER_E m, k, j, lij, ls, i, l,lt, lvo, ltij, inci,incj,inck
     
      REAL_E dx_inv, dy_inv, dz_inv, b, diag, diff_i, diff_j, diff_k

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

! on blinde si pas assez de travail pour tous les threads

       inci      =1
       incj      =param_int(NIJK)
       inck      =param_int(NIJK)*param_int(NIJK+1)

      if(ind_loop(2).lt.ind_loop(1)) return

#include "FastC/HPC_LAYER/loop_begin.for"
            ! soustraire vitesse!!!
            dx_inv = ti(lt,1)/vol(lt)
            dy_inv = tj(lt,2)/vol(lt)
            dz_inv = tk(lt,1)/vol(lt)
 
            diff_i = ( src( l+inci, 2)-rop(l+inci,2) )*src( l+inci, 6)
     &              -( src( l-inci, 2)-rop(l-inci,2) )*src( l-inci, 6)
            diff_j = ( src( l+incj, 3)-rop(l+incj,3) )*src( l+incj, 6)
     &              -( src( l-incj, 3)-rop(l-incj,3) )*src( l-incj, 6)
            diff_k = ( src( l+inck, 4)-rop(l+inck,4) )*src( l+inck, 6)
     &              -( src( l-inck, 4)-rop(l-inck,4) )*src( l-inck, 6)
           
            rhs_out(l) =  0.5* diff_i*dx_inv
     &                  + 0.5* diff_j*dy_inv
     &                  + 0.5* diff_k*dz_inv

c      if (k.eq.5.and.rhs_out(l).ne.0)write(*,*)'zz',rhs_out(l),k,j,l-lij
c          if(k.eq.5.and.j.eq.5.and.l-lij.le.10)
c          if(k.eq.5.and.j.eq.3.and.l-lij.eq.40)
c     & write(*,*)'verif', diff_i, diff_j, diff_k !k,j,l-lij

#include "FastC/HPC_LAYER/loop_end.for"

      end


c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 40 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine vec_add(ndo, ithread,
     &                   param_int,
     &                   ind_loop, 
     &                   vec1, vec2, alpha, beta)
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndo, ithread, ind_loop(6), param_int(0:*)

      REAL_E vec1(param_int(NDIMDX)), vec2(param_int(NDIMDX)), alpha,
     & beta

c     Var loc
      INTEGER_E m, k, j, lij, ls, i, l,lt, lvo, ltij, inci,incj,inck
     

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

! on blinde si pas assez de travail pour tous les threads

      if(ind_loop(2).lt.ind_loop(1)) return

#include "FastC/HPC_LAYER/loop_begin.for"
         vec1(l)=vec1(l)*beta + alpha*vec2(l)
#include "FastC/HPC_LAYER/loop_end.for"

      end

c***********************************************************************
      subroutine dot_product(ndo, ithread, param_int, 
     &                      ind_loop, 
     &                     vec1, vec2, residu)
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndo, ithread, ind_loop(6), param_int(0:*)

      REAL_E vec1(param_int(NDIMDX)), vec2(param_int(NDIMDX))
      REAL_E residu(*)

c     Var loc
      INTEGER_E m, k, j, lij, ls, i, l,lt, lvo, ltij, inci,incj,inck
     

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

! on blinde si pas assez de travail pour tous les threads

       inci      =1
       incj      =param_int(NIJK)
       inck      =param_int(NIJK)*param_int(NIJK+1)

      if(ind_loop(2).lt.ind_loop(1)) return

#include "FastC/HPC_LAYER/loop_begin.for"
          residu(ithread)=residu(ithread)+ vec1(l)*vec2(l)
#include "FastC/HPC_LAYER/loop_end.for"

      end
