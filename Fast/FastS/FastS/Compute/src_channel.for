c***********************************************************************
c     $Date: 2014-03-19 20:08:08 +0100 (mer. 19 mars 2014) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine src_channel(ndom, ithread, param_int, param_real, 
     &                 ind_loop, nitcfg, nitrun, mode,
     &                 forcage, drodm, rop, tk, coe)
c***********************************************************************
c     _U   USER : PECHIER
c     
c     ACT
c     _A    Appel des routines de calcul des flux frontieres aux interfaces
c     _A    du domaine.
c     
c     VAL
c     _V    processeur domaine
c     _V    unsteady/steady
c     
c     OUT
c     _O    drodm
C     
c     CONDITION LIMITES PAR MODIFICATION DE FLUX FRONTIERE
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom,ithread,nitcfg, nitrun, ind_loop(6), param_int(0:*)
      INTEGER_E mode

      REAL_E rop(*), drodm(*), tk(*), param_real(0:*), forcage(2),coe(*)

C     Var loc
      REAL_E surf, f4, vol_dt
      INTEGER_E incmax, l, i,j,k,ne,lij,ltij,lt, v4, lvo

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      v4 = param_int(NDIMDX)*3

      !!! calcul debit
      if (mode.eq.1) then
      !!! 3d cart uniquement
#include "FastC/HPC_LAYER/loop_begin.for"
            surf       = abs( tk(lt))
            forcage(1) = forcage(1) + rop(l)*rop(l+v4)*surf
            forcage(2) = forcage(2) + surf
#include "FastC/HPC_LAYER/loop_end.for"

      !!! ajout forcage
      else

#include "FastC/HPC_LAYER/loop_begin.for"

           vol_dt     = 1./max(coe(l),1.e-15)

           drodm(l+v4)= drodm( l +v4) + param_real(FORCAGE_DEBN)*vol_dt
#include "FastC/HPC_LAYER/loop_end.for"

      endif
      end
