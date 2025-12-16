/*    
    Copyright 2013-2024 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

# include "FastS/fastS.h"
# include "FastC/fastc.h"
# include "FastS/param_solver.h"
# include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _MPI
#include <mpi.h>
#endif

using namespace K_FLD;
using namespace std;

#undef Conservatif
//#define Conservatif
#undef TimeShow
//#define TimeShow

#ifdef TimeShow

#include <iostream>
#include <fstream>
#include <string> 
#include <iomanip>
#include <sstream>

E_Float time_COM=0.0;
E_Float time_init;
#endif

//E_Float time_COM=0.0;
//E_Float time_init;


void K_FASTS::laplace( 
  E_Int**& param_int  , E_Float**& param_real , E_Int& nidom        , E_Int& nssiter   ,  E_Int& threadmax_sdm, 
  E_Int* ipt_ijkv_sdm , E_Int* ipt_ind_dm_omp , E_Int* ipt_topology , E_Int* ipt_ind_CL,
  E_Int** ipt_ind_dm  , E_Int* iptdtloc       , 
  E_Float**  iptx     , E_Float**  ipty       , E_Float**    iptz   ,
  E_Float**& iptro    , E_Float**& iptro_p1   , 
  E_Float*  iptphi1   , E_Float*   iptphi2    , E_Float*   iptb     ,
  E_Float**  ipti     , E_Float**  iptj       , E_Float** iptk      , E_Float** iptvol, E_Float** iptsrc, E_Float* residu )

 {
      E_Int omp_mode = iptdtloc[8];
      E_Int shift_omp= iptdtloc[11];
      E_Int* ipt_omp = iptdtloc + shift_omp;

      E_Int nitcfg=1;

      E_Int npass         = 0;
      E_Int ibord_ale     = 1;      // on autorise un calcul optimisee des vitesse entrainement en explicit
      E_Int nptpsi        = 1;
      E_Int balance       = 0;

      //
      // choix du tableau de travail en fonction du schema et sous-iteration
      //
      ////  Verifier iptro_CL pour nitcfg> 1 en implicite
      //
      //
      E_Float** iptro_ssiter;
      E_Float** iptro_CL; 

      E_Int nbtask = ipt_omp[nitcfg-1]; 
      E_Int ptiter = ipt_omp[nssiter+ nitcfg-1];

/****************************************************
----- Debut zone // omp
****************************************************/

#ifdef _OPENMP
  E_Int Nthread_max  = omp_get_max_threads();
#else
  E_Int Nthread_max  = 1;
#endif


#pragma omp parallel default(shared)
  {
#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num() +1;
    E_Int  Nbre_thread_actif = omp_get_num_threads();
    E_Float rhs_begin        = omp_get_wtime();
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
    E_Float rhs_begin       = 0;
#endif
   //E_Int Nbre_socket   = NBR_SOCKET;                       // nombre de proc (socket) sur le noeud a memoire partagee
   E_Int Nbre_socket   = 1;                       // nombre de proc (socket) sur le noeud a memoire partagee
   if( Nbre_thread_actif < Nbre_socket) Nbre_socket = 1;

   E_Int Nbre_thread_actif_loc, ithread_loc;

   E_Int thread_parsock  =  Nbre_thread_actif/Nbre_socket;
   E_Int socket          = (ithread-1)/thread_parsock +1;
   E_Int  ithread_sock   = ithread-(socket-1)*thread_parsock;

   E_Int* ipt_topology_socket    = ipt_topology       + (ithread-1)*3;
   E_Int* ipt_ijkv_sdm_thread    = ipt_ijkv_sdm       + (ithread-1)*3;
   E_Int* ipt_ind_CL_thread      = ipt_ind_CL         + (ithread-1)*6;
   E_Int* ipt_ind_CL119          = ipt_ind_CL         + (ithread-1)*6 +  6*Nbre_thread_actif;
   E_Int* ipt_ind_CLgmres        = ipt_ind_CL         + (ithread-1)*6 + 12*Nbre_thread_actif;
   E_Int* ipt_ind_dm_socket      = ipt_ind_dm_omp     + (ithread-1)*12;
   E_Int* ipt_ind_dm_omp_thread  = ipt_ind_dm_socket  + 6;

   E_Int* ipt_nidom_loc, nb_subzone;
   /****************************************************
      -----Boucle sous-iteration
    ****************************************************/
    E_Int nbtask = ipt_omp[nitcfg-1]; 
    E_Int ptiter = ipt_omp[nssiter+ nitcfg-1];
   
    // -----Init BICGSTAB
    // ---------------------------------------------------------------------
    E_Int shift_zone=0;

    E_Float rho1=1;
    E_Float alpha1=1;
    E_Float w1=1;
    E_Int ndimt=0;
    for (E_Int ntask = 0; ntask < nbtask; ntask++)
      {
        E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
        E_Int nd = ipt_omp[ pttask ];
        ndimt = ndimt + param_int[nd][ NDIMDX ];
      }

   printf("ndilt: %d \n", ndimt);
        
    for (E_Int ntask = 0; ntask < nbtask; ntask++)
      {
        E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
        E_Int nd = ipt_omp[ pttask ];

        shift_zone=0;
        for (E_Int n = 0; n < nd; n++)
        {
         shift_zone = shift_zone + param_int[n][ NDIMDX ];
        }

       ithread_loc              = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
       Nbre_thread_actif_loc    = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
       E_Int* ipt_ind_dm_thread = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;

       if (ithread_loc == -1) {continue;}

       // phi1=0
       E_Float* phi_in  = iptphi1 +  shift_zone;
       E_Float* iptr0   = iptb    +  shift_zone;
       E_Float* iptr    = iptb    +  shift_zone + ndimt;
       E_Float* iptp    = iptb    +  shift_zone + ndimt*2;
       E_Float* iptv    = iptb    +  shift_zone + ndimt*3;

       E_Float beta=0;
       E_Float alpha=0;
       vec_add_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread ,  phi_in, phi_in , alpha, beta);
       cl_vec_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , phi_in);
       //initvec_(nd,  ithread_loc, param_int[nd] , ipt_ind_dm_thread , phi_in );
       
       // iptb=div(f)
       rhs_(nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , 
            ipti[nd] , iptj[nd] , iptk[nd]  , iptvol[nd],
            iptro[nd], iptr0 , iptsrc[nd]);  

       //r = b -Ax = b
       beta=0;
       alpha=1;
       vec_add_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptr , iptr0 , alpha, beta);


       //p=r
       beta=0;
       alpha=1;
       //alpha=0;
       vec_add_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptp , iptr , alpha, beta);

       //v=0
       beta=0;
       alpha=0;
       vec_add_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptv , iptv , alpha, beta);

       cl_vec_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptp);
       cl_vec_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptv);
       //initvec_(nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptp );
       //initvec_(nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptv );

       E_Float* res_inout = residu +nd*Nbre_thread_actif;
       res_inout[ithread_loc-1]=0.;

       dot_product_(nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptr0 , iptr0 ,   res_inout);
      }
   #pragma omp barrier
   rho1=0.;
   //E_Float* res_inout = residu +nd*Nbre_thread_actif;
   // valid pour une zone
   E_Float* res_inout = residu;
   for (E_Int ith = 0; ith < Nbre_thread_actif; ith++){rho1=rho1+res_inout[ith];}
   printf("residu Jac init: %g \n", rho1);

    for (E_Int it_jacobi = 1; it_jacobi < 25; it_jacobi++)
    {
       E_Int shift_zone=0;

       for (E_Int ntask = 0; ntask < nbtask; ntask++)
         {
             E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
             E_Int nd = ipt_omp[ pttask ];


             shift_zone=0;
             for (E_Int n = 0; n < nd; n++)
             {
              shift_zone = shift_zone + param_int[n][ NDIMDX ];
             }

            ithread_loc              = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
            Nbre_thread_actif_loc    = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
            E_Int* ipt_ind_dm_thread = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;

            E_Float* res_inout = residu +nd*Nbre_thread_actif;
            res_inout[ithread_loc-1]=0.;

            if (ithread_loc == -1) {continue;}
 
            //calcul de rho(it_jacobi) : r(it_jacobi=0).r(it_jacobi) --> stocké dans rho2

            E_Float* phi_in= iptphi1 +  shift_zone;
            E_Float* iptsol= phi_in;
            E_Float* iptr0 = iptb +  shift_zone;
            E_Float* iptr  = iptb +  shift_zone + ndimt;
            E_Float* iptp  = iptb +  shift_zone + ndimt*2;
            E_Float* iptv  = iptb +  shift_zone + ndimt*3;

        /*
            dot_product_(nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread, 
                    iptr0, iptr, res_inout);

            #pragma omp barrier
            E_Float rho2=0.;
            // valid pour une zone
            for (E_Int ith = 0; ith < Nbre_thread_actif; ith++){rho2=rho2+res_inout[ith];}

            E_Float beta2 = rho2 * alpha1/(rho1*w1);
            printf("rho2 beta2: %g %g \n", rho2, beta2);

            // stockahe de p(j)- w(j)v(j) dans p(j)
            E_Float alpha = -w1 ;
            E_Float beta  = 1;
            vec_add_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptp, iptv, alpha, beta);   //vec1(l)=vec1(l)*beta + alpha*vec2(l)

            //p(j+1) = r(j)+ beta(j+1)(p(j) -  w(j)v(j)
            alpha = beta2 ;
            beta  = 1;
            vec_add_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptp, iptr, alpha, beta);
         */

            // 1.) 
            cl_vec_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptp);
            //v(j+1) = (A D_inv) p(j) 
            matvec_(nd, ithread_loc, param_int[nd], 
                    ipt_ind_dm_thread,
                    ipti[nd], iptj[nd], iptk[nd], iptvol[nd],
                    iptp , iptv);

            // 2.)
            // r(j=0).v(j+1)
            res_inout[ithread_loc-1]=0.;
            dot_product_(nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptv, iptr0,  res_inout);
            #pragma omp barrier
            E_Float tmp=0.;
            // valid pour une zone
            for (E_Int ith = 0; ith < Nbre_thread_actif; ith++){tmp=tmp+res_inout[ith];}
            //E_Float alpha2 = rho2/tmp;
            E_Float alpha2 = rho1/tmp;

            // 3.)   h = y(j-1) + alpha2 p(j-1)
            E_Float alpha = alpha2;
            E_Float beta= 1;
            vec_add_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptsol, iptp, alpha, beta);

            // 4.)   s = r(j-1) - alpha2 v(j)  (s stockee dans r)
            alpha = -alpha2;
            beta= 1;
            vec_add_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptr, iptv, beta, alpha);

            // 6;) t= As
            cl_vec_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptr);
            //t(j+1) = (A D_inv) s(j+1)  (t dans v, s dans r)
            matvec_(nd, ithread_loc, param_int[nd],
                    ipt_ind_dm_thread, 
                    ipti[nd], iptj[nd], iptk[nd], iptvol[nd],
                    iptr , iptv);

            //7.) w = t(j+1).s(j+1)/t(j+1).t(j+1)
            #pragma omp barrier
            // t(j+1).t(j+1)
            res_inout[ithread_loc-1]=0.;
            dot_product_(nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptv, iptv,  res_inout);
            #pragma omp barrier
            tmp=0.;
            for (E_Int ith = 0; ith < Nbre_thread_actif; ith++)
              {tmp=tmp+res_inout[ith];
              }

            //t(j+1).s(j+1)
            #pragma omp barrier
            res_inout[ithread_loc-1]=0.;
            dot_product_(nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread, iptr, iptv,  res_inout);
            #pragma omp barrier
            E_Float tmp2=0.;
            for (E_Int ith = 0; ith < Nbre_thread_actif; ith++){tmp2=tmp2+res_inout[ith];}

            E_Float w2 = tmp2/tmp;

            //8.) y(j)= h(j) + w s(j)
            alpha = w2;
            beta= 1;
            vec_add_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptsol, iptr, alpha, beta);

            //9.) r(j)= s(j) - w t(j)
            alpha = -w2;
            beta= 1;
            vec_add_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptr, iptv, alpha, beta);

            alpha1=alpha2;

            //11.) rho2 = rO.r
            #pragma omp barrier
            res_inout[ithread_loc-1]=0.;
            dot_product_(nd, ithread_loc, param_int[nd], ipt_ind_dm_thread , iptr, iptr0,  res_inout);
            #pragma omp barrier
            E_Float rho2 =0.;
            for (E_Int ith = 0; ith < Nbre_thread_actif; ith++){rho2=rho2+res_inout[ith];}

            //12.) rho2 = rO.r
            E_Float beta2 = rho2/rho1*alpha2/w2;
            rho1  = rho2;
 
            //13.)  stockahe de p(j)- w(j)v(j) dans p(j)
            alpha = -w2 ;
            beta  = 1;
            vec_add_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptp, iptv, alpha, beta);   //vec1(l)=vec1(l)*beta + alpha*vec2(l)

            //p(j+1) = r(j)+ beta(j+1)(p(j) -  w(j)v(j)
            alpha = 1 ;
            beta  = beta2;
            vec_add_( nd, ithread_loc, param_int[nd] , ipt_ind_dm_thread , iptp, iptr, alpha, beta);
            
            #pragma omp barrier
            res_inout[ithread_loc-1]=0.;
            dot_product_(nd, ithread_loc, param_int[nd], ipt_ind_dm_thread , iptr, iptr,  res_inout);


       } //task
   #pragma omp barrier
   #pragma omp single
     {
      E_Float res=0.;
      //E_Float* res_inout = residu +nd*Nbre_thread_actif;
      /// valid pour une zone
      /// valid pour une zone
       /// valid pour une zone
      E_Float* res_inout = residu;
      for (E_Int ith = 0; ith < Nbre_thread_actif; ith++){res=res+res_inout[ith];}
      printf("residu Jac: %g %d \n", res, it_jacobi);
     }
 }// it jacobi 

} // Fin zone // omp


}
