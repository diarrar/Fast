/*    
    Copyright 2013-2025 Onera.

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
// Mettre a 1 pour un CPU timer
#define TIMER 0

# include "FastC/fastc.h"
# include "FastS/fastS.h"
# include "FastS/param_solver.h"
# include <string.h>
#include <cstdlib>
//# include <omp.h>
#if TIMER == 1
# include <ctime>
E_Float timein;
E_Float timeout;
#endif
using namespace std;
using namespace K_FLD;
#ifdef _MPI
#include <mpi.h>
#endif

#include <iostream>

//=============================================================================
// Compute pour l'interface pyTree
//=============================================================================
PyObject* K_FASTS::_computePT_laplacien(PyObject* self, PyObject* args)
{
  PyObject* zones; PyObject* metrics; PyObject* work; 

#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOO" , &zones , &metrics, &work)) return NULL;
#else 
  if (!PyArg_ParseTuple(args, "OOO" , &zones , &metrics, &work)) return NULL;
#endif

#if TIMER == 1
#ifdef _OPENMP
  clock_t c_start,c_end;
  c_start = clock();
  timein = omp_get_wtime();
#endif
#endif
  
  //* tableau pour stocker dimension sous-domaine omp *//
  E_Int threadmax_sdm  = __NUMTHREADS__;

  PyObject* tmp; 


  PyObject* dtlocArray  = PyDict_GetItemString(work,"dtloc"); FldArrayI* dtloc;
  K_NUMPY::getFromNumpyArray(dtlocArray, dtloc); E_Int* iptdtloc  = dtloc->begin();
  E_Int nssiter = iptdtloc[0];
  E_Int shift_omp= iptdtloc[11];
  E_Int* ipt_omp = iptdtloc + shift_omp;

  E_Int lcfl= 0;
 // 
 // 
 //
 // partie code specifique au Transfer Data
 // 
 // 
 //
 PyObject* pyParam_int_tc; PyObject* pyParam_real_tc; PyObject* iskipArray;
 PyObject* pyLinlets_int; PyObject* pyLinlets_real; 
 FldArrayI* param_int_tc; FldArrayI* iskip_lu; FldArrayI* linelets_int;
 FldArrayF* param_real_tc; FldArrayF* linelets_real;
 E_Int* ipt_iskip_lu; E_Int* ipt_param_int_tc;
 E_Float* ipt_param_real_tc; E_Float* ipt_linelets_real; E_Int* ipt_linelets_int;

 E_Int lssiter_loc;
 E_Int it_target = iptdtloc[4];

  ipt_linelets_int  = NULL;
  ipt_linelets_real = NULL; 
  ipt_param_int_tc  = NULL;
  ipt_param_real_tc = NULL;

// Fin partie code specifique au Transfer Data

  E_Int nidom    = PyList_Size(zones);

  FldArrayI n0_flt(nidom); E_Int* ipt_n0_flt = n0_flt.begin();


  E_Int**   ipt_param_int;  E_Int** ipt_ind_dm; E_Int** ipt_it_lu_ssdom; E_Int** ipt_degen;
  
  E_Float** ipt_param_real;
  E_Float** iptx;          E_Float** ipty;         E_Float** iptz;    E_Float** iptro; E_Float** iptro_p1;
  E_Float** ipti;          E_Float** iptj;         E_Float** iptk;    E_Float** iptvol;
  E_Float** ipti0;         E_Float** iptj0;        E_Float** iptk0;   E_Float** iptsrc;
  E_Float** ipti_df;       E_Float** iptj_df;      E_Float** iptk_df ; E_Float** iptvol_df;
  E_Float** iptventi;      E_Float** iptventj;     E_Float** iptventk;
  E_Float** iptrdm;        E_Float** iptkrylov;    E_Float** iptkrylov_transfer;
  E_Float** iptCellN  ;    E_Float** iptCellN_IBC; E_Float** ipt_cfl_zones;

  ipt_param_int     = new E_Int*[nidom*4];
  ipt_ind_dm        = ipt_param_int   + nidom;
  ipt_it_lu_ssdom   = ipt_ind_dm      + nidom;
  ipt_degen         = ipt_it_lu_ssdom + nidom;

  iptx              = new E_Float*[nidom*36];
  ipty              = iptx               + nidom;
  iptz              = ipty               + nidom;
  iptro             = iptz               + nidom;
  iptro_p1          = iptro              + nidom;
  ipti              = iptro_p1           + nidom;
  iptj              = ipti               + nidom;
  iptk              = iptj               + nidom;
  iptvol            = iptk               + nidom;
  ipti0             = iptvol             + nidom;
  iptj0             = ipti0              + nidom;
  iptk0             = iptj0              + nidom;
  ipti_df           = iptk0              + nidom;
  iptj_df           = ipti_df            + nidom;
  iptk_df           = iptj_df            + nidom;
  iptvol_df         = iptk_df            + nidom;
  iptventi          = iptvol_df          + nidom;
  iptventj          = iptventi           + nidom;
  iptventk          = iptventj           + nidom;
  ipt_param_real    = iptventk           + nidom;
  iptsrc            = ipt_param_real     + nidom;

  vector<PyArrayObject*> hook;
    
  E_Int nb_pulse      = 0;                  // =1 => pulse acoustic centree a l'origine
  E_Int ndimdx_max    = 0;
  E_Int ndimdx_sa     = 0;
  E_Int nisdom_lu_max = 0;
  E_Int ndimt         = 0;
  E_Int ndim_drodm    = 0;
  E_Int ndimt_flt     = 0;
  E_Int kles          = 0;
  E_Int kimpli        = 0;
  E_Int kfiltering    = 0;
  E_Int neq_max       = 5;
  E_Int iorder_flt    =10;


  for (E_Int nd = 0; nd < nidom; nd++)
  {
    // check zone
    PyObject* zone = PyList_GetItem(zones, nd); // domaine i

    /* Get numerics from zone */
    PyObject*   numerics  = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
    PyObject*          o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
    ipt_param_int[nd]     = K_PYTREE::getValueAI(o, hook);
                       o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_real"); 
    ipt_param_real[nd]    = K_PYTREE::getValueAF(o, hook);

    PyObject* metric = PyList_GetItem(metrics, nd); // metric du domaine i

    if( ipt_param_int[nd][ NDIMDX ] > ndimdx_max ) ndimdx_max = ipt_param_int[nd][ NDIMDX ];
    
    //
    //
    //Pointeur maillage
    //
    //
    if(ipt_param_int[nd][ LALE ]== 0 || ipt_param_int[nd][ LALE ]== 2 || ipt_param_int[nd][ LALE ]== 3)
    { GET_XYZ( "GridCoordinates", zone, iptx[nd], ipty[nd], iptz[nd])}
    else { GET_XYZ( "GridCoordinates#Init", zone, iptx[nd], ipty[nd], iptz[nd])}

    //
    //Pointeur var primitive + visco + distance paroi + cellN
    //
    //
    PyObject* sol_center; PyObject* t;
    sol_center   = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
    t            = K_PYTREE::getNodeFromName1(sol_center, "Density");
    iptro[nd]    = K_PYTREE::getValueAF(t, hook);
    t            = K_PYTREE::getNodeFromName1(sol_center, "Density_P1");
    iptro_p1[nd] = K_PYTREE::getValueAF(t, hook);
    t            = K_PYTREE::getNodeFromName1(sol_center, "Density_src");
    if (t == NULL) iptsrc[nd] = NULL;
    else iptsrc[nd]   = K_PYTREE::getValueAF(t, hook);

    //
    //
    // Pointeur metric
    //
    E_Float* dummy;
    GET_TI(METRIC_TI  , metric, ipt_param_int[nd], ipti [nd]  , iptj [nd]  , iptk [nd]  , iptvol[nd]   )
    GET_TI(METRIC_TI0 , metric, ipt_param_int[nd], ipti0[nd]  , iptj0[nd]  , iptk0[nd]  , dummy )
    GET_TI(METRIC_TIDF, metric, ipt_param_int[nd], ipti_df[nd], iptj_df[nd], iptk_df[nd], iptvol_df[nd])

    GET_VENT( METRIC_VENT, metric, ipt_param_int[nd], iptventi[nd], iptventj[nd], iptventk[nd] )

    ipt_ind_dm[ nd ]      =  K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_INDM) );
    ipt_it_lu_ssdom[ nd ] =  K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_ITLU) );

    ndimt      += ipt_param_int[nd][ NDIMDX ];
    ndim_drodm += ipt_param_int[nd][ NDIMDX ]*ipt_param_int[nd][ NEQ ];

  } // boucle zone
  
//  // Reservation tableau travail temporaire pour calcul du champ N+1

  //tableau pour jacobi
  FldArrayF  phi1(ndimt); E_Float* iptphi1 = phi1.begin();
  FldArrayF  phi2(ndimt); E_Float* iptphi2 = phi2.begin();    //optimisation memoire possible pour roflt2 (voir Funk)
  FldArrayF   b(ndimt*4); E_Float* iptb    = b.begin();    //optimisation memoire possible pour roflt2 (voir Funk)

  //printf("threadmax = %d  , mx_nidom= %d \n",threadmax_sdm, mx_nidom );
  FldArrayF  res(threadmax_sdm*nidom); E_Float* iptres = res.begin();    
  FldArrayI ijkv_sdm(         3*threadmax_sdm); E_Int* ipt_ijkv_sdm   =  ijkv_sdm.begin();
  FldArrayI topology(         3*threadmax_sdm); E_Int* ipt_topology   =  topology.begin();
  FldArrayI ind_CL(          18*threadmax_sdm); E_Int* ipt_ind_CL     =  ind_CL.begin();
  FldArrayI ind_dm_omp(      12*threadmax_sdm); E_Int* ipt_ind_dm_omp =  ind_dm_omp.begin();

  laplace(  ipt_param_int, ipt_param_real , nidom, nssiter, threadmax_sdm ,
            ipt_ijkv_sdm , ipt_ind_dm_omp , ipt_topology     , ipt_ind_CL        ,
            ipt_ind_dm   , iptdtloc       ,
            iptx         , ipty           , iptz             ,
            iptro        , iptro_p1       , 
            iptphi1      , iptphi2        , iptb,
            ipti         , iptj           , iptk             , iptvol        , 
            iptsrc, iptres);

  delete [] iptx; delete [] ipt_param_int;
  
  RELEASESHAREDN( dtlocArray     , dtloc);
  RELEASEHOOK(hook)

  Py_INCREF(Py_None);
  return Py_None;
 
}
