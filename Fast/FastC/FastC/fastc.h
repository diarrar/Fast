# ifndef _FASTC_FASTC_H_
# define _FASTC_FASTC_H_

# define SIZECF(type,meshtype,sizecf){                                  \
    if      (type == 1                  ) sizecf=1;                     \
    else if (type == 2   && meshtype == 1) sizecf=8;                    \
    else if (type == -2  && meshtype == 1) sizecf=8;                    \
    else if (type == 44  && meshtype == 1) sizecf=12;                   \
    else if (type == -44 && meshtype == 1) sizecf=12;                   \
    else if (type == 22  && meshtype == 1) sizecf=4;                    \
    else if (type == -22 && meshtype == 1) sizecf=4;                    \
    else if (type == 3   && meshtype == 1) sizecf=9;                    \
    else if (type == 4   && meshtype == 1) sizecf=8;                    \
    else if (type == 4   && meshtype == 2) sizecf=4;                    \
    else if (type == 5   && meshtype == 1) sizecf=15;                   \
    else sizecf=-1;                                                     \
  }

#ifdef _MPI
#include "CMP/include/pending_message_container.hpp"
#include "CMP/include/recv_buffer.hpp"
#include "CMP/include/send_buffer.hpp"
#include "TRANSFERT/setInterpTransfersD.h"
typedef typename CMP::PendingMsgContainer<CMP::RecvBuffer> RecvQueue;
typedef typename CMP::PendingMsgContainer<CMP::SendBuffer> SendQueue;
#endif

# include "kcore.h"
# include "FastC/Fortran.h"

using namespace K_FLD;

namespace K_FASTC
{ 

  // Check a value in Numerics Dictionary
  E_Int checkNumericsValue(PyObject* numerics, const char* value,
                           E_Int& retInt, E_Float& retFloat, char*& retChar);


  PyObject* _motionlaw(PyObject* self, PyObject* args);
  PyObject* PygetRange(PyObject* self, PyObject* args);
  PyObject* souszones_list(PyObject* self, PyObject* args);
  PyObject* distributeThreads(PyObject* self, PyObject* args);
  PyObject* _updateNatureForIBMGhost(PyObject* self, PyObject* args);

  PyObject* init_metric(PyObject* self, PyObject* args);

  //===========
  // - init Numa (optimisation placement memoire sur DRAM)
  //===========
  PyObject* initNuma(PyObject* self, PyObject* args);
  PyObject* initNumaTransfer(PyObject* self, PyObject* args);

  //===== Distrib OMP Test

 void distributeThreads_c( E_Int**& ipt_param_int,  E_Float**& ipt_param_real, E_Int**& ipt_ind_dm,
                           E_Int& nidom          ,  E_Int* ipt_dtloc         , E_Int& mx_sszone   , E_Int& nstep, E_Int& nitrun, E_Int& display);

 E_Int topo_test( E_Int* topo, E_Int* nijk, E_Int& cells_tg, E_Int& lmin, E_Int& dim_i,  E_Int& dim_j, E_Int& dim_k);

  //=======
  // - BC -
  //=======
  /* For BCcompact */

  E_Int getRange( E_Int* ind_cgns,  E_Int* ind_fen, E_Int* param_int);
  E_Int getDir(E_Int* ijkv,  E_Int* ind_fen);


  void souszones_list_c( E_Int**& ipt_param_int, E_Float**& ipt_param_real, E_Int**& ipt_ind_dm, E_Int**& ipt_it_lu_ssdom,
                         E_Int* dtloc          , E_Int* ipt_iskip_lu      , E_Int lssiter_loc       , E_Int nidom    , 
                         E_Int nitrun          , E_Int nstep              , E_Int flag_res          , E_Int& lexit_lu, E_Int& lssiter_verif);

  //==============================
  // - Transfer with Python -
  //==============================a
  PyObject* __setInterpTransfersD(PyObject* self, PyObject* args);
  PyObject* ___setInterpTransfers(PyObject* self, PyObject* args);
  PyObject* ___setInterpTransfers4GradP(PyObject* self, PyObject* args);
 

  E_Int setIBCTransfersCommonVar2(E_Int bctype,
                                  E_Int* rcvPtsI, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin, E_Int& ithread,
                                  E_Float* xPC, E_Float* yPC, E_Float* zPC,
                                  E_Float* xPW, E_Float* yPW, E_Float* zPW,
                                  E_Float* xPI, E_Float* yPI, E_Float* zPI,
                                  E_Float* densPtr,
                                  E_Float* tmp, E_Int&  size,  E_Int& nvars,
                                  E_Float* param_real,
                                  E_Float** WIn, E_Float** WOut,
                                  //std::vector<E_Float*>& WIn, std::vector<E_Float*>& WOut,
                                  E_Int nbptslinelets=0, E_Float* linelets=NULL, E_Int* indexlinelets=NULL);

 //==============================
  // - Transfer with CMP library -
  //==============================

  /* Call to transfers from FastS */

  void setInterpTransfersFast(
  E_Float**& iptro_tmp    , E_Int& vartype             , E_Int*& param_int_tc, E_Float*& param_real_tc , E_Int**& param_int      , E_Float**& param_real      , E_Int*& ipt_omp,
  E_Int*& ipt_linelets_int, E_Float*& ipt_linelets_real, E_Int& it_target    , E_Int& nidom            , E_Float*& ipt_timecount, E_Int& mpi           ,
  E_Int& nitcfg           , E_Int& nssiter             , E_Int& rk           , E_Int& exploc           , E_Int& numpassage, E_Int& Nopass );
  
  /* Transferts FastS Intra process */
  void setInterpTransfersIntra(E_Float**& ipt_ro, E_Int& vartype         , E_Int*& ipt_param_int   , E_Float*& ipt_param_real, 
                              E_Int**& param_int, E_Float**& param_real  , E_Int*& ipt_omp         , E_Int*& ipt_linelets_int, E_Float*& ipt_linelets_real, E_Int& nitrun, E_Int& nidom,
                              E_Int& NoTransfert, E_Float*& ipt_timecount,
                              E_Int& nitcfg     , E_Int& nssiter         , E_Int& rk, E_Int& exploc, E_Int& numpassage );

  /*---
  //LBM
  ---*/

  /* Transfers FastLBM */
  void setInterpTransfersFastLBM(E_Float**& iptro_tmp       , E_Int& vartype,
                E_Int*& param_int_tc       , E_Float*& param_real_tc    ,
                E_Int**& param_int         , E_Float**& param_real      ,
                E_Int*& ipt_linelets_int   , E_Float*& ipt_linelets_real,
                E_Int& it_target           , E_Int& nidom               , 
                E_Float*& ipt_timecount    , E_Int& mpi                 ,
                E_Int& nitcfg              , E_Int& nssiter             ,
                E_Int& rk                  , E_Int& exploc              ,
                E_Int& numpassage          ,
                E_Float**& ipt_macros_local, E_Float**& ipt_Qneq_local  );
  
  /* Intra process */
  void setInterpTransfersIntraLBM(E_Float**& ipt_ro          , E_Int& vartype,
                E_Int*& ipt_param_int      , E_Float*& ipt_param_real   ,
                E_Int**& param_int         , E_Float**& param_real      ,
                E_Int*& ipt_linelets_int   , E_Float*& ipt_linelets_real,
                E_Int& TypeTransfert       , E_Int& nitrun              ,
                E_Int& nidom               , E_Int& NoTransfert         ,
                E_Float*& ipt_timecount    , E_Int& nitcfg              ,
                E_Int& nssiter             , E_Int& rk                  ,
                E_Int& exploc              , E_Int& numpassage          ,
                E_Float**& ipt_macros_local, E_Float**& ipt_Qneq_local  );

  #ifdef _MPI
  /* Transferts FastS Inter process */
  void setInterpTransfersInter(E_Float**& ipt_ro , E_Int& vartype        , E_Int*& ipt_param_int   , E_Float*& ipt_param_real ,
                               E_Int**& param_int, E_Float**& param_real , E_Int*& ipt_omp, E_Int*& ipt_linelets_int, E_Float*& ipt_linelets_real, E_Int& TypeTransfert, E_Int& nitrun, E_Int& nidom,
                               E_Int& NoTransfert,
                               std::pair<RecvQueue*, SendQueue*>*& pair_of_queue, E_Int& etiquette,
                               E_Float*& ipt_timecount,
                               E_Int& nitcfg, E_Int& nssiter, E_Int& rk, E_Int& exploc, E_Int& numpassage , E_Int& nb_send_buffer);

  /* Get Transfert Inter process */
  void getTransfersInter(E_Int& nbcom, E_Float**& ipt_roD, E_Int**& param_int, E_Float**& param_real, E_Int*& param_int_tc , std::pair<RecvQueue*, SendQueue*>*& pair_of_queue, E_Float*& ipt_timecount);

  /* Init Transfert Inter process */
  void init_TransferInter(std::pair<RecvQueue*, SendQueue*>*& pair_of_queue);

  /* Delete Transfert Inter process */
  void del_TransferInter(std::pair<RecvQueue*, SendQueue*>*& pair_of_queue);
  #endif

  /* Fonctions pour les lois de paroi */
  E_Float logf(E_Float x, E_Float a, E_Float b, E_Float c, E_Float d);
  E_Float logfprime(E_Float x, E_Float a, E_Float d) ;
  E_Float musker(E_Float x, E_Float a, E_Float b);
  E_Float muskerprime(E_Float x, E_Float a, E_Float b);
  E_Float fnutilde(E_Float nutilde, E_Float nut, E_Float rho, E_Float xmu);
  E_Float fnutildeprime(E_Float nutilde, E_Float nut, E_Float rho, E_Float xmu);


}
#endif
