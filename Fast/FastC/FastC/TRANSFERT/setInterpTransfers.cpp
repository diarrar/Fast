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

# include "fastc.h"
# include "param_solver.h"

#include <stdlib.h>

using namespace std;
using namespace K_FLD;

//=============================================================================
// Idem: in place + from zone + tc compact au niveau base. Valid pour FastS uniquememnt
//=============================================================================
PyObject* K_FASTC::___setInterpTransfers(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zonesD;
  PyObject *pyVariables, *pydtloc;
  PyObject *pyParam_int, *pyParam_real;
  E_Int vartype, no_transfert, It_target, nstep, nitmax;
  E_Int rk, exploc, num_passage;
  
  if (!PYPARSETUPLE_(args, OOOO_ OO_ IIII_ IIII_,  
                    &zonesR, &zonesD, &pyVariables, &pydtloc, &pyParam_int, &pyParam_real,\
                    &It_target, &vartype, &no_transfert, &nstep,\
                    &nitmax, &rk, &exploc, &num_passage))
  {
    return NULL;
  }
  E_Int it_target=  E_Int(It_target);
  /* varType :
     1  : conservatives,
     11 : conservatives + ronutildeSA
     2  : (ro,u,v,w,t)
     21 : (ro,u,v,w,t) + ronutildeSA
     22 : (ro,u,v,w,t) + ronutildeSA + (gradxRo, gradyRo, gradzRo) + (gradxT, gradyT, gradzT)
     3  : (ro,u,v,w,p)
     31 : (ro,u,v,w,p) + ronutildeSA
     4  : (ro,u,v,w,t) + (Q1,..., QN) LBM
     41 : (ro,u,v,w,t) + (Sxx,...) + (corr_xx,...) + (Q1,...,QN) LBM OVERSET
     ---------------------------------------------
     5  : (ro,u,v,w,t) + (Sxx,...) Couplage NS LBM
     51 : Couplage NS LBM improved ( a coder ) */

  E_Int varType     = E_Int(vartype);
  E_Int NoTransfert = E_Int(no_transfert);

  E_Int kmd, cnNfldD, nvars, meshtype;

  if     ( vartype <= 3 &&  vartype >= 1) nvars =5;
  else if( vartype == 4 ) nvars =32; // LBM transfer, 19 or 27 Qs and 5 macros (32 max in total)
                                     // on majore pour la LBM, car nvar sert uniquememnt a dimensionner taille vector
  else if (vartype == 41) nvars =44; //38;    // LBM Overset : 19 or 27 Q, 5 macro and 6 gradients
  else if( vartype == 5 ) nvars =30;    // Hybrid NSLBM transfer, 5 macros, 19 Qs + 6 gradients
  else                    nvars =6;

  E_Int nidomR   = PyList_Size(zonesR);
  E_Int nidomD   = PyList_Size(zonesD);

  //pointeur pour stocker solution au centre ET au noeud
  E_Int* ipt_ndimdxD; E_Int** ipt_param_intR; E_Int** ipt_cnd;
  E_Float** ipt_roR; E_Float** ipt_roD; E_Float** ipt_roR_vert;  E_Float** ipt_roD_vert; E_Float** ipt_param_realR;
  E_Float** ipt_roR_Pnt2;
  E_Float** ipt_roD_Pnt2;
  
  E_Float** ipt_qR; E_Float** ipt_qD; E_Float** ipt_qR_vert; E_Float** ipt_qD_vert;
  E_Float** ipt_SR; E_Float** ipt_SD; E_Float** ipt_SR_vert; E_Float** ipt_SD_vert;
  E_Float** ipt_psiGR; E_Float** ipt_psiGD; E_Float** ipt_psiGR_vert; E_Float** ipt_psiGD_vert;

  //ipt_ndimdxR      = new E_Int*[nidomR*3];   // on stocke ndimdx  en centre et vertexe
  ipt_param_intR   = new E_Int*[nidomR];

  ipt_roR            = new E_Float*[nidomR*10]; //1
  ipt_roR_vert       = ipt_roR         + nidomR;//2
  ipt_param_realR    = ipt_roR_vert    + nidomR;//3
  ipt_roR_Pnt2       = ipt_param_realR + nidomR;//4
  ipt_qR             = ipt_roR_Pnt2    + nidomR;//5
  ipt_qR_vert        = ipt_qR          + nidomR;//6
  ipt_SR             = ipt_qR_vert     + nidomR;//7
  ipt_SR_vert        = ipt_SR          + nidomR;//8
  ipt_psiGR          = ipt_SR_vert     + nidomR;//9
  ipt_psiGR_vert     = ipt_psiGR       + nidomR;//10

  ipt_ndimdxD        = new E_Int[nidomD*8];  //on stocke ndimdx, imd, jmd, en centre et vertexe, meshtype et cnDfld
  ipt_cnd            = new E_Int*[nidomD];

  ipt_roD            = new E_Float*[nidomD*9];  //1
  ipt_qD             = ipt_roD         + nidomD;//2
  ipt_SD             = ipt_qD          + nidomD;//3
  ipt_psiGD          = ipt_SD          + nidomD;//4
  ipt_roD_vert       = ipt_psiGD       + nidomD;//5
  ipt_qD_vert        = ipt_roD_vert    + nidomD;//6
  ipt_SD_vert        = ipt_qD_vert     + nidomD;//7
  ipt_psiGD_vert     = ipt_SD_vert     + nidomD;//8
  ipt_roD_Pnt2       = ipt_psiGD_vert  + nidomD;//9

  vector<PyArrayObject*> hook;

  FldArrayI* dtloc;
  E_Int res_donor = K_NUMPY::getFromNumpyArray(pydtloc, dtloc);
  E_Int* iptdtloc = dtloc->begin();
  /*-------------------------------------*/
  /* Extraction tableau int et real de tc*/
  /*-------------------------------------*/
  FldArrayI* param_int;
  res_donor = K_NUMPY::getFromNumpyArray(pyParam_int, param_int);
  E_Int* ipt_param_int = param_int->begin();
  FldArrayF* param_real;
  res_donor = K_NUMPY::getFromNumpyArray(pyParam_real, param_real);
  E_Float* ipt_param_real = param_real->begin();

  /*------------------------------------------------*/
  /* RECUPERATION DU NOM DES VARIABLES A TRANSFERER */
  /*------------------------------------------------*/
  char* varname = NULL; char* varname1 = NULL; char* varname2 = NULL; char* varname3 = NULL;
  char* vartmp  = NULL; 
  E_Int nbvar_inlist   = PyList_Size(pyVariables);

  for (E_Int ivar = 0; ivar < nbvar_inlist; ivar++)
  {
    PyObject* tpl0= PyList_GetItem(pyVariables, ivar);
    if (PyString_Check(tpl0)) vartmp = PyString_AsString(tpl0);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(tpl0)) vartmp = (char*)PyUnicode_AsUTF8(tpl0);
#endif
    // printf("varname %s \n", vartmp);
    if     (ivar==0) {varname  = vartmp;}
    else if(ivar==1) {varname1 = vartmp;}  //En LBM, on a besoin d'echanger les macro !!ET!! les Q donc deux varname
    else if(ivar==2) {varname2 = vartmp;}  //Pour le couplage NS-LBM, on a besoin d'echanger les Q !!ET!! les macro !!ET!! les gradients
    else if(ivar==3) {varname3 = vartmp;}
    else {printf("Warning: souci varname setInterpTransfers \n"); }
  }


  //on recupere sol et solcenter ainsi que connectivite et taille zones Donneuses (tc)
  for (E_Int nd = 0; nd < nidomD; nd++)
  {
    PyObject* zoneD = PyList_GetItem(zonesD, nd);
#    include "FastC/TRANSFERT/getfromzoneDcompact_all.h"
  }
  //on recupere sol et solcenter taille zones receuveuses, param_int et param_real (t)
  for (E_Int nd = 0; nd < nidomR; nd++)
  {
    PyObject* zoneR = PyList_GetItem(zonesR, nd);
#     include "FastC/TRANSFERT/getfromzoneRcompact_all.h"
  }

  E_Int nbcomID     = ipt_param_int[2];
  E_Int shift_graph = nbcomID + 2;

  E_Int threadmax_sdm  = __NUMTHREADS__;
  E_Int ech            = ipt_param_int[ NoTransfert +shift_graph];
  E_Int nrac           = ipt_param_int[ ech +1 ];          //nb total de raccord
  E_Int nrac_inst      = ipt_param_int[ ech +2 ];          //nb total de raccord instationnaire
  E_Int timelevel      = ipt_param_int[ ech +3 ];          //nb de pas de temps stocker pour chaque raccord instationnaire
  E_Int nrac_steady    = nrac - nrac_inst;                 //nb total de raccord stationnaire

  //gestion nombre de pass pour raccord instationnaire
  E_Int pass_inst_deb=0;
  E_Int pass_inst_fin=1;
  if (nrac_inst > 0) pass_inst_fin=2;

  //
  //on optimise les transfert pour implicit local
  //
  E_Int impli_local[nidomR];
  E_Int nssiter  = iptdtloc[0];
  E_Int shift_omp= iptdtloc[11];
  E_Int* ipt_omp = iptdtloc + shift_omp;
  E_Int nbtask   = ipt_omp[nstep-1]; 
  E_Int ptiter   = ipt_omp[nssiter+ nstep-1];


  for (E_Int nd = 0; nd < nidomR; nd++) {impli_local[nd]=0;}//par defaut pas de transfert
  for (E_Int ntask = 0; ntask < nbtask; ntask++)            //transfert sur les zones modifiees � la ssiter nstep
  {
    E_Int pttask = ptiter + ntask*(6+threadmax_sdm*7);
    E_Int nd = ipt_omp[ pttask ];
    impli_local[nd] = 1;
  }
  E_Int maxlevel      =  iptdtloc[ 9];  //transfert sur les zones qui recupere leur valeur interpolees en LBM
  E_Int it_cycl_lbm   =  iptdtloc[10];
  E_Int level_it      =  iptdtloc[13+it_cycl_lbm];
  E_Int max_it        = pow(2, maxlevel-1);

  E_Int level_next_it =  maxlevel;
  if (it_cycl_lbm != max_it -1 ) { level_next_it = iptdtloc[13 +it_cycl_lbm +1];}

  for (E_Int nd = 0; nd < nidomR; nd++)   
    {
       if (  ipt_param_intR[nd][IFLOW] == 4)
       {
          E_Int level = ipt_param_intR[nd][LEVEL];
          if (level <=level_next_it +1 ){impli_local[nd]=1;}
          else                          {impli_local[nd]=0;}
       }
       //printf("implilocal %d %d %d \n", impli_local[nd], nd, it_cycl_lbm);
    }


  E_Int size_autorisation = nrac_steady+1;
  size_autorisation = K_FUNC::E_max(size_autorisation , nrac_inst+1);
  E_Int autorisation_transferts[pass_inst_fin][size_autorisation];

  E_Int ntab_int =18;
  E_Float cutoff_coef=1.e-12;


  // printf("nrac = %d, nrac_inst = %d, level= %d, it_target= %d , nitrun= %d \n",  nrac, nrac_inst, timelevel,it_target, NitRun);
  //on dimension tableau travail pour IBC
  E_Int nbRcvPts_mx =0; E_Int ibcTypeMax=0;
  for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
  {
    E_Int irac_deb= 0; E_Int irac_fin= nrac_steady;
    if (pass_inst == 1) { irac_deb = ipt_param_int[ ech + 4 + it_target ]; irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];}


    for  (E_Int irac=irac_deb; irac< irac_fin; irac++)
    {
      E_Int shift_rac =  ech + 4 + timelevel*2 + irac;
      if (ipt_param_int[ shift_rac+ nrac*10 ] > nbRcvPts_mx) nbRcvPts_mx = ipt_param_int[ shift_rac+ nrac*10 ];

      if (ipt_param_int[shift_rac+nrac*3] > ibcTypeMax) ibcTypeMax =  ipt_param_int[shift_rac+nrac*3];

      E_Int irac_auto= irac-irac_deb;
      autorisation_transferts[pass_inst][irac_auto]=0;

      if (exploc == 1)  //if(rk==3 && exploc == 2) // Si on est en explicit local, on va autoriser les transferts entre certaines zones seulement en fonction de la ss-ite courante
      {
        E_Int debut_rac = ech + 4 + timelevel*2 + nrac*ntab_int + 27*irac;
        E_Int levelD = ipt_param_int[debut_rac + 25];
        E_Int levelR = ipt_param_int[debut_rac + 24];
        E_Int cyclD  = nitmax/levelD;

        // Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse
        if (levelD > levelR && num_passage == 1)
        {
          if ( nstep%cyclD==cyclD-1 || (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==1) )
          {
            autorisation_transferts[pass_inst][irac_auto]=1;
          }
          else {continue;}
        }
        // Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
        else if (levelD < levelR && num_passage == 1)
        {
          if (nstep%cyclD==1 || nstep%cyclD==cyclD/4 || nstep%cyclD== cyclD/2-1 || nstep%cyclD== cyclD/2+1 || nstep%cyclD== cyclD/2+cyclD/4 || nstep%cyclD== cyclD-1)
            { autorisation_transferts[pass_inst][irac_auto]=1; }
          else {continue;}
        }
        // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse
        else if (levelD == levelR && num_passage == 1)
        {
          if (nstep%cyclD==cyclD/2-1 || (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==0) || nstep%cyclD==cyclD-1)
            { autorisation_transferts[pass_inst][irac_auto]=1; }
          else {continue;}
        }
        // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse (cas du deuxieme passage)
        else if (levelD == levelR && num_passage == 2)
        {
          if (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==1)
            { autorisation_transferts[pass_inst][irac_auto]=1; }
          else {continue;}
        }
        else {continue;}
      }
      // Sinon, on autorise les transferts  si la zone donneuse a ete modifiee a l'iteration nstep
      else 
      {
        E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ];
        if (impli_local[NoD]==1) autorisation_transferts[pass_inst][irac_auto]=1;
        //autorisation_transferts[pass_inst][irac_auto]=1;
      }
        
    }
    }


  E_Int size = (nbRcvPts_mx/threadmax_sdm)+1; // on prend du gras pour gerer le residus
  E_Int r =  size % 8;
  if (r != 0) size  = size + 8 - r;           // on rajoute du bas pour alignememnt 64bits
  if (ibcTypeMax <=1 ) size = 0;              // tableau inutile : SP voir avec Ivan

  FldArrayF  tmp(size*17*threadmax_sdm);
  E_Float* ipt_tmp = tmp.begin();
  
  E_Float** RcvFields = new E_Float*[ (nvars)*threadmax_sdm];
  E_Float** DnrFields = new E_Float*[              nvars*threadmax_sdm];

  //# pragma omp parallel default(shared)  num_threads(1)
# pragma omp parallel default(shared)
  {
      
    E_Float gamma, cv, muS, Cs, Ts, Pr;

#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num()+1;
    E_Int  Nbre_thread_actif = omp_get_num_threads(); // nombre de thread actif dans cette zone
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif

    E_Int indR, type;
    E_Int indD0, indD, i, j, k, ncfLoc/*, nocf*/, indCoef, noi, sizecoefs, /*Nbchunk,*/ imd, jmd, imdjmd;

    E_Float** vectOfRcvFields = RcvFields + nvars*(ithread-1);
    E_Float** vectOfDnrFields = DnrFields + nvars*(ithread-1);

    //1ere pass_inst: les raccords fixes
    //2eme pass_inst: les raccords instationnaires
    for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
    {
        //printf("pass_inst = %d, level= %d \n",  pass_inst, nrac_inst_level );
        E_Int irac_deb= 0; E_Int irac_fin= nrac_steady;
        if(pass_inst == 1)
        {
          irac_deb = ipt_param_int[ ech + 4 + it_target             ];
          irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];
        }

        for  (E_Int irac=irac_deb; irac< irac_fin; irac++)
        {
          E_Int irac_auto= irac-irac_deb;
          if (autorisation_transferts[pass_inst][irac_auto]==1)
          {
            E_Int shift_rac =  ech + 4 + timelevel*2 + irac;
            //printf("pass_inst= %d, irac=  %d, ithread= %d \n",pass_inst,irac , ithread );
            // ibc: -1 : raccord ID, 0 : wallslip, 1 slip, 2 : log , 3 : Musker, 4 : outpress
            E_Int ibcType =  ipt_param_int[shift_rac+nrac*3];
            E_Int ibc = 1;
            if (ibcType < 0) ibc = 0;

            E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5  ];
            E_Int loc      =  ipt_param_int[ shift_rac + nrac*9  ]; 
            E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 ];
            E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 ]; //neq fonction raccord rans/LES
            E_Int rotation =  ipt_param_int[ shift_rac + nrac*14 ]; //flag pour periodicite azimutale

            // COUPLAGE NS LBM - Recupere les solveurs des zones R et D
            E_Int solver_D=2; E_Int solver_R=2;
            if (nvars_loc == 11) {solver_R =4;}
            if (nvars_loc == -5) {solver_D =4; nvars_loc = 5;}
            if (nvars_loc == 19) {solver_D =4; solver_R=4;}
            if (nvars_loc == 27) {solver_D =4; solver_R=4;}

            E_Int overset  =  ipt_param_intR[NoD][LBM_OVERSET];        //flag pour overset en LBM
            if      (nvars_loc==19 && overset==0) nvars_loc = nvars_loc + 5;
            else if (nvars_loc==27 && overset==0) nvars_loc = nvars_loc + 5;
            else if (nvars_loc==19 && overset==1) nvars_loc = nvars_loc + 5 + 6 + 6;
            else if (nvars_loc==27 && overset==1) nvars_loc = nvars_loc + 5 + 6 + 6;
            // cout << nvars_loc << endl;
            // printf("irac=  %d, nvar/nvar_loc= %d %d,  ithread= %d  solverDR %d %d \n",irac , nvars,nvars_loc, ithread,  solver_D, solver_R);

            E_Int levelD   = ipt_param_intR[NoD][LEVEL];
            E_Int levelR   = ipt_param_intR[NoR][LEVEL];

            E_Int meshtype = ipt_ndimdxD[NoD + nidomD*6];
            E_Int cnNfldD  = ipt_ndimdxD[NoD + nidomD*7];
            E_Int* ptrcnd  = ipt_cnd[    NoD           ];

            if (loc == 0)
            {
              printf("Error: transferts optimises non code en vextex " SF_D3_ "\n", shift_rac + nrac*9  +1, NoD, NoR );
              //imd= ipt_ndimdxD[ NoD+ nidomD*4]; jmd= ipt_ndimdxD[ NoD + nidomD*5];
              imd = 0; jmd = 0;
            }
            else
            {
                /*--------------------------------------------------------------------*/
                /*                GESTION DES TRANSFERTS EN CENTER                    */
                /* 2 cas: - si pas de couplage NSLBM, transferts habituels,           */
                /*        - si couplage NSLBM, raccords a adapter                     */
                /*--------------------------------------------------------------------*/
                if (nvars_loc == 5 || nvars_loc == 6) // Transferts NS classiques ou LBM -> NS
                {
                 for (E_Int eq = 0; eq < nvars_loc; eq++)
                 {
                   vectOfRcvFields[eq] = ipt_roR[ NoR] + eq*ipt_param_intR[ NoR][ NDIMDX ];
                   vectOfDnrFields[eq] = ipt_roD[ NoD] + eq*ipt_param_intR[ NoD][ NDIMDX ];
                 }
                }
                else if (nvars_loc == 24 || nvars_loc == 32) // Transferts LBM classiques
                {
                  // On commence par copier les 5 variables macros
                  for (E_Int eq = 0; eq < 5; eq++)
                  {
                    vectOfRcvFields[eq] = ipt_roR[ NoR] + eq*ipt_param_intR[ NoR][ NDIMDX ];
                    vectOfDnrFields[eq] = ipt_roD[ NoD] + eq*ipt_param_intR[ NoD][ NDIMDX ];
                  }
                  // Puis on copie les fonctions de distribution
                  for (E_Int eq = 5; eq < nvars_loc; eq++)
                  {
                    vectOfRcvFields[eq] = ipt_qR[ NoR] + (eq-5)*ipt_param_intR[ NoR][ NDIMDX ];
                    vectOfDnrFields[eq] = ipt_qD[ NoD] + (eq-5)*ipt_param_intR[ NoD][ NDIMDX ];
                  }
                }
                else if (nvars_loc == 11 ) // //Transfert NS -> LBM    
                {
                  // On commence par copier les 5 variables macros
                  for (E_Int eq = 0; eq < 5; eq++)
                  {
                    vectOfRcvFields[eq] = ipt_roR[ NoR] + eq*ipt_param_intR[ NoR][ NDIMDX ];
                    vectOfDnrFields[eq] = ipt_roD[ NoD] + eq*ipt_param_intR[ NoD][ NDIMDX ];
                  }
                  // Puis on copie les gradients
                  for (E_Int eq = 5; eq < nvars_loc; eq++)
                  {
                    vectOfRcvFields[eq] = ipt_SR[ NoR] + (eq-5)*ipt_param_intR[ NoR][ NDIMDX ];
                    vectOfDnrFields[eq] = ipt_SD[ NoD] + (eq-5)*ipt_param_intR[ NoD][ NDIMDX ];
                  }
                }
                else if (nvars_loc == 36 || nvars_loc == 44) // //Transfert LBM  overset   
                {
                  // On commence par copier les 5 variables macros
                  for (E_Int eq = 0; eq < 5; eq++)
                  {
                    vectOfRcvFields[eq] = ipt_roR[ NoR] + eq*ipt_param_intR[ NoR][ NDIMDX ];
                    vectOfDnrFields[eq] = ipt_roD[ NoD] + eq*ipt_param_intR[ NoD][ NDIMDX ];
                  }
                  // Puis on copie les gradients
                  for (E_Int eq = 5; eq < 11; eq++)
                  {
                    vectOfRcvFields[eq] = ipt_SR[ NoR] + (eq-5)*ipt_param_intR[ NoR][ NDIMDX ];
                    vectOfDnrFields[eq] = ipt_SD[ NoD] + (eq-5)*ipt_param_intR[ NoD][ NDIMDX ];
                  }
                  for (E_Int eq =11; eq < 17; eq++)
                  {
                    vectOfRcvFields[eq] = ipt_psiGR[ NoR] + (eq-11)*ipt_param_intR[ NoR][ NDIMDX ];
                    vectOfDnrFields[eq] = ipt_psiGD[ NoD] + (eq-11)*ipt_param_intR[ NoD][ NDIMDX ];
                  }
                  for (E_Int eq =17; eq < nvars_loc; eq++)
                  {
                    vectOfRcvFields[eq] = ipt_qR[ NoR] + (eq-17)*ipt_param_intR[ NoR][ NDIMDX ];
                    vectOfDnrFields[eq] = ipt_qD[ NoD] + (eq-17)*ipt_param_intR[ NoD][ NDIMDX ];
                  }
                }
               
              imd= ipt_param_intR[ NoD ][ NIJK ]; jmd= ipt_param_intR[ NoD ][ NIJK+1];
            }

            imdjmd = imd*jmd;


            ////
            //  Interpolation parallele
            ////
            ////

            E_Int nbRcvPts = ipt_param_int[ shift_rac +  nrac*10 ];

            E_Int pos;
            pos  = ipt_param_int[ shift_rac + nrac*7 ]; E_Int* ntype      = ipt_param_int  + pos;
            pos  = pos +1 + ntype[0]                  ; E_Int* types      = ipt_param_int  + pos;
            pos  = ipt_param_int[ shift_rac + nrac*6 ]; E_Int* donorPts   = ipt_param_int  + pos;
            pos  = ipt_param_int[ shift_rac + nrac*12]; E_Int* rcvPts     = ipt_param_int  + pos;   // donor et receveur inverser car storage donor
            pos  = ipt_param_int[ shift_rac + nrac*8 ]; E_Float* ptrCoefs = ipt_param_real + pos;

            E_Int nbInterpD = ipt_param_int[ shift_rac +  nrac ]; E_Float* xPC=NULL; E_Float* xPI=NULL; E_Float* xPW=NULL; E_Float* densPtr=NULL;
            if (ibc == 1)
            {
              xPC     = ptrCoefs + nbInterpD;
              xPI     = ptrCoefs + nbInterpD +3*nbRcvPts;
              xPW     = ptrCoefs + nbInterpD +6*nbRcvPts;
              densPtr = ptrCoefs + nbInterpD +9*nbRcvPts;
            }

            E_Int ideb      = 0;
            E_Int ifin      = 0;
            E_Int shiftCoef = 0;
            E_Int shiftDonor  = 0;
            
            for (E_Int ndtyp = 0; ndtyp < ntype[0]; ndtyp++)
            {
              type      = types[ifin];

              SIZECF(type, meshtype, sizecoefs);
              ifin =  ifin + ntype[ 1 + ndtyp];

              E_Int pt_deb, pt_fin;

              // Calcul du nombre de champs a traiter par chaque thread
              E_Int size_bc =  ifin-ideb;
              E_Int chunk   =  size_bc/Nbre_thread_actif;
              E_Int r       =  size_bc - chunk*Nbre_thread_actif;
              // pts traitees par thread
              if (ithread <= r)
              { pt_deb = ideb + (ithread-1)*(chunk+1);           pt_fin = pt_deb + (chunk+1); }
              else { pt_deb = ideb + (chunk+1)*r+(ithread-r-1)*chunk; pt_fin = pt_deb + chunk; }

              //Si type 0, calcul sequentiel
              if  ( type == 0 )
              { if (ithread ==1 ){ pt_deb = ideb; pt_fin = ifin;}
                else             { pt_deb = ideb; pt_fin = ideb;}
              }

              noi       = shiftDonor;                             // compteur sur le tableau d indices donneur
              indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
    
              E_Int shiftv  =0;

              if ( (nvars_loc==5 || (ibc==1 && solver_R==4)) )
              {  
#             include "FastC/TRANSFERT/commonInterpTransfers_reorder_5eq.h"
              }
              else if (nvars_loc==6 )
              {
#             include "commonInterpTransfers_reorder_6eq.h"
              }
              else if(nvars_loc==19)
              {
#             include "FastC/TRANSFERT/commonInterpTransfers_reorder_19eq.h" 
              }
              else
              {
#             include "FastC/TRANSFERT/commonInterpTransfers_reorder_neq.h"
              }
            
              // =============================================================
              // RAFFINEMENT DE MAILLAGE LBM/LBM
              // =============================================================
              if (solver_D==4 && solver_R==4 && levelD > levelR)
              {
#               include "FastC/TRANSFERT/includeTransfers_coarse2fine_LBM.h"
              }
              else if (solver_D == 4 && solver_R == 4 && levelD < levelR)
              {
#               include "FastC/TRANSFERT/includeTransfers_fine2coarse_LBM.h"
              }
              // =============================================================

              // =============================================================
              // COUPLAGE NS-LBM: changement d'unite
              // =============================================================
              // Code que pour le cas D3Q19 pour le moment, a adapter au 27
              if (solver_D==4 && solver_R<4)
              {
               // Transfert LBM vers NS: repasse dans unites SI
#              include "FastC/TRANSFERT/includeTransfers_dimLBMtoNS.h"
              }
              else if (solver_D<4 && solver_R==4)
              {
               // Transfert NS vers LBM : adimensionnement
#              include "FastC/TRANSFERT/includeTransfers_dimNStoLBM.h"
              }

              // Prise en compte de la periodicite par rotation
              if (rotation == 1)
              {
               E_Float* angle = ptrCoefs + nbInterpD;
#              include "FastC/TRANSFERT/includeTransfers_rotation.h"
              }
            
              // ibc
              if (ibc == 1)
              {
                  K_FASTC::setIBCTransfersCommonVar2(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                  xPC, xPC+nbRcvPts, xPC+nbRcvPts*2,
                                  xPW, xPW+nbRcvPts, xPW+nbRcvPts*2,
                                  xPI, xPI+nbRcvPts, xPI+nbRcvPts*2,
                                  densPtr, 
                                  ipt_tmp, size, nvars,
                                  ipt_param_realR[ NoR ],
                                  vectOfDnrFields, vectOfRcvFields);
                  // Si calcul LBM => reconstruction a l'equilibre des Qs
                  if (solver_R==4)
                   {
#                   include "FastC/TRANSFERT/includeTransfers_LBM_feq.h"
                   }
              }//ibc
            
              ideb       =  ideb + ntype[ 1 + ndtyp];
              shiftCoef  = shiftCoef  +  ntype[1+ndtyp]*sizecoefs; //shift coef   entre 2 types successif
              shiftDonor = shiftDonor +  ntype[1+ndtyp];           //shift donor entre 2 types successif
            }// type
          }// autorisation transfert
          }//irac
    }//pass_inst
  }// omp


  delete [] ipt_param_intR; delete [] ipt_roR; delete [] ipt_ndimdxD; delete [] ipt_roD; delete [] ipt_cnd;
  delete [] RcvFields; delete [] DnrFields;

  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  RELEASESHAREDN(pydtloc        , dtloc        );
  RELEASESHAREDN(pyParam_int    , param_int    );
  RELEASESHAREDN(pyParam_real   , param_real   );

  Py_INCREF(Py_None);

  return Py_None;
}

