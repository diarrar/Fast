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
# include <math.h>
#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif
using namespace std;
using namespace K_FLD;

# include "TRANSFERT/IBC/commonLaws.h"
# define NUTILDE_FERRARI 2

//=============================================================================
//Retourne -2: incoherence entre meshtype et le type d'interpolation
//         -1: type invalide
//          1: ok
// Entree/Sortie:  (ro,u,v,w,t) ( + nutildeSA )
//=============================================================================
E_Int K_FASTC::setIBCTransfersCommonVar2(
					     E_Int bctype,
					     E_Int* rcvPts, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin, E_Int& ithread,
					     E_Float* xPC, E_Float* yPC, E_Float* zPC,
					     E_Float* xPW, E_Float* yPW, E_Float* zPW,
					     E_Float* xPI, E_Float* yPI, E_Float* zPI, 
					     E_Float* densPtr, 
					     E_Float* tmp, E_Int& size, E_Int& nvars, 
					     E_Float*  param_real,
					     E_Float** vectOfDnrFields, E_Float** vectOfRcvFields,
					     E_Int nbptslinelets, E_Float* linelets, E_Int* indexlinelets)
{
  E_Float Pr           = param_real[ PRANDT ];
  E_Float Ts           = param_real[ TEMP0 ];
  E_Float Tinf         = param_real[ TINF ];
  E_Float Pinf         = param_real[ PINF ];
  E_Float Roinf        = param_real[ ROINF ];
  E_Float Cs           = param_real[ CS ];
  E_Float muS          = param_real[ XMUL0 ];
  E_Float cv           = param_real[ CVINF ];
  E_Float gamma        = param_real[ GAMMA ];
  E_Float K_wire       = param_real[ KWire ];
  E_Float Delta_V_wire = param_real[ DeltaVWire ];
  E_Float Diam_wire    = param_real[ DiameterWire ];
  E_Float Ct_WM        = param_real[ CtWire ];
  E_Float R_gas        = Pinf/(Roinf*Tinf);


  E_Int bctypeLocal; 
  int   motionType      = (int) param_real[MotionType];
  //[AJ] Keep for now
  //E_Float transpeed[3]    = {param_real[TransSpeed],param_real[TransSpeed+1],param_real[TransSpeed+2]};
  //E_Float axispnt[3]      = {param_real[AxisPnt],param_real[AxisPnt+1],param_real[AxisPnt+2]};
  //E_Float axisvec[3]      = {param_real[AxisVec],param_real[AxisVec+1],param_real[AxisVec+2]};
  //E_Float omg             = param_real[OMG];  

  E_Float cmx,cmy,cmz;
  E_Float kvcmx,kvcmy,kvcmz;
  E_Float tmp_x,tmp_y,tmp_z;
  E_Float uGrid_local,vGrid_local,wGrid_local;
  E_Float normalVelGrid_local;
  E_Int c_ale = max(0,min(1,motionType));
    
  E_Float* pressPtr = densPtr + 1*nbRcvPts;
  E_Float* vxPtr    = densPtr + 2*nbRcvPts;
  E_Float* vyPtr    = densPtr + 3*nbRcvPts; 
  E_Float* vzPtr    = densPtr + 4*nbRcvPts;

  E_Float* utauPtr = NULL;
  E_Float* yplusPtr = NULL;
  E_Float* kcurvPtr = NULL;

  E_Float* d1 = NULL;
  E_Float* d2 = NULL;
  E_Float* d3 = NULL;
  E_Float* d4 = NULL;
  E_Float* d5 = NULL;

  E_Float* tempPtr      = NULL;
  E_Float* tempExtraPtr = NULL;

  E_Float* gradxPressPtr = NULL;
  E_Float* gradyPressPtr = NULL;
  E_Float* gradzPressPtr = NULL;

  E_Float* gradxUPtr = NULL;
  E_Float* gradyUPtr = NULL;
  E_Float* gradzUPtr = NULL;

  E_Float* gradxVPtr = NULL;
  E_Float* gradyVPtr = NULL;
  E_Float* gradzVPtr = NULL;

  E_Float* gradxWPtr = NULL;
  E_Float* gradyWPtr = NULL;
  E_Float* gradzWPtr = NULL;

  E_Float* motionPtr   = NULL;
  E_Float* transpeedPtrX = NULL;
  E_Float* transpeedPtrY = NULL;
  E_Float* transpeedPtrZ = NULL;
  E_Float* axispntPtrX = NULL;
  E_Float* axispntPtrY = NULL;
  E_Float* axispntPtrZ = NULL;
  E_Float* axisvecPtrX = NULL;
  E_Float* axisvecPtrY = NULL;
  E_Float* axisvecPtrZ = NULL;
  E_Float* omgPtr   = NULL;

  E_Float* y_linePtr = NULL;
  E_Float* u_linePtr = NULL;
  E_Float* nutilde_linePtr = NULL;
  E_Float* psi_linePtr = NULL;
  E_Float* matm_linePtr = NULL;
  E_Float* mat_linePtr = NULL;
  E_Float* matp_linePtr = NULL;
  E_Float* alphasbeta_linePtr = NULL;
  E_Float* index_linePtr = NULL;

  // bctype = 3 for all Musker, SA, MuskerLin, & SALin to avoid adding
  // more bctype conditions in if statements.
  // bctypeLocal will be kept for a flag switch for SA (32), MuskerLin (331), & SALin (332).
  // These are in development and will be added in the near future.
  bctypeLocal = bctype;
  if ( bctypeLocal == 32 || bctypeLocal == 331 || bctypeLocal == 332){
    bctype=3;
  }
  
  if (motionType==3){
    E_Int shift_var=0;
    // log, Musker, TBLE, MuskerMob, Pohlhausen, Thwaites - also have utau & yplus - need the shift
    if (bctype == 2 || bctype == 3 || bctype == 6 || bctype == 7 || bctype == 8 || bctype == 9) shift_var=2;
      
    motionPtr    =densPtr + (14+shift_var)*nbRcvPts;

    transpeedPtrX=densPtr + (15+shift_var)*nbRcvPts;
    transpeedPtrY=densPtr + (16+shift_var)*nbRcvPts;
    transpeedPtrZ=densPtr + (17+shift_var)*nbRcvPts;

    axispntPtrX=densPtr + (18+shift_var)*nbRcvPts;
    axispntPtrY=densPtr + (19+shift_var)*nbRcvPts;
    axispntPtrZ=densPtr + (20+shift_var)*nbRcvPts;

    axisvecPtrX=densPtr + (21+shift_var)*nbRcvPts;
    axisvecPtrY=densPtr + (22+shift_var)*nbRcvPts;
    axisvecPtrZ=densPtr + (23+shift_var)*nbRcvPts;

    omgPtr     =densPtr + (24+shift_var)*nbRcvPts;
  }

  if (bctype == 11)
    {
      nbptslinelets = param_real[ NBPTS_LINELETS ];
    }

  if (nbptslinelets == 0 && bctype == 11) //TBLE_FULL -> Musker
    {
      bctype = 3; //Musker
    }

  if ( bctype == 0 || bctype == 1 || bctype == 4 || bctype == 140 || bctype == 141)
    {;}
  else if (bctype==100)//slip + curvature radius
    {
      kcurvPtr = densPtr+5*nbRcvPts;
    }
  else if (bctype == 2 || bctype == 3 || bctype == 6 || bctype == 7 || bctype == 8 || bctype == 9)// log, Musker, TBLE, MuskerMob, Pohlhausen, Thwaites
    {
      utauPtr  = densPtr+5*nbRcvPts;
      yplusPtr = densPtr+6*nbRcvPts;
    }
  else if (bctype == 5)// injection
    {
      d1 = densPtr+5*nbRcvPts;
      d2 = densPtr+6*nbRcvPts;
      d3 = densPtr+7*nbRcvPts;
      d4 = densPtr+8*nbRcvPts;
      d5 = densPtr+9*nbRcvPts;
    }
  else if (bctype == 10) // Mafzal
    {
      utauPtr       = densPtr+5*nbRcvPts;
      yplusPtr      = densPtr+6*nbRcvPts;

      gradxPressPtr = densPtr+7*nbRcvPts;
      gradyPressPtr = densPtr+8*nbRcvPts;
      gradzPressPtr = densPtr+9*nbRcvPts;

      //E_Int   mafzalMode    = param_real[ MAFZAL_MODE ];
      //E_Float alphaGradP    = param_real[ ALPHAGRADP ];
      //nbptslinelets         = param_real[ NBPTS_LINELETS ];
      // std::cout << "mafzalMode = " << mafzalMode << " alpha = " << alphaGradP << " nbpts linelets = " << nbptslinelets << std::endl;
    }
  else if (bctype == 11) // TBLE-FULL
    {
      utauPtr       = densPtr+5*nbRcvPts;
      yplusPtr      = densPtr+6*nbRcvPts;

      gradxPressPtr = densPtr+7*nbRcvPts;
      gradyPressPtr = densPtr+8*nbRcvPts;
      gradzPressPtr = densPtr+9*nbRcvPts;

      gradxUPtr = densPtr+10*nbRcvPts;
      gradyUPtr = densPtr+11*nbRcvPts;
      gradzUPtr = densPtr+12*nbRcvPts;

      gradxVPtr = densPtr+13*nbRcvPts;
      gradyVPtr = densPtr+14*nbRcvPts;
      gradzVPtr = densPtr+15*nbRcvPts;

      gradxWPtr = densPtr+16*nbRcvPts;
      gradyWPtr = densPtr+17*nbRcvPts;
      gradzWPtr = densPtr+18*nbRcvPts;

      y_linePtr          = densPtr+nbRcvPts*(19+nbptslinelets*0);
      u_linePtr          = densPtr+nbRcvPts*(19+nbptslinelets*1);
      nutilde_linePtr    = densPtr+nbRcvPts*(19+nbptslinelets*2);
      psi_linePtr        = densPtr+nbRcvPts*(19+nbptslinelets*3);
      matm_linePtr       = densPtr+nbRcvPts*(19+nbptslinelets*4);
      mat_linePtr        = densPtr+nbRcvPts*(19+nbptslinelets*5);
      matp_linePtr       = densPtr+nbRcvPts*(19+nbptslinelets*6);
      alphasbeta_linePtr = densPtr+nbRcvPts*(19+nbptslinelets*7);
      index_linePtr      = densPtr+nbRcvPts*(19+nbptslinelets*7+1);
    }
  else if (bctype == 12 || bctype == 13 )// isothermal or heat flux
    {
      tempPtr      = densPtr+5*nbRcvPts;
      tempExtraPtr = densPtr+6*nbRcvPts;
    }
  else 
    {
      printf("Warning !!! setIBCTransfersCommonVar2: bcType " SF_D_ " not implemented.\n", bctype);
      return 0;
    }

  /* lois de paroi */
  E_Float roext, uext, pext, text, muext, yext, yplus, yibc, eta, delta;
  E_Float uext_wall,uext_image;
  E_Float uscaln_wall,uscaln_image;
  E_Float gradxPext, gradyPext, gradzPext;
  E_Float gradxUext, gradyUext, gradzUext;
  E_Float gradxVext, gradyVext, gradzVext;
  E_Float gradxWext, gradyWext, gradzWext;
  E_Float uscaln, un, vn, wn, ut, vt, wt, utau, utauv, utau0, umod;  
  E_Float utinf, vtinf, wtinf;
  E_Float t_thwaites, a_thwaites, b_thwaites, c_thwaites, lambda_thwaites, m_thwaites;
  E_Float ut_thwaites, vt_thwaites, wt_thwaites, utau_thwaites, x_T, alpha_T;
  E_Float ut_transition, vt_transition, wt_transition;
  E_Float aa, bb, dd, fp, tp, f1v;
  E_Float expy, denoml10, ax, px, kx, y1, y2, l1, l2, l3, l4, l5, l11, l12, l13;
  E_Float ag11, ag12, ag13, bg10, bg11, bg12, bg13, bg14, cg10, cg11, cg12, cg13, cg14;
  E_Float ag22, ag23, bg20, bg21, bg22, bg23, bg24, cg20, cg21, cg22, cg23, cg24;
  E_Float ucible0, ucible, vcible, wcible, tcible, nutilde, signibc, twall, rowall, muwall;
  E_Int npass;
  // Lois de paroi: criteres d'arret pour estimer le frottement par Newton
  E_Float newtoneps = 1.e-7; // critere d'arret pour u+
  E_Float newtonepsnutilde = 1.e-10; // critere d arret pour nutilde
  E_Float newtonepsprime = 1.e-12;// critere d'arret pour la derivee
  E_Float cvgam = cv*(gamma-1.);
  E_Float cvgaminv = 1./(cvgam);
  E_Float coefSuth = muS * (1.+Cs/Ts);
  E_Float Tsinv = 1./Ts;
  E_Float kappa = 0.4; // Constante de Von Karman
  E_Float kappainv = 1./kappa;
  E_Float cc = 5.2; //pour la loi log
  E_Float one_third = 1./3.;
  /* fin parametres loi de parois*/

  //E_Int nvars    = vectOfDnrFields.size();
  //E_Int nvarsRcv = vectOfRcvFields.size();

  E_Float a0,a1,a2,b0,b1,b2,n0,n1,n2;
  E_Float t0,t1,t2;
  E_Float normb, ro, u, v, w, t;
  E_Float vnc, alpha, beta, alphasbeta;
  E_Float* roOut = vectOfRcvFields[0];// ro
  E_Float* uOut  = vectOfRcvFields[1];// u
  E_Float* vOut  = vectOfRcvFields[2];// v
  E_Float* wOut  = vectOfRcvFields[3];// w
  E_Float* tOut  = vectOfRcvFields[4];// temperature


  E_Float* varSAOut = NULL;
  
  //---------------------------------
  // ODE-based wall model
  //---------------------------------
  E_Float Cv1cube = pow(7.1,3);
  E_Int nmax      = 20;
  
  E_Float L2norm ;
  E_Float L2norm0;

  E_Float ynm,ynp,dy,dym,dyp,nutm,nutp;
  E_Float xim,xi,xip,m;
  E_Float nutrm,nutrp,nutr;

  E_Float* ipt_u1dold      = NULL;
  E_Float* yline           = NULL;
  E_Float* u_line          = NULL;
  E_Float* nutilde_line    = NULL;
  E_Float* matm_line       = NULL;
  E_Float* mat_line        = NULL;
  E_Float* matp_line       = NULL;
  E_Float* alphasbeta_line = NULL;

  FldArrayF u1dold(nbptslinelets);
  ipt_u1dold = u1dold.begin();

  if (nbptslinelets != 0 && (bctype == 6))
    {
      yline           = linelets;
      u_line          = linelets + nbRcvPts*nbptslinelets;
      nutilde_line    = linelets + nbRcvPts*nbptslinelets*2;
      matm_line       = linelets + nbRcvPts*nbptslinelets*3;
      mat_line        = linelets + nbRcvPts*nbptslinelets*4;
      matp_line       = linelets + nbRcvPts*nbptslinelets*5;
      alphasbeta_line = linelets + nbRcvPts*nbptslinelets*6;
    }
  //---------------------------------

  if (nvars == 6) varSAOut = vectOfRcvFields[5]; // nutildeSA

  // if ( (bctype==2 || (bctype==3)) && nvars < 6)
  // {
  //   printf("Warning: setIBCTransfersCommonVar2: number of variables (<6) inconsistent with bctype.\n");
  //   return 0;
  // }

  if (bctype == 100)//wallslip + curvature radius
    {
#ifdef _OPENMP4
#pragma omp simd
#endif 
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
	  E_Int indR = rcvPts[noind+ideb];

	  // vitesse
	  u = uOut[indR];
	  v = vOut[indR];
	  w = wOut[indR];
# include "TRANSFERT/IBC/commonBCType0.h"
	  uOut[indR] = ucible; 
	  vOut[indR] = vcible; 
	  wOut[indR] = wcible;
	  if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

	  pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam; 
	  densPtr[ noind + ideb] = roOut[indR];

	  vxPtr[noind+ideb] = uOut[indR];
	  vyPtr[noind+ideb] = vOut[indR];
	  vzPtr[noind+ideb] = wOut[indR];
        }
    }
  else if (bctype == 0) // wallslip
    {
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
	  E_Int indR = rcvPts[noind+ideb];

	  // vitesse
	  u = uOut[indR];
	  v = vOut[indR];
	  w = wOut[indR];

	  //[AJ]
# include "TRANSFERT/IBC/commonIBCmotionAbs2Rel.h"  
	  
# include "TRANSFERT/IBC/commonBCType0.h"

	  uOut[indR] = ucible;
	  vOut[indR] = vcible;
	  wOut[indR] = wcible;

	  //[AJ]
# include "TRANSFERT/IBC/commonIBCmotionRel2Abs.h"
	  
	  if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

	  pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam;
	  densPtr[ noind + ideb] = roOut[indR];

	  vxPtr[noind+ideb] = uOut[indR];
	  vyPtr[noind+ideb] = vOut[indR];
	  vzPtr[noind+ideb] = wOut[indR];
        }
    }
  else if (bctype == 1) // adherence (lineaire)
    {
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
	  E_Int indR = rcvPts[noind+ideb];

	  // vitesse
	  u = uOut[indR];
	  v = vOut[indR];
	  w = wOut[indR];

	  //[AJ]
# include "TRANSFERT/IBC/commonIBCmotionAbs2Rel.h"  	    

# include "TRANSFERT/IBC/commonBCType1.h"
	  
	  uOut[indR] = ucible;
	  vOut[indR] = vcible;
	  wOut[indR] = wcible;

	  //[AJ]
# include "TRANSFERT/IBC/commonIBCmotionRel2Abs.h"

	  if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

	  pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam;
	  densPtr[ noind + ideb] = roOut[indR];

	  vxPtr[noind+ideb] = uOut[indR];
	  vyPtr[noind+ideb] = vOut[indR];
	  vzPtr[noind+ideb] = wOut[indR];

        }
    }
  else if (bctype == 2) // loi de paroi log
    {
#   include "TRANSFERT/IBC/pointer.h"

      E_Int err  = 0;
      E_Int skip = 0;
      //initialisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
	{
	  //E_Int indR = rcvPts[noind];
	  E_Int indR = rcvPts[noind+ideb];

	  roext = roOut[indR]; // densite du point interpole
	  text  = tOut[indR];  // pression du point interpole
	  pext  = text*roext*cvgam;

	  // vitesse du pt ext
	  u = uOut[indR];
	  v = vOut[indR];
	  w = wOut[indR];

#       include "TRANSFERT/IBC/commonLogLaw_init.h"
	  // out= utau  et err
	}

      // Newton pour utau
#    include "TRANSFERT/IBC/commonLogLaw_Newton.h"

      //initialisation Newton SA  + vitesse cible
#if NUTILDE_FERRARI == 0
#    include "TRANSFERT/IBC/commonLogLaw_cible.h"
#elif NUTILDE_FERRARI == 1
#    include "TRANSFERT/IBC/nutilde_Ferrari.h"
#else
#    include "TRANSFERT/IBC/nutilde_Ferrari_adim.h"
#endif

      if (nvars == 6)
	{
	  // Newton pour mut
#if NUTILDE_FERRARI == 0
#       include "TRANSFERT/IBC/nutildeSA_Newton.h"
#endif
	  // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
	  for (E_Int noind = 0; noind < ifin-ideb; noind++)
	    {
	      E_Int indR = rcvPts[noind+ideb];


	      // For Post (tOut temperature du point image en entree, pt corrige en sortie)

	      twall = tOut[indR]  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
	      densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
	      pressPtr[noind+ideb]= press_vec[noind ];

	      roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
	      uOut[indR]     = ucible_vec[noind]; vOut[indR] = vcible_vec[noind]; wOut[indR]  = wcible_vec[noind];
	      tOut[indR]     = tcible_vec[noind];
	      varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind];                                         //nutilde*signibc

	      vxPtr[noind+ideb] = uOut[indR];
	      vyPtr[noind+ideb] = vOut[indR];
	      vzPtr[noind+ideb] = wOut[indR];
	    }
	}
      else //5eq
	{
	  // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif
	  for (E_Int noind = 0; noind < ifin-ideb; noind++)
	    {
	      E_Int indR = rcvPts[noind+ideb];

	      roOut[indR]= press_vec[noind ]/tcible_vec[noind]*cvgaminv;

	      uOut[indR] = ucible_vec[noind]; vOut[indR] = vcible_vec[noind]; wOut[indR] = wcible_vec[noind];
	      tOut[indR] = tcible_vec[noind];

	      vxPtr[noind+ideb] = uOut[indR];
	      vyPtr[noind+ideb] = vOut[indR];
	      vzPtr[noind+ideb] = wOut[indR];
	    }
	}
    }
  else if (bctype == 3) // loi de paroi Musker
    {
#   include "TRANSFERT/IBC/pointer.h" 

      E_Int err  = 0;
      E_Int skip = 0; 
      //initialisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif 
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
	{
	  //E_Int indR = rcvPts[noind];
	  E_Int indR = rcvPts[noind+ideb];
 
	  roext = roOut[indR]; // densite @ Image Pnt
	  text  = tOut[indR];  // pression @ Image Pnt
	  pext  = text*roext*cvgam;

	  // vitesse @ Image Pnt
	  u = uOut[indR];
	  v = vOut[indR]; 
	  w = wOut[indR];
	  
# include "TRANSFERT/IBC/commonIBCmotionAbs2Rel.h"
	  //Tangential and Normal velocities: need relative velocity 
#       include "TRANSFERT/IBC/commonMuskerLaw_init.h"
	  // out= utau  et err
	}  

      // Newton pour utau
#    include "TRANSFERT/IBC/commonMuskerLaw_Newton.h" 

      //initialisation Newton SA  + vitesse cible
#if NUTILDE_FERRARI == 0
#    include "TRANSFERT/IBC/commonMuskerLaw_cible.h"
#elif NUTILDE_FERRARI == 1
#    include "TRANSFERT/IBC/nutilde_Ferrari.h"
#else
#    include "TRANSFERT/IBC/nutilde_Ferrari_adim.h"
#endif
      if (nvars == 6)
	{
	  // Newton pour mut
#if NUTILDE_FERRARI == 0
#       include "TRANSFERT/IBC/nutildeSA_Newton.h" 
#endif
	  // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif 
	  for (E_Int noind = 0; noind < ifin-ideb; noind++)
	    {
	      E_Int indR = rcvPts[noind+ideb];

	      // For Post (tOut temperature du point image en entree, pt corrige en sortie)
	      twall = tOut[indR] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
	      densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
	      pressPtr[noind+ideb]= press_vec[noind ];

	      // Mise a jour pt corrige
	      roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;       
	      uOut[indR]     = ucible_vec[noind];
	      vOut[indR]     = vcible_vec[noind];
	      wOut[indR]     = wcible_vec[noind];
	      tOut[indR]     = tcible_vec[noind];
	      varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind];  //nutilde*signibc

# include "TRANSFERT/IBC/commonIBCmotionRel2Abs.h"  

	      vxPtr[noind+ideb] = uOut[indR];
	      vyPtr[noind+ideb] = vOut[indR];
	      vzPtr[noind+ideb] = wOut[indR];

	      

	      // printf("OUT WALL LAW: %f %f %f %f\n",uOut[indR],vOut[indR],wOut[indR],varSAOut[indR]);
	    }
	}
      else //5eq 
	{
	  // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif 
	  for (E_Int noind = 0; noind < ifin-ideb; noind++)
	    {
	      E_Int indR = rcvPts[noind+ideb];

	      // For Post (tOut temperature du point image)
	      twall = tOut[indR]  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
	      densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
	      pressPtr[noind+ideb]= press_vec[noind ];

	      roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;   
	      uOut[indR]     = ucible_vec[noind];
	      vOut[indR]     = vcible_vec[noind];
	      wOut[indR]     = wcible_vec[noind];
	      tOut[indR]     = tcible_vec[noind];

# include "TRANSFERT/IBC/commonIBCmotionRel2Abs.h"  

	      vxPtr[noind+ideb] = uOut[indR];
	      vyPtr[noind+ideb] = vOut[indR];
	      vzPtr[noind+ideb] = wOut[indR];
	    }
	}

    }
  else if (bctype == 4) // outpres
    {
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
	  E_Int indR = rcvPts[noind+ideb];

	  if (densPtr[noind+ideb] > 0.) {
            roOut[indR] = densPtr[noind+ideb];
	  }
	  else{
	    densPtr[noind+ideb] = -roOut[indR];
	  }
	  if (pressPtr[noind+ideb] > 0.) {
            tOut[indR] = pressPtr[noind+ideb]/(roOut[indR]*cvgam);//pext/(roext*cvgam)
	  }
	  
	  //tOut[indR] = pressPtr[noind+ideb]/(roOut[indR]*cvgam);//pext/(roext*cvgam)
	  //densPtr[noind+ideb] = roOut[indR];
	  
	  
	  vxPtr[noind+ideb] = uOut[indR];
	  vyPtr[noind+ideb] = vOut[indR];
	  vzPtr[noind+ideb] = wOut[indR];

	  
        }
    }
  else if (bctype == 5) // inj
    {
      //printf("injection\n");
      if (d1 == NULL)
	{
	  printf("Warning: IBC inj with no data.\n"); return 0; // no data
	}
      E_Int indR, newtonmax, nitnwt;
      E_Float d0x, d0y, d0z, gam, gam1, gam2, rgp, tolnewton, c4, c5, c6, c0, c1, c2, c3;
      E_Float roi, ui, vi, wi, Ti, pi, ro0, u0, v0, w0, T0, p0, wni, ri, residug;
      E_Float tnx, tny, tnz, norm, vni, roc0, usd0n, usd0n2, qen, wn0, wng;
      E_Float rog, vg, Tg, pg, ha, pa, b, vng, dwng, f, df;
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
	{
	  indR = rcvPts[noind+ideb];
	  // std::cout<< "********************************************************************" << std::endl;
	  // printf("ro,p stored %f %f\n", densPtr[noind+ideb], pressPtr[noind+ideb]);
	  // printf("u, v, w: %f %f %f\n", vxPtr[noind+ideb], vyPtr[noind+ideb], vzPtr[noind+ideb]);
	  // printf("u,y stored %f %f\n", utauPtr[noind+ideb], yplusPtr[noind+ideb]);
	  // printf("state: %f %f %f %f %f\n", d1[noind+ideb], d2[noind+ideb], d3[noind+ideb], d4[noind+ideb], d5[noind+ideb]);
	  // std::cout<< "********************************************************************" << std::endl;
	  d0x = d3[noind+ideb];
	  d0y = d4[noind+ideb];
	  d0z = d5[noind+ideb];
	  
	  ha = d1[noind+ideb];
	  pa = d2[noind+ideb];

	  gam = gamma;
	  gam1 = gam - 1.;
	  //gam1_1 = 1. / gam1;
	  gam2 = gam / gam1;
	  rgp = cv*gam1;
	  newtonmax = 40;
	  tolnewton = 1.e-6;
	  //c4 = 5. / 6.;
	  //c5 = 2. / 6.;
	  //c6 = -1. / 6.;
	  //c0 = 1./c6;
	  //c1 =-(c4 + c5)*c0;
	  //c2 =- c6*c0;
	  //c3 = (2.- c5- c4)*c0;

	  // Init newton
	  roi = roOut[indR];
	  ui  = uOut[indR];
	  vi  = vOut[indR];
	  wi  = wOut[indR];
	  Ti  = tOut[indR];
	  pi  = roi*rgp*Ti;

	  ro0 = roOut[indR];
	  u0  = uOut[indR];
	  v0  = vOut[indR];
	  w0  = wOut[indR];
	  T0  = tOut[indR];
	  p0  = ro0*rgp*T0;

	  // normale ext normalisee, a inverser en Euler
	  tnx = xPI[noind+ideb]-xPW[noind+ideb];
	  tny = yPI[noind+ideb]-yPW[noind+ideb];
	  tnz = zPI[noind+ideb]-zPW[noind+ideb];
	  tnx = -tnx; tny = -tny; tnz = -tnz; // sure?
	  norm = sqrt(tnx*tnx+tny*tny+tnz*tnz);
	  norm = 1./K_FUNC::E_max(norm, 1.e-12);
	  tnx = tnx*norm;
	  tny = tny*norm;
	  tnz = tnz*norm;

	  vni  = ui*tnx + vi*tny + wi*tnz;
	  roc0 = sqrt(ro0*gam*p0); // rho*c

	  // sans dimension: produit scalaire direction vitesse . normale
	  usd0n  = 1./(d0x*tnx + d0y*tny + d0z*tnz);
	  usd0n2 = usd0n*usd0n;

	  qen = 0.; // a ajouter en ALE

	  // ...   Inner caracteristic variable
	  // ...   Relative normal velocity
	  wni = vni - qen;
	  ri  = pi + roc0*wni;

	  // ...   Newton Initialization for the relative normal velocity
	  wn0  = u0*tnx + v0*tny + w0*tnz  - qen;
	  wng  = wn0;

	  // resolution Newton
	  residug  = 1.e+20;
	  nitnwt   = 0;

	  while (residug > tolnewton && nitnwt < newtonmax)
	    {
	      nitnwt  += 1;
	      residug = 0.;

	      // LINEARISATION a ajouter??
	      b = 1. - ((wng+qen)*(wng+qen))*usd0n2/(2.*ha);

	      b = K_FUNC::E_max(0.2,b); //cutoff robustesse
	      // p = Pi (1 -U^2/(2CpTi))^(gamma/(gamma-1))
	      pg  = pa*pow(b,gam2);

	      //      nan = isnan(pg)
	      //      if(nan)   write(*,*)'fuck Nan inflow_newton',pa(li),b,gam2

	      rog = gam2*pg/(ha*b);

	      f    = pg + roc0*wng - ri;
	      df   = roc0 - rog*(wng+qen)*usd0n2;

	      dwng = -f/df;
	      wng = wng + dwng;

	      residug = K_FUNC::E_max(residug, K_FUNC::E_abs(dwng/wng));
	    }

	  // LINEARISATION A ajouter??
	  vng = (wng+qen);
	  //...      Absolute velocity module
	  vg  = vng*usd0n;

	  b = 1. - (vng*vng)*usd0n2/(2.*ha);
	  b = K_FUNC::E_max(0.2,b); //cutoff robustesse

	  pg   = pa* pow(b,gam2);
	  rog  = gam2*pg/(ha*b);

	  Tg   = pg/(rog*rgp);

	  roOut[indR] = rog;
	  uOut[indR] = vg*d0x;
	  vOut[indR] = vg*d0y;
	  wOut[indR] = vg*d0z;
	  tOut[indR] = Tg;

	  // update des grandeurs pour post
	  densPtr[noind+ideb] = rog;
	  pressPtr[noind+ideb] = pg;
	  vxPtr[noind+ideb] = uOut[indR];
	  vyPtr[noind+ideb] = vOut[indR];
	  vzPtr[noind+ideb] = wOut[indR];
	}
    }
  else if (bctype == 6) // TBLE
    {
	      if (ithread==1)
		{
		  std::cout << "TBLE pas codee !  " << std::endl;
		}
    }      
  else if (bctype == 7) // loi de paroi adh paroi rotation
    {
#   include "TRANSFERT/IBC/pointer.h" 

      E_Int err  = 0;
      E_Int skip = 0; 

      E_Float teta_out = param_real[ROT_TETA];
      E_Float tetap    = param_real[ROT_TETAP];
      //E_Float teta     = teta_out;
      E_Float teta     = 0;


    E_Float cay,caz,ctheta, stheta,vx,vy,vz,vn_paroi;
    stheta = sin(teta);
    ctheta = cos(teta);
#ifdef _OPENMP4
       #pragma omp simd
#endif
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
         E_Int indR = rcvPts[noind+ideb];

         cay = -stheta*(zPW[noind+ideb] - param_real[ROT_CENTER+2]) + ctheta*(yPW[noind+ideb] - param_real[ROT_CENTER+1]);
         caz =  stheta*(yPW[noind+ideb] - param_real[ROT_CENTER+1]) + ctheta*(zPW[noind+ideb] - param_real[ROT_CENTER+2]);
# include "TRANSFERT/IBC/commonGeom.h"
         vx  =  0;
         vy  = -tetap*caz;
         vz  =  tetap*cay;

         //composante normale de la vitesse paroi
         vn_paroi = vx*n0 + vy*n1 + vz*n2;
  
         //composante tangentielle de la vitesse paroi au pt interpole
         vx = vx-vn_paroi*n0;
         vy = vy-vn_paroi*n1;
         vz = vz-vn_paroi*n2;


         // vitesse relative paroi
         u = uOut[indR]-vx;
         v = vOut[indR]-vy; 
         w = wOut[indR]-vz;
         //if (noind == 0){printf("avt %f %f %f %f \n",vOut[indR],vy,wOut[indR],vz );}
         //
         vn_paroi = u*n0 + v*n1 + w*n2;         

         ucible = (u-vn_paroi*n0)*alphasbeta;// u du pt corrige
         vcible = (v-vn_paroi*n1)*alphasbeta;// v du pt corrige
         wcible = (w-vn_paroi*n2)*alphasbeta;// w du pt corrige


         //E_Float uext2    = vx*vx+vy*vy+vz*vz;
         //E_Float pressure = tOut[indR]*roOut[indR];
         //tOut[indR]    = tOut[indR] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext2); // Crocco-Busemann
         //roOut[indR]   = pressure/tOut[indR];

         uOut[indR] = ucible+vx;
         vOut[indR] = vcible+vy;
         wOut[indR] = wcible+vz;
         //printf("apr %f %f \n",vOut[indR],wOut[indR] );
         if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

         pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam;
         densPtr[ noind + ideb] = roOut[indR];

         vxPtr[noind+ideb] = uOut[indR];
         vyPtr[noind+ideb] = vOut[indR];
         vzPtr[noind+ideb] = wOut[indR];

        }

    }//bctype 
  else if (bctype == 8) // loi de paroi Pohlhausen
    {
#   include "TRANSFERT/IBC/pointer.h"

      E_Int err  = 0;
      E_Int skip = 0;
      //initialisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
	{
	  //E_Int indR = rcvPts[noind];
	  E_Int indR = rcvPts[noind+ideb];

	  roext = roOut[indR]; // densite du point interpole
	  text  = tOut[indR];  // pression du point interpole
	  pext  = text*roext*cvgam;

	  // vitesse du pt ext
	  u = uOut[indR];
	  v = vOut[indR];
	  w = wOut[indR];
#       include "TRANSFERT/IBC/commonPohlhausenLaw_init.h"
	}

#    include "TRANSFERT/IBC/commonPohlhausenLaw_cible.h"

      if (nvars == 6)
	{
	  // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif
	  for (E_Int noind = 0; noind < ifin-ideb; noind++)
	    {
	      E_Int indR = rcvPts[noind+ideb];

	      // For Post (tOut temperature du point image en entree, pt corrige en sortie)
	      twall = tOut[indR] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
	      densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
	      pressPtr[noind+ideb]= press_vec[noind ];

	      // Mise a jour pt corrige
	      roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
	      uOut[indR]     = ucible_vec[noind];
	      vOut[indR]     = vcible_vec[noind];
	      wOut[indR]     = wcible_vec[noind];
	      tOut[indR]     = tcible_vec[noind];
	      varSAOut[indR] = 0.; // laminar

	      vxPtr[noind+ideb] = uOut[indR];
	      vyPtr[noind+ideb] = vOut[indR];
	      vzPtr[noind+ideb] = wOut[indR];

	    }
	}
      else //5eq
	{
	  // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
	  for (E_Int noind = 0; noind < ifin-ideb; noind++)
	    {
	      E_Int indR = rcvPts[noind+ideb];

	      // For Post (tOut temperature du point image)
	      twall = tOut[indR]  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
	      densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
	      pressPtr[noind+ideb]= press_vec[noind ];

	      roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
	      uOut[indR]     = ucible_vec[noind];
	      vOut[indR]     = vcible_vec[noind];
	      wOut[indR]     = wcible_vec[noind];
	      tOut[indR]     = tcible_vec[noind];

	      vxPtr[noind+ideb] = uOut[indR];
	      vyPtr[noind+ideb] = vOut[indR];
	      vzPtr[noind+ideb] = wOut[indR];
	    }
	}

    }//bctype
  else if (bctype == 9) // loi de paroi Thwaites
    {
#   include "TRANSFERT/IBC/pointer.h"

      E_Int err  = 0;
      E_Int skip = 0;
      E_Float* matp;

      //initialisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
	{
	  //E_Int indR = rcvPts[noind];
	  E_Int indR = rcvPts[noind+ideb];

	  roext = roOut[indR]; // densite du point interpole
	  text  = tOut[indR];  // pression du point interpole
	  pext  = text*roext*cvgam;

	  // vitesse du pt ext
	  u = uOut[indR];
	  v = vOut[indR];
	  w = wOut[indR];
#       include "TRANSFERT/IBC/commonThwaitesLaw_init.h"
	  // out= utau  et err
	}

#    include "TRANSFERT/IBC/commonThwaitesLaw_cible.h"

      if (nvars == 6)
	{
	  // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif
	  for (E_Int noind = 0; noind < ifin-ideb; noind++)
	    {
	      E_Int indR = rcvPts[noind+ideb];

	      // For Post (tOut temperature du point image en entree, pt corrige en sortie)
	      twall = tOut[indR] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
	      densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
	      pressPtr[noind+ideb]= press_vec[noind ];

	      // Mise a jour pt corrige
	      roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
	      uOut[indR]     = ucible_vec[noind];
	      vOut[indR]     = vcible_vec[noind];
	      wOut[indR]     = wcible_vec[noind];
	      tOut[indR]     = tcible_vec[noind];
	      varSAOut[indR] = 0.;

	      vxPtr[noind+ideb] = uOut[indR];
	      vyPtr[noind+ideb] = vOut[indR];
	      vzPtr[noind+ideb] = wOut[indR];

	    }
	}
      else //5eq
	{
	  // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
	  for (E_Int noind = 0; noind < ifin-ideb; noind++)
	    {
	      E_Int indR = rcvPts[noind+ideb];

	      // For Post (tOut temperature du point image)
	      twall = tOut[indR]  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
	      densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
	      pressPtr[noind+ideb]= press_vec[noind ];

	      roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
	      uOut[indR]     = ucible_vec[noind];
	      vOut[indR]     = vcible_vec[noind];
	      wOut[indR]     = wcible_vec[noind];
	      tOut[indR]     = tcible_vec[noind];

	      vxPtr[noind+ideb] = uOut[indR];
	      vyPtr[noind+ideb] = vOut[indR];
	      vzPtr[noind+ideb] = wOut[indR];

	    }
	}
    }//bctype
  else if (bctype == 10) // loi de paroi Mafzal
    {
#   include "TRANSFERT/IBC/pointer.h"

      E_Float MafzalMode = 3; // param_real[ MAFZAL_MODE ];

      E_Int err  = 0;
      E_Int skip = 0;
      E_Float tgradU = 0.;
      E_Float ngradU = 0.;
      E_Float unext = 0.;
      //initialisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
	{
	  //E_Int indR = rcvPts[noind];
	  E_Int indR = rcvPts[noind+ideb];

	  roext = roOut[indR]; // densite du point interpole
	  text  = tOut[indR];  // pression du point interpole
	  pext  = text*roext*cvgam;

	  // vitesse du pt ext
	  u = uOut[indR];
	  v = vOut[indR];
	  w = wOut[indR];

	  // gradient de pression au point ext
	  gradxPext = gradxPressPtr[noind+ideb];
	  gradyPext = gradyPressPtr[noind+ideb];
	  gradzPext = gradzPressPtr[noind+ideb];

#       include "TRANSFERT/IBC/commonMafzalLaw_init.h"
	}

      // PREMIERE PASSE MUSKER POUR UTAU_ORI #####################################
#    include "TRANSFERT/IBC/commonMuskerLaw_Newton.h"
      err  = 0;
      skip = 0;
      E_Float alphaMafMus;

#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
	{
	  utau0 = utauOri_vec[noind];
	  utauOri_vec[noind] = utau_vec[noind];
	  utau_vec[noind] = utau0;

#  include "TRANSFERT/IBC/mafzal_vec.h"
	}
      // #########################################################################

      // Newton pour utau -> Mafzal
#    include "TRANSFERT/IBC/commonMafzalLaw_Newton.h"

      //initialisation Newton SA  + vitesse cible
#if NUTILDE_FERRARI == 0
#    include "TRANSFERT/IBC/commonMafzalLaw_cible.h"
#elif NUTILDE_FERRARI == 1
#    include "TRANSFERT/IBC/nutilde_Ferrari_Mafzal.h"
#else
#    include "TRANSFERT/IBC/nutilde_Ferrari_adim_Mafzal.h"
#endif
      if (nvars == 6)
	{
	  // Newton pour mut
#if NUTILDE_FERRARI == 0
#       include "TRANSFERT/IBC/nutildeSA_Newton.h"
#endif
	  // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif
	  for (E_Int noind = 0; noind < ifin-ideb; noind++)
	    {
	      E_Int indR = rcvPts[noind+ideb];

	      // For Post (tOut temperature du point image en entree, pt corrige en sortie)
	      twall = tOut[indR] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
	      densPtr[noind+ideb]  = press_vec[noind]/twall*cvgaminv;
	      pressPtr[noind+ideb] = press_vec[noind];

	      // Mise a jour pt corrige
	      roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
	      uOut[indR]     = ucible_vec[noind];
	      vOut[indR]     = vcible_vec[noind];
	      wOut[indR]     = wcible_vec[noind];
	      tOut[indR]     = tcible_vec[noind];
	      varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind];  //nutilde*signibc

	      vxPtr[noind+ideb] = uOut[indR];
	      vyPtr[noind+ideb] = vOut[indR];
	      vzPtr[noind+ideb] = wOut[indR];
	    }
	}
      else //5eq
	{
	  // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
	  for (E_Int noind = 0; noind < ifin-ideb; noind++)
	    {
	      E_Int indR = rcvPts[noind+ideb];

	      // For Post (tOut temperature du point image)
	      twall = tOut[indR]  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
	      densPtr[noind+ideb]  = press_vec[noind]/twall*cvgaminv;
	      pressPtr[noind+ideb] = press_vec[noind];

	      roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
	      uOut[indR]     = ucible_vec[noind];
	      vOut[indR]     = vcible_vec[noind];
	      wOut[indR]     = wcible_vec[noind];
	      tOut[indR]     = tcible_vec[noind];

	      vxPtr[noind+ideb] = uOut[indR];
	      vyPtr[noind+ideb] = vOut[indR];
	      vzPtr[noind+ideb] = wOut[indR];
	    }
	}

    }//bctype
  else if (bctype == 11) // TBLE_FULL
    {
      printf("Warning !!! setIBCTransfersCommonVar2: bcType %d TBLE_FULL  not implemented.\n", bctype);

    }
  else if (bctype == 12) // isothermal - prescribed wall temp
    {
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
	  E_Int indR = rcvPts[noind+ideb];

	  // get values at image point (stored in FlowSolution#Centers)
	  u = uOut[indR];     
	  v = vOut[indR];     
	  w = wOut[indR];     
	  t = tOut[indR]; 

	  //Calc values at target points
# include "TRANSFERT/IBC/commonBCType1.h"
	  tcible = t*alphasbeta + (1-alphasbeta)*tempExtraPtr[ noind + ideb];

	  // set values at target point(stored in FlowSolution#Centers)
	  uOut[indR] = ucible; 
	  vOut[indR] = vcible; 
	  wOut[indR] = wcible; 
	  tOut[indR] = tcible; 
	  if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

	  // set values in tc
	  pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam; 
	  densPtr[ noind + ideb] = roOut[indR];                   
	  tempPtr[ noind + ideb] = tOut[indR];                    
	  vxPtr[noind+ideb] = uOut[indR];                         
	  vyPtr[noind+ideb] = vOut[indR];                         
	  vzPtr[noind+ideb] = wOut[indR];                         
        }
    }
  else if (bctype == 13) // heat flux - prescribed heat flux
    {
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
	  E_Int indR = rcvPts[noind+ideb];

	  // get values at image point (stored in FlowSolution#Centers)
	  u = uOut[indR];     
	  v = vOut[indR];     
	  w = wOut[indR];     
	  t = tOut[indR]; 

	  //Calc values at target points
# include "TRANSFERT/IBC/commonBCType1.h"
	  if (int(copysign(1.0,alpha)) == int(copysign(1.0,beta)))
	    {      
	      //Target points are inside the fluid
	      // second order one sided difference w/ first order one sided different
	      // ∂T/∂n=q=(T_T - T_W)/(x_T-x_W) -> T_W=T_T-q*(x_T-x_W)
	      // use above in
	      //  ∂T/∂n=q=(-T_IP*(x_T-x_W)^2+T_T*(x_IP-x_W)^2-T_W*[(x_IP-x_W)^2-(x_T-x_W)^2])/[(x_T-x_W)(x_I-x_W)(x_I-x_T)]
	      E_Float absquared=pow((alpha+beta),2);
	      E_Float asquared =pow(alpha,2);
	      E_Float tmp_var  =absquared-asquared;
	      tcible = (tempExtraPtr[ noind + ideb]*alpha*beta*(alpha+beta)+asquared*t-alpha*tempExtraPtr[ noind + ideb]*tmp_var)/(absquared-tmp_var);
	    }
	  else
	    {
	      //Target points are inside the solid
	      // second order central difference
	      // ∂T/∂n=q=(T_IP - T_T)/(x_IP-x_T) -> T_T=T_IP-q*(x_IP-x_T) 
	      tcible = t-(abs(alpha)+abs(beta))*tempExtraPtr[ noind + ideb];
	    }
	  

	  // set values at target point(stored in FlowSolution#Centers)
	  uOut[indR] = ucible; 
	  vOut[indR] = vcible; 
	  wOut[indR] = wcible; 
	  tOut[indR] = tcible; 
	  if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

	  // set values in tc
	  pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam; 
	  densPtr[ noind + ideb] = roOut[indR];                   
	  tempPtr[ noind + ideb] = tOut[indR];                    
	  vxPtr[noind+ideb] = uOut[indR];                         
	  vyPtr[noind+ideb] = vOut[indR];                         
	  vzPtr[noind+ideb] = wOut[indR];                         
        }
    }
  else if (bctype == 140) // Wire Model - M. Terracol & E. Monaha 2021 - Numerical Wire Mesh Model for the Simulation of Noise Reduction Devices
    {
#   include "TRANSFERT/IBC/pointer.h"
      E_Int err  = 0;
      E_Int skip = 0;
      
      E_Float ro_Pnt2, u_Pnt2, v_Pnt2, w_Pnt2, t_Pnt2;
      E_Float p_w,p_t,s_w,t_t;
      E_Float ro_up,u_up,v_up,w_up;
      E_Float u_n_w,u_t_w;
      E_Float v_n_w,v_t_w;
      E_Float w_n_w,w_t_w;

      E_Float* roOut_Pnt2    = densPtr +  5*nbRcvPts;// ro          
      E_Float* uOut_Pnt2     = densPtr +  6*nbRcvPts;// u           
      E_Float* vOut_Pnt2     = densPtr +  7*nbRcvPts;// v            
      E_Float* wOut_Pnt2     = densPtr +  8*nbRcvPts;// w           
      E_Float* tOut_Pnt2     = densPtr +  9*nbRcvPts;// temperature 
      E_Float* varSAOut_Pnt2 = densPtr + 10*nbRcvPts;// pseudoviscosity nu_tilde 
        
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
	{
	  E_Int indR    = rcvPts[noind+ideb];
	  	  
	  // get values at image point (stored in FlowSolution#Centers)
	  // These image point values are those for tc
	  // 1pnt
	  ro            = roOut[indR]; 
	  u             = uOut[indR];     
	  v             = vOut[indR];     
	  w             = wOut[indR];     
	  t             = tOut[indR];

	  // These image point values are those for tc2
	  // Value at A2 for the wire model (second image points)
	  // 2pnt
	  ro_Pnt2 = roOut_Pnt2[noind+ ideb]; 
	  u_Pnt2  = uOut_Pnt2[noind+ ideb];     
	  v_Pnt2  = vOut_Pnt2[noind+ ideb];     
	  w_Pnt2  = wOut_Pnt2[noind+ ideb];     
	  t_Pnt2  = tOut_Pnt2[noind+ ideb];
	  
# include "TRANSFERT/IBC/commonGeom.h"
	  s_w = u*n0+v*n1+w*n2;
	  s_w = copysign(1,s_w);
	  // Flow -->
	  //           |
	  // s_w=-1    |  s_w=1
	  //           |
	  //   x   n<--|
	  //   ↑
	  //  target pnt
	  if (s_w<0){
	    ro_up=ro;
	    u_up=u;
	    v_up=v;
	    w_up=w;
	  }
	  else{
	    ro_up=ro_Pnt2;
	    u_up=u_Pnt2;
	    v_up=v_Pnt2;
	    w_up=w_Pnt2;
	  }

	  p_w    = 0.5*R_gas*(ro_Pnt2*t_Pnt2+ro*t);
	  p_t    = p_w-s_w*0.25*K_wire*ro_up*(u_up*u_up+v_up*v_up+w_up*w_up);
	  tcible = p_t/(ro_up*R_gas);
	  twall  = p_w/(ro_up*R_gas);

	  roOut[indR]= ro_up;
	  tOut[indR] = tcible;

	  //save in tc 
	  pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam;	  
	  densPtr[ noind + ideb] = roOut[indR];

	  // vec is to save the value for the following loop of the target points only
	  // save wall values as interpolation is of interpolated & wall are used to determine
	  // target values
	  ro_vec[noind]=ro_up;	  
	  mu_vec[noind]=coefSuth * sqrt(K_FUNC::E_abs(twall)*Tsinv) / (1.+Cs/twall);

	  // UP TO HERE CORRECT	  
	  
	  //n0=-n0;n1=-n1;n2=-n2; // needed to get normals in correct direction
# include "TRANSFERT/IBC/normalTangentVelocity.h"
	  t0=ut/uext;t1=vt/uext;t2=wt/uext;
	  
	  uext_image  =uext;
	  uscaln_image=uscaln;
	  
	  u=u_up;v=v_up;w=w_up;	
# include "TRANSFERT/IBC/normalTangentVelocity.h"
	  uext_wall  = uext;
	  uscaln_wall= uscaln;
	  uext_wall  = Delta_V_wire*uext_wall;
	  

	  //uext      = alphasbeta*uext_image  +(1.-alphasbeta)*uext_wall;
	  //uscaln    = alphasbeta*uscaln_image+(1.-alphasbeta)*uscaln_wall;
	  
	  //ucible = uscaln*n0+uext*t0;
	  //vcible = uscaln*n1+uext*t1;
	  //wcible = uscaln*n2+uext*t2;

	  ucible      = alphasbeta*uOut[indR]  +(1.-alphasbeta)*(uscaln_wall*n0+uext_wall*t0);
	  vcible      = alphasbeta*vOut[indR]  +(1.-alphasbeta)*(uscaln_wall*n1+uext_wall*t1);
	  wcible      = alphasbeta*wOut[indR]  +(1.-alphasbeta)*(uscaln_wall*n2+uext_wall*t2);

	  //if (yPC[noind+ideb]>0 && yPC[noind+ideb]<0.0007){
	  //  printf("xPC[noind+ideb] uscaln_wall uscaln_target:: %g %g %g %g %g\n",xPC[noind+ideb],uscaln_wall,uscaln_target,ucible,vcible);
	  //}	  
	  
	  uOut[indR] = ucible; 
	  vOut[indR] = vcible; 
	  wOut[indR] = wcible; 	  

	  //save in tc 
	  vxPtr[noind+ideb] = uOut[indR];                         
	  vyPtr[noind+ideb] = vOut[indR];                         
	  vzPtr[noind+ideb] = wOut[indR];

	  //vec is to save the value for the following loop of the target points only	  
	  ut_vec[noind]=uext_wall*t0;
	  vt_vec[noind]=uext_wall*t1;	
	  wt_vec[noind]=uext_wall*t2;

	  ucible_vec[noind]=uscaln_wall*n0;
	  vcible_vec[noind]=uscaln_wall*n1;	
	  wcible_vec[noind]=uscaln_wall*n2;
	}
	  
      if (nvars==6){
	E_Int count = 0;
	for (E_Int noind = 0; noind < ifin-ideb; noind++)
	{
	  E_Int indR    = rcvPts[noind+ideb];

	  //|u_w|=sqrt(un,w²+ut,w²)
	  E_Float norm_wall_vel = sqrt(ucible_vec[noind]*ucible_vec[noind]+vcible_vec[noind]*vcible_vec[noind]+wcible_vec[noind]*wcible_vec[noind]+
				       ut_vec[noind]    *ut_vec[noind]    +vt_vec[noind]    *vt_vec[noind]    +wt_vec[noind]    *wt_vec[noind]);
	  
	  E_Float nu_t_local    = Ct_WM*Diam_wire*norm_wall_vel; // mu_t = rho*C_t*d*|u_w| & nu_t=mu_t/rho

	  nutcible_vec[noind] = nu_t_local;
	  nutilde             = K_FUNC::E_abs( nutcible_vec[noind] );
	  aa_vec[noind]       = nutilde; //to start first iteration
	  ut_vec[noind]       = nutilde; // Sauvegarde du nutilde, au cas ou newton non convergent.
	                                 // save  ṽ(O)=v_t (v=nu)  	  
# include "TRANSFERT/IBC/fnutilde_vec.h"
	  
	}
	// Newton pour ṽ
#       include "TRANSFERT/IBC/nutildeSA_Newton.h"
	
	for (E_Int noind = 0; noind < ifin-ideb; noind++)
	{
	  E_Int indR    = rcvPts[noind+ideb];
# include "TRANSFERT/IBC/commonGeom.h"	  
	  varSAOut_Pnt2[noind+ ideb]=nutcible_vec[noind];
	  varSAOut[indR]=alphasbeta*varSAOut[indR]  +(1.-alphasbeta)*(aa_vec[noind]);
	}
      }
    }
  else if (bctype == 141) // Wire Model prt2 - just interpolation of image points & placed at control points
    {;}      
  else
    {
      printf("Warning !!! setIBCTransfersCommonVar2: bcType " SF_D_ " not implemented.\n", bctype);
      return 0;
    }

 WireMeshSkip:
  return 1;
}

