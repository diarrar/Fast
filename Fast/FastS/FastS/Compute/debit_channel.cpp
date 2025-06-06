     E_Float surf=0; E_Float deb =0;
     E_Int Nbre_thread_actif = __NUMTHREADS__;
     //calcul debit channel flow
    E_Int nbtask_loc = ipt_omp[0]; 
     for (E_Int ntask = 0; ntask < nbtask_loc; ntask++)
      {
       E_Int ptiter_loc = ipt_omp[nssiter];
       E_Int pttask = ptiter_loc + ntask*(6+Nbre_thread_actif*7);
       E_Int nd = ipt_omp[ pttask ];
       if (param_int[nd][CHANNEL_BODYFORCE]==1)
         {
            E_Int Nbre_thread_actif_loc = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
            for (E_Int ith = 0; ith < Nbre_thread_actif_loc; ith++)
             {
               E_Float* ipt_forcage  = forcage  + ith*2+ nd*2*Nbre_thread_actif;
               deb  += ipt_forcage[0];
               //printf("surf %f %f %d \n",surf, ipt_forcage[1], ith);
               surf += ipt_forcage[1];
             }
         }
      }
      for (E_Int ntask = 0; ntask < nbtask_loc; ntask++)
      {
       E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
       E_Int nd = ipt_omp[ pttask ];
       if (param_int[nd][CHANNEL_BODYFORCE]==1)
         { param_real[nd][FORCAGE_DEBN] = -1*(deb/surf - param_real[nd][FORCAGE_DEBN+1]);
           if(nitcfg==1 && ntask==0){printf("debit: %15.10f debit0: %15.10f , surf: %f , forcage: %15.10f it: %d %d \n", deb/surf,param_real[0][FORCAGE_DEB0], surf, param_real[0][FORCAGE_DEBN], nitcfg, nitrun);}
         }
      }
