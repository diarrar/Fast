    E_Int nbtask_loc = ipt_omp[0]; 
    for (E_Int ntask = 0; ntask < nbtask_loc; ntask++)
      {
        E_Int ptiter_loc = ipt_omp[nssiter];
        E_Int pttask     = ptiter_loc + ntask*(6+Nbre_thread_actif*7);
        E_Int nd         = ipt_omp[ pttask ];
        ithread_loc      = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
        E_Int nd_subzone      = ipt_omp[ pttask + 1 ];
        Nbre_thread_actif_loc = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
        //ipt_inddm_omp         = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;

        if (ithread_loc == -1) {continue;}

        if( param_int[nd][CHANNEL_BODYFORCE]==1)
           { 
             E_Float* ipt_forcage     = forcage   + (ithread_loc-1)*2+ nd*2*Nbre_thread_actif; 
             E_Int* ipt_ind_dm_thread = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;
 
             //initialisation debit et surface pour channel flow
             ipt_forcage[0]=0; ipt_forcage[1]=0;

             E_Int mode =1;
             src_channel_(nd,ithread, param_int[nd], param_real[nd],
                          ipt_ind_dm_thread , nitcfg, nitrun, mode,
                          ipt_forcage, iptdrodm ,iptro_ssiter[nd],  iptk[nd],  iptcoe);
           }
      }
