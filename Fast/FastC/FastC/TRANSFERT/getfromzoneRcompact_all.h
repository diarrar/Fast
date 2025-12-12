       //
       // Calcul dimension zone receveuse
       //
       PyObject* sol;  PyObject* t;
       ipt_roR[nd]= NULL; ipt_qR[nd]= NULL; ipt_SR[nd]= NULL; ipt_psiGR[nd]= NULL;
       sol = K_PYTREE::getNodeFromName1(zoneR , "FlowSolution#Centers");
       if(sol != NULL)
       {  
         for (E_Int ivar = 0; ivar < nbvar_inlist; ivar++)
          {
             if (ivar==0)
               { t           = K_PYTREE::getNodeFromName1(sol, varname );
                 ipt_roR[nd] = K_PYTREE::getValueAF(t, hook);
               }
             else if (ivar==1)
               { t = K_PYTREE::getNodeFromName1(sol, varname1);
                 if (t != NULL) ipt_qR[nd] = K_PYTREE::getValueAF(t, hook);
               }
             else if (ivar==2)
               { t = K_PYTREE::getNodeFromName1(sol, varname2);
                 if (t != NULL) ipt_SR[nd] = K_PYTREE::getValueAF(t, hook);
               }
             else if (ivar==3)
               { t = K_PYTREE::getNodeFromName1(sol, varname3);
                 if (t != NULL) ipt_psiGR[nd] = K_PYTREE::getValueAF(t, hook);
               }
          }
       }

       ipt_roD_vert[nd]= NULL; ipt_qD_vert[nd]= NULL; ipt_SD_vert[nd]= NULL; ipt_psiGD_vert[nd]=NULL;

       PyObject* own      = K_PYTREE::getNodeFromName1(zoneR , ".Solver#ownData");
               t          = K_PYTREE::getNodeFromName1(own, "Parameter_int");
       ipt_param_intR[nd] = K_PYTREE::getValueAI(t, hook);

               t          = K_PYTREE::getNodeFromName1(own, "Parameter_real");
       ipt_param_realR[nd]= K_PYTREE::getValueAF(t, hook);

