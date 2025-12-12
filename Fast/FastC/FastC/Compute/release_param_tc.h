    for (E_Int nopass = 1; nopass <= npass_transfer; nopass++)
    {
      if      (nopass == 1)
       { 
         RELEASESHAREDN( pyParam_int_tc1, param_int_tc1); 
         if ( pyParam_real_tc1 != Py_None) { RELEASESHAREDN( pyParam_real_tc1, param_real_tc1);}
       }
      else if (nopass == 2)
       { 
         RELEASESHAREDN( pyParam_int_tc2, param_int_tc2); 
         if ( pyParam_real_tc2 != Py_None) { RELEASESHAREDN( pyParam_real_tc2, param_real_tc2);}
       }
      else if (nopass == 3)
       { 
         RELEASESHAREDN( pyParam_int_tc3, param_int_tc3); 
         if ( pyParam_real_tc3 != Py_None) { RELEASESHAREDN( pyParam_real_tc3, param_real_tc3);}
       }
      else if (nopass == 4) 
       { 
         RELEASESHAREDN( pyParam_int_tc4, param_int_tc4); 
         if ( pyParam_real_tc4 != Py_None) { RELEASESHAREDN( pyParam_real_tc4, param_real_tc4);}
       }
      else {printf("transfert : npass > 4 pas codee."); return NULL;}
    }
