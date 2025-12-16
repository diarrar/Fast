    for (E_Int nopass = 1; nopass <= npass_transfer; nopass++)
    {
      //printf("transfert : npass %d %d \n",nopass, npass_transfer ); fflush(0);
      char str_real[16];char str_int[15];
      if      (nopass == 1)
       { strcpy(str_int, "param_int_tc1"); strcpy(str_real, "param_real_tc1");
         pyParam_int_tc1  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc1 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc1 , param_int_tc1 );  int_tc[ nopass-1 ]= param_int_tc1 -> begin();
         if ( pyParam_real_tc1 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc1, param_real_tc1); real_tc[ nopass-1 ]= param_real_tc1-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 2)
       { strcpy(str_int, "param_int_tc2"); strcpy(str_real, "param_real_tc2");
         pyParam_int_tc2  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc2 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc2 , param_int_tc2 );  int_tc[ nopass-1 ]= param_int_tc2 -> begin();
         if ( pyParam_real_tc2 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc2, param_real_tc2); real_tc[ nopass-1 ]= param_real_tc2-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 3)
       { strcpy(str_int, "param_int_tc3"); strcpy(str_real, "param_real_tc3");
         pyParam_int_tc3  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc3 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc3 , param_int_tc3 );  int_tc[ nopass-1 ]= param_int_tc3 -> begin();
         if ( pyParam_real_tc3 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc3, param_real_tc3); real_tc[ nopass-1 ]= param_real_tc3-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 4) 
       { strcpy(str_int, "param_int_tc4"); strcpy(str_real, "param_real_tc4");
         pyParam_int_tc4  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc4 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc4 , param_int_tc4 );  int_tc[ nopass-1 ]= param_int_tc4 -> begin();
         if ( pyParam_real_tc4 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc4, param_real_tc4); real_tc[ nopass-1 ]= param_real_tc4-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else {printf("transfert : npass > 4 pas codee."); return NULL;}
    }
    if (npass_transfer == 0){int_tc[0] = NULL; real_tc[0] = NULL; }
