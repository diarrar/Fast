         pos = param_int_tc[shift_rac + nrac*17];
         E_Int* fluD = param_int_tc + pos + 4*nbflu;
         E_Int idir    = fluD[0];
         E_Int sizefluD= fluD[1];
         E_Int nobcD   = fluD[2];
         E_Int nobcR   = fluD[3];
         E_Int ratio_i = fluD[4];
         E_Int ratio_j = fluD[5];
         E_Int ratio_k = fluD[6];

         E_Int pt_bcs= param_int[ NoD ][PT_BC];
         E_Int nb_bc = param_int[ NoD ][pt_bcs];
         E_Int pt_bc = param_int[ NoD ][pt_bcs + 1 + nobcD];
         E_Int adrFlu= param_int[ NoD ][pt_bcs + 1 + nobcD + nb_bc] +3;
         E_Int* fen  = param_int[ NoD ] + pt_bc + BC_FEN;
         E_Float* fluxD = param_real[ NoD ] + adrFlu;

         pt_bcs= param_int[ NoR ][PT_BC];
         nb_bc = param_int[ NoR ][pt_bcs];
         //adrFlu= param_int[ NoR ][pt_bcs + 1 + nobcR + nb_bc];
         adrFlu= param_int[ NoR ][pt_bcs + 1 + nobcR + nb_bc] +3;
         E_Float* fluxR = param_real[ NoR ] + adrFlu;
         
        //if(NoR==50) {printf("NoR:, %d,NoD: %d, idir: %d , nobcD: %d , nobcR: %d , ratioijk: %d %d %d  szfluD: %d \n", NoR, NoD, idir,nobcD,nobcR, ratio_i,ratio_j,ratio_k, sizefluD) ;}
         // Cas 2D
         if(param_int[ NoR ][NIJK+4]==0)
           { 
            E_Int sizefluR= sizefluD/2;
            if(idir<=2)
              {
               E_Int jmax= (fen[3]-fen[2]+1)/2;
               #pragma omp for nowait
               for (E_Int j = 0; j < jmax; j++)
                {
                  E_Int lr  = j;
                  E_Int ld0 = j*2;
                  E_Int ld1 = j*2 +1;
                 //if(NoR==48) {printf("jr %d, kr %d, lr %d , ld0 %d , ld1 %d \n", j, k, lr, ld0, ld1);}
                  fluxR[lr            ] = fluxD[ld0            ];
                  fluxR[lr            ]+= fluxD[ld1            ];
                  fluxR[lr +sizefluR  ] = fluxD[ld0 +sizefluD  ];
                  fluxR[lr +sizefluR  ]+= fluxD[ld1 +sizefluD  ];
                  fluxR[lr +sizefluR*2] = fluxD[ld0 +sizefluD*2];
                  fluxR[lr +sizefluR*2]+= fluxD[ld1 +sizefluD*2];
                  fluxR[lr +sizefluR*3] = fluxD[ld0 +sizefluD*3];
                  fluxR[lr +sizefluR*3]+= fluxD[ld1 +sizefluD*3];
                  fluxR[lr +sizefluR*4] = fluxD[ld0 +sizefluD*4];
                  fluxR[lr +sizefluR*4]+= fluxD[ld1 +sizefluD*4];
                  //printf("fluR %f %f %f %f %f \n", fluxR[lr],fluxR[lr +sizefluR],fluxR[lr +sizefluR*2], fluxR[lr +sizefluR*3],fluxR[lr +sizefluR*4]);
                }
              }
            else if(idir<=4)
              {
               E_Int imax= (fen[1]-fen[0]+1)/2;
               #pragma omp for nowait
               for (E_Int i = 0; i < imax; i++)
                {
                  E_Int lr  = i;
                  E_Int ld0 = i*2;
                  E_Int ld1 = i*2 +1;
                 //if(NoR==28) {printf("ir %d, lr %d , ld0 %d , ld1 %d , sizeRD %d %d, nobcR %d\n", i, lr, ld0, ld1,sizefluR,sizefluD, nobcR);}
                
                  fluxR[lr            ] = fluxD[ld0            ];
                  fluxR[lr            ]+= fluxD[ld1            ];
                  fluxR[lr +sizefluR  ] = fluxD[ld0 +sizefluD  ];
                  fluxR[lr +sizefluR  ]+= fluxD[ld1 +sizefluD  ];
                  fluxR[lr +sizefluR*2] = fluxD[ld0 +sizefluD*2];
                  fluxR[lr +sizefluR*2]+= fluxD[ld1 +sizefluD*2];
                  fluxR[lr +sizefluR*3] = fluxD[ld0 +sizefluD*3];
                  fluxR[lr +sizefluR*3]+= fluxD[ld1 +sizefluD*3];
                  fluxR[lr +sizefluR*4] = fluxD[ld0 +sizefluD*4];
                  fluxR[lr +sizefluR*4]+= fluxD[ld1 +sizefluD*4];
      //if(NoR==28)  printf("fluR %f %f %f %f %f %f \n", fluxR[lr],fluxR[lr +sizefluR],fluxR[lr +sizefluR*2],fluxR[lr +sizefluR*4],fluxD[ld0],fluxD[ld1] );
                }
              }
           }
         else  //3D
           { 
            //E_Int sizefluR= sizefluD/4;
            if(idir<=2)
             {
              E_Int sizefluR= sizefluD/(ratio_j*ratio_k);
              E_Int jmax= (fen[3]-fen[2]+1)/ratio_j;
              E_Int kmax= (fen[5]-fen[4]+1)/ratio_k;
              E_Int jmaxD = jmax*ratio_j;
              E_Int incjD = 1;
              E_Float coefj = 1;
              E_Float coefk = 1;
              if(ratio_j==1){incjD=0; coefj = 0;}
              E_Int inckD = jmax*ratio_j;
              if(ratio_k==1){inckD=0; coefk = 0;}
              
              #pragma omp for collapse(2) nowait
              for (E_Int k = 0; k < kmax; k++)
               {
               for (E_Int j = 0; j < jmax; j++)
                {
                 E_Int lr  = j   + k*jmax;
                 //E_Int ld0 = j*2 + k*2*inck;
                 E_Int ld0 = j*ratio_j + k*ratio_k*jmaxD;
                 E_Int ld1 = ld0 + incjD;
                 E_Int ld2 = ld0 + inckD;
                 E_Int ld3 = ld1 + inckD;
            //if(NoR==50) printf("indice0 %d %d %d %d %d ratiojk: %d %d \n", lr, ld0, ld3, j,k, ratio_j, ratio_k);
            //fflush(stdout);

                 fluxR[lr            ] = fluxD[ld0];
                 fluxR[lr            ]+= fluxD[ld1]*coefj;
                 fluxR[lr            ]+= fluxD[ld2]*coefk;
                 fluxR[lr            ]+= fluxD[ld3]*coefj*coefk;
                 fluxR[lr +sizefluR  ] = fluxD[ld0 +sizefluD  ];
                 fluxR[lr +sizefluR  ]+= fluxD[ld1 +sizefluD  ]*coefj;
                 fluxR[lr +sizefluR  ]+= fluxD[ld2 +sizefluD  ]*coefk;;
                 fluxR[lr +sizefluR  ]+= fluxD[ld3 +sizefluD  ]*coefj*coefk;
                 fluxR[lr +sizefluR*2] = fluxD[ld0 +sizefluD*2];
                 fluxR[lr +sizefluR*2]+= fluxD[ld1 +sizefluD*2]*coefj;
                 fluxR[lr +sizefluR*2]+= fluxD[ld2 +sizefluD*2]*coefk;;
                 fluxR[lr +sizefluR*2]+= fluxD[ld3 +sizefluD*2]*coefj*coefk;
                 fluxR[lr +sizefluR*3] = fluxD[ld0 +sizefluD*3];
                 fluxR[lr +sizefluR*3]+= fluxD[ld1 +sizefluD*3]*coefj;
                 fluxR[lr +sizefluR*3]+= fluxD[ld2 +sizefluD*3]*coefk;;
                 fluxR[lr +sizefluR*3]+= fluxD[ld3 +sizefluD*3]*coefj*coefk;
                 fluxR[lr +sizefluR*4] = fluxD[ld0 +sizefluD*4];
                 fluxR[lr +sizefluR*4]+= fluxD[ld1 +sizefluD*4]*coefj;
                 fluxR[lr +sizefluR*4]+= fluxD[ld2 +sizefluD*4]*coefk;;
                 fluxR[lr +sizefluR*4]+= fluxD[ld3 +sizefluD*4]*coefj*coefk;
            //printf("indice1 %d %d %d %d %d \n", lr, ld0, ld3, j,k);
      //if(NoR==57)  printf("fluR %f %f %f %f %f %f \n", fluxR[lr],fluxR[lr +sizefluR],fluxR[lr +sizefluR*2],fluxR[lr +sizefluR*4],fluxD[ld0],fluxD[ld1] );
      //if(NoR==50)  printf("fluRi %f %f %f %f %f %d %d  \n", fluxR[lr],fluxD[ld0],fluxD[ld1]*coefj, fluxD[ld2]*coefk,fluxD[ld3]*coefj*coefk, j,k );
      //      fflush(stdout);
                }
               }
             }
            else if(idir<=4)
             {
              E_Int sizefluR= sizefluD/(ratio_i*ratio_k);
              E_Int imax= (fen[1]-fen[0]+1)/ratio_i;
              E_Int kmax= (fen[5]-fen[4]+1)/ratio_k;
              //E_Int inck = imax*2;
              E_Int imaxD = imax*ratio_i;
              E_Int inciD = 1;
              E_Float coefi = 1;
              E_Float coefk = 1;
              if(ratio_i==1){inciD=0; coefi = 0;}
              E_Int inckD = imax*ratio_i;
              if(ratio_k==1){inckD=0;coefk = 0;}

              #pragma omp for collapse(2) nowait
              for (E_Int k = 0; k < kmax; k++)
               {
               for (E_Int i = 0; i < imax; i++)
                {
                 E_Int lr  = i   + k*imax;
                 E_Int ld0 = i*ratio_i + k*ratio_k*imaxD;
                 E_Int ld1 = ld0 +inciD ;
                 E_Int ld2 = ld0 +inckD;
                 E_Int ld3 = ld1 +inckD;

            //if(NoR==50) printf("indice0j %d %d %d %d %d ratiojk: %d %d \n", lr, ld0, ld3, j,k, ratio_i, ratio_k);
                 fluxR[lr            ] = fluxD[ld0];
                 fluxR[lr            ]+= fluxD[ld1]*coefi;
                 fluxR[lr            ]+= fluxD[ld2]*coefk;
                 fluxR[lr            ]+= fluxD[ld3]*coefi*coefk;
                 fluxR[lr +sizefluR  ] = fluxD[ld0 +sizefluD  ];
                 fluxR[lr +sizefluR  ]+= fluxD[ld1 +sizefluD  ]*coefi;
                 fluxR[lr +sizefluR  ]+= fluxD[ld2 +sizefluD  ]*coefk;;
                 fluxR[lr +sizefluR  ]+= fluxD[ld3 +sizefluD  ]*coefi*coefk;
                 fluxR[lr +sizefluR*2] = fluxD[ld0 +sizefluD*2];
                 fluxR[lr +sizefluR*2]+= fluxD[ld1 +sizefluD*2]*coefi;
                 fluxR[lr +sizefluR*2]+= fluxD[ld2 +sizefluD*2]*coefk;;
                 fluxR[lr +sizefluR*2]+= fluxD[ld3 +sizefluD*2]*coefi*coefk;
                 fluxR[lr +sizefluR*3] = fluxD[ld0 +sizefluD*3];
                 fluxR[lr +sizefluR*3]+= fluxD[ld1 +sizefluD*3]*coefi;
                 fluxR[lr +sizefluR*3]+= fluxD[ld2 +sizefluD*3]*coefk;;
                 fluxR[lr +sizefluR*3]+= fluxD[ld3 +sizefluD*3]*coefi*coefk;
                 fluxR[lr +sizefluR*4] = fluxD[ld0 +sizefluD*4];
                 fluxR[lr +sizefluR*4]+= fluxD[ld1 +sizefluD*4]*coefi;
                 fluxR[lr +sizefluR*4]+= fluxD[ld2 +sizefluD*4]*coefk;;
                 fluxR[lr +sizefluR*4]+= fluxD[ld3 +sizefluD*4]*coefi*coefk;

      //if(NoR==50)  printf("fluRj %f %f %f %f %f %d %d  \n", fluxR[lr],fluxD[ld0],fluxD[ld1]*coefi, fluxD[ld2]*coefk,fluxD[ld3]*coefi*coefk, i,k );
        //    fflush(stdout);
                }
               }
             }
            else 
             {
              E_Int sizefluR= sizefluD/(ratio_i*ratio_j);
              E_Int imax= (fen[1]-fen[0]+1)/ratio_i;
              E_Int jmax= (fen[3]-fen[2]+1)/ratio_j;
              E_Int imaxD = imax*ratio_i;
              E_Int inciD = 1;
              E_Float coefi = 1;
              E_Float coefj = 1;
              if(ratio_i==1){inciD=0; coefi = 0;}
              E_Int incjD = imax*ratio_i;
              if(ratio_j==1){incjD=0; coefj = 0;}

              #pragma omp for collapse(2) nowait
              for (E_Int j = 0; j < jmax; j++)
               {
               for (E_Int i = 0; i < imax; i++)
                {
                 E_Int lr  = i   + j*imax;
                 E_Int ld0 = i*ratio_i + j*ratio_j*imaxD;
                 E_Int ld1 = ld0 +inciD ;
                 E_Int ld2 = ld0 +incjD;
                 E_Int ld3 = ld1 +incjD;


                 fluxR[lr            ] = fluxD[ld0];
                 fluxR[lr            ]+= fluxD[ld1]*coefi;
                 fluxR[lr            ]+= fluxD[ld2]*coefj;
                 fluxR[lr            ]+= fluxD[ld3]*coefi*coefj;
                 fluxR[lr +sizefluR  ] = fluxD[ld0 +sizefluD  ];
                 fluxR[lr +sizefluR  ]+= fluxD[ld1 +sizefluD  ]*coefi;
                 fluxR[lr +sizefluR  ]+= fluxD[ld2 +sizefluD  ]*coefj;;
                 fluxR[lr +sizefluR  ]+= fluxD[ld3 +sizefluD  ]*coefi*coefj;
                 fluxR[lr +sizefluR*2] = fluxD[ld0 +sizefluD*2];
                 fluxR[lr +sizefluR*2]+= fluxD[ld1 +sizefluD*2]*coefi;
                 fluxR[lr +sizefluR*2]+= fluxD[ld2 +sizefluD*2]*coefj;;
                 fluxR[lr +sizefluR*2]+= fluxD[ld3 +sizefluD*2]*coefi*coefj;
                 fluxR[lr +sizefluR*3] = fluxD[ld0 +sizefluD*3];
                 fluxR[lr +sizefluR*3]+= fluxD[ld1 +sizefluD*3]*coefi;
                 fluxR[lr +sizefluR*3]+= fluxD[ld2 +sizefluD*3]*coefj;;
                 fluxR[lr +sizefluR*3]+= fluxD[ld3 +sizefluD*3]*coefi*coefj;
                 fluxR[lr +sizefluR*4] = fluxD[ld0 +sizefluD*4];
                 fluxR[lr +sizefluR*4]+= fluxD[ld1 +sizefluD*4]*coefi;
                 fluxR[lr +sizefluR*4]+= fluxD[ld2 +sizefluD*4]*coefj;;
                 fluxR[lr +sizefluR*4]+= fluxD[ld3 +sizefluD*4]*coefi*coefj;
                }
               }
             }
           }//2d/3d

