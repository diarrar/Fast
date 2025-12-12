         flag_correct_flu=0
          
         !Wall model boundary condition
            mobile_coef = 1.
            if (nbdata.ne.0) mobile_coef = param_real( iptdata )

            if(idir.le.2) THEN
             incijk     =1
             call bflwallmodel(ndom,ithread, idir, nitcfg,
     &                    param_int(NEQ_IJ), mobile_coef,
     &                    param_int, param_real,
     &                    incijk, ind_CL,rop,drodm, ti,venti,xmut,x,y,z)
            elseif(idir.le.4) THEN
             incijk     = param_int(NIJK)
             call bflwallmodel(ndom,ithread, idir, nitcfg,
     &                    param_int(NEQ_IJ), mobile_coef,
     &                    param_int, param_real,
     &                    incijk, ind_CL, rop,drodm,tj,ventj,xmut,x,y,z)
            else
             incijk     = param_int(NIJK+1)*param_int(NIJK)
             call bflwallmodel(ndom,ithread, idir, nitcfg,
     &                    param_int(NEQ_K), mobile_coef,
     &                    param_int, param_real,
     &                    incijk, ind_CL, rop,drodm,tk,ventk,xmut,x,y,z)
            endif

            !Mise a jour flag pour calcul correct des efforts
            flag_correct_flu=1
