c     Produit [ +/- A(l1) + rho(l1) I ]/2 . dW(l1) pour une direction
c     (transcription de FastS/Compute/LU/lu_i_3dfull.for).
c     Entrees : l1 (cellule voisine), lt1, lt2 (ses deux faces),
c               idir (2,3,4 : rayon spectral coe(l1,idir)), signe,
c               macro TN (ti, tj ou tk), drodm_out, rop, coe
c     Sorties : b1..b5
        tcx = 0.5*(TN(lt2,1)+TN(lt1,1))
        tcy = 0.5*(TN(lt2,2)+TN(lt1,2))
        tcz = 0.5*(TN(lt2,3)+TN(lt1,3))

        u     = rop(l1,2)
        v     = rop(l1,3)
        w     = rop(l1,4)
        t     = rop(l1,5)
        q2    = .5*(u*u+v*v+w*w)
        h     = cp*t + q2
        ph2   = gamm1*q2
        qn    = tcx*u+tcy*v+tcz*w
        rad   = coe(l1,idir)*signe

        b11=rad
        b12=tcx
        b13=tcy
        b14=tcz
        b21=(tcx*ph2-u*qn)
        b22=(qn-tcx*gam2*u)+rad
        b23=(tcy*u-tcx*gamm1*v)
        b24=(tcz*u-tcx*gamm1*w)
        b25=tcx*gamm1
        b31=(tcy*ph2-v*qn)
        b32=(tcx*v-tcy*gamm1*u)
        b33=(qn-tcy*gam2*v)+rad
        b34=(tcz*v-tcy*gamm1*w)
        b35=tcy*gamm1
        b41=(tcz*ph2-w*qn)
        b42=(tcx*w-tcz*gamm1*u)
        b43=(tcy*w-tcz*gamm1*v)
        b44=(qn-tcz*gam2*w)+rad
        b45=tcz*gamm1
        b51=qn*(ph2 - h)
        b52=(tcx*h-gamm1*u*qn)
        b53=(tcy*h-gamm1*v*qn)
        b54=(tcz*h-gamm1*w*qn)
        b55=(gam1*qn)+rad

        b1= b11*drodm_out(l1,1)+b12*drodm_out(l1,2)+b13*drodm_out(l1,3)
     &     +b14*drodm_out(l1,4)
        b2= b21*drodm_out(l1,1)+b22*drodm_out(l1,2)+b23*drodm_out(l1,3)
     &     +b24*drodm_out(l1,4)+b25*drodm_out(l1,5)
        b3= b31*drodm_out(l1,1)+b32*drodm_out(l1,2)+b33*drodm_out(l1,3)
     &     +b34*drodm_out(l1,4)+b35*drodm_out(l1,5)
        b4= b41*drodm_out(l1,1)+b42*drodm_out(l1,2)+b43*drodm_out(l1,3)
     &     +b44*drodm_out(l1,4)+b45*drodm_out(l1,5)
        b5= b51*drodm_out(l1,1)+b52*drodm_out(l1,2)+b53*drodm_out(l1,3)
     &     +b54*drodm_out(l1,4)+b55*drodm_out(l1,5)
