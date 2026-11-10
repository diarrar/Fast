c     Flux de Roe, reconstruction o3, sur la face (lf - 1/2) dans la
c     direction d increment inc. Metrique: macro TN (ti, tj ou tk).
c     Transcription des includes FastS:
c       Normale/normale_3dfull_*, Slope/o3_slope_var, etat_roe_GD,
c       Vit_ent/qn_3dfull_*, Vit_ent/fludiffer_3dfull_*
c     Entrees : lf, ltf, inc, rop, TN, c4,c5,c6, gam, gam1, rgp, psiroe
c     Sorties : flu1..flu5
        tcx = TN(ltf +v1mtr)
        tcy = TN(ltf +v2mtr)
        tcz = TN(ltf +v3mtr)
        si  = sqrt (tcx*tcx + tcy*tcy + tcz*tcz)

        nm  = lf -  inc
        nm2 = lf -2*inc
        np  = lf +  inc

c       pentes: qm etat droit, qp etat gauche
        qm1=  c4*rop(lf +v1) +  c5*rop(nm  +v1) +  c6*rop(np +v1)
        qp1=  c4*rop(nm +v1) +  c6*rop(nm2 +v1) +  c5*rop(lf +v1)
        qm2=  c4*rop(lf +v2) +  c5*rop(nm  +v2) +  c6*rop(np +v2)
        qp2=  c4*rop(nm +v2) +  c6*rop(nm2 +v2) +  c5*rop(lf +v2)
        qm3=  c4*rop(lf +v3) +  c5*rop(nm  +v3) +  c6*rop(np +v3)
        qp3=  c4*rop(nm +v3) +  c6*rop(nm2 +v3) +  c5*rop(lf +v3)
        qm4=  c4*rop(lf +v4) +  c5*rop(nm  +v4) +  c6*rop(np +v4)
        qp4=  c4*rop(nm +v4) +  c6*rop(nm2 +v4) +  c5*rop(lf +v4)
        qm5=  c4*rop(lf +v5) +  c5*rop(nm  +v5) +  c6*rop(np +v5)
        qp5=  c4*rop(nm +v5) +  c6*rop(nm2 +v5) +  c5*rop(lf +v5)

c       etat gauche / droit et moyenne de Roe (etat_roe_GD)
        r1=qp1
        p1=r1*qp5*rgp
        h1=gam1*qp5*rgp + .5*(qp2*qp2+qp3*qp3+qp4*qp4)

        r2 =qm1
        p2=r2*qm5*rgp
        h2=gam1*qm5*rgp + .5*(qm2*qm2+qm3*qm3+qm4*qm4)

        qp1=sqrt(r1)
        qm1=sqrt(r2)
        r  = qp1*qm1
        r_1= 1./(qp1+qm1)
        u =(qp1*qp2+qm1*qm2)*r_1
        v =(qp1*qp3+qm1*qm3)*r_1
        w =(qp1*qp4+qm1*qm4)*r_1
        h =(qp1*h1 +qm1*h2 )*r_1
        q =.5*(u*u+v*v+w*w)
        c =sqrt((gam-1.)*abs(h-q))

c       vitesse normale (qn_3dfull)
        qn1=qp2*tcx+qp3*tcy+qp4*tcz
        qn2=qm2*tcx+qm3*tcy+qm4*tcz

c       dissipation de Roe + partie centree (fludiffer_3dfull)
        du =(qm2-qp2)*si
        dv =(qm3-qp3)*si
        dw =(qm4-qp4)*si
        dp = p2-p1
        dqn=qn2-qn1

        s_1 =1.0/max(si,1.E-30)
        nx = tcx*s_1
        ny = tcy*s_1
        nz = tcz*s_1

        qn =u*nx+ v*ny+ w*nz

        f1 = qn
        f4 = psiroe*(abs(f1) + c)

        f2 =0.5 + SIGN(0.5,ABS(f1)-f4)
        f3 =(f2-1.0)*(f1*f1/f4 + f4)*0.5  -f2*abs(f1)

        f1  = (r2-r1 - dp/(c*c))*si

        flu1= f3*( f1                                   )
        flu2= f3*( f1 * u + r * (du - nx*dqn)           )
        flu3= f3*( f1 * v + r * (dv - ny*dqn)           )
        flu4= f3*( f1 * w + r * (dw - nz*dqn)           )
        flu5= f3*( f1 * q + r * (u*du+v*dv+w*dw-qn*dqn) )

        f1 = qn+c
        f2 =0.5+sign(0.5,abs(f1)-f4)
        f3 =(f2-1.)*(f1*f1/f4 + f4)*0.5  -f2*abs(f1)

        f1 = (dp+r*c*dqn*s_1) / (2.*c*c)*si*f3

        flu1 = flu1 + f1
        flu2 = flu2 + f1 * (u + nx* c)
        flu3 = flu3 + f1 * (v + ny* c)
        flu4 = flu4 + f1 * (w + nz* c)
        flu5 = flu5 + f1 * (h + qn *c)

        f1 = qn-c
        f2 =0.5+sign(0.5,abs(f1)-f4)
        f3 =(f2-1.)*(f1*f1/f4 + f4)*0.5  -f2*abs(f1)

        f1 = (dp-r*c*dqn*s_1) / (2.*c*c)*f3*si

        flu1 = flu1 + f1
        flu2 = flu2 + f1 * (u - nx* c)
        flu3 = flu3 + f1 * (v - ny* c)
        flu4 = flu4 + f1 * (w - nz* c)
        flu5 = flu5 + f1 * (h - qn* c)

        rou1=r1*qn1
        rou2=r2*qn2
        p1p2=p1+p2

        flu1 = (flu1 + rou1     + rou2                 )*0.5
        flu2 = (flu2 + rou1*qp2 + rou2*qm2 + tcx * p1p2)*0.5
        flu3 = (flu3 + rou1*qp3 + rou2*qm3 + tcy * p1p2)*0.5
        flu4 = (flu4 + rou1*qp4 + rou2*qm4 + tcz * p1p2)*0.5
        flu5 = (flu5 + rou1*h1  + rou2*h2              )*0.5
