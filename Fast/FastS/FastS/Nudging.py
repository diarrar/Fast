import Converter.Internal as Internal
import Generator.PyTree as G
import Geom.PyTree as D
import Converter.PyTree as C
import Connector.PyTree as X
import numpy as numpy
import math

#transfert du terme de forcage calculé sur grille capteur vers la grille ns
def _computeNudging(nudging, zones, nstep):

  (t_cap, data_caps, it_loc, dt_scale, tc_cap2ns, tc_ns2cap, zone_tg, nud_mod, verbose) = nudging
  varsNud=['VelocityX_P1','VelocityY_P1']
  if nstep==1: varsNud =['VelocityX','VelocityY']
  if nud_mod ==0:
      #Etape I: interp sol NS vers position capteur de reference
      if tc_ns2cap is not None: _interpDataNS2SensorGrid( zones , t_cap, tc_ns2cap, vars=varsNud)
      ## Etape II
      #interpolation temporelle (eventuelle) de la valeur de reference et calcul terme forcage= Vns -vref sur la grille capteur
      if tc_ns2cap is not None: _nudgingInterpTargetValues(t_cap, data_caps, it_loc, dt_scale, vars=varsNud, verbose=verbose)
      else:
        Dim=2 
        if param_int_firstZone[4] !=0: Dim=3
        _nudgingTargetValues(zones, zone_tg, data_caps, it_loc, dt_scale, nstep=nstep, mode='gauss2', Dim=Dim)
      ## Etape III
      # transfert du terme de forcage calculé sur grille capteur vers la grille ns
      if tc_cap2ns is not None: _interpForcage2grilleNS(zones, t_cap, tc_cap2ns, nstep=nstep)
  else:
     if zone_tg is not None:
       tps_10=Time.time()
       Dim=2 
       if param_int_firstZone[4] !=0: Dim=3
       _nudgingTargetValues( zones, zone_tg, zone_tg , it_loc, dt_scale, nstep=nstep, mode='gauss2', Dim=Dim)

       #fasts._computePT_laplacien(zones, metrics, hook1)
       if nstep==1 : tps_filtre=0.
       #Nud.filter(t, zone_tg, nstep)
       #Nud.filterSquareDuct(t, zone_tg, nstep)
       tps_filtre+=Time.time()-tps_10
       if nstep == nitmax: print("cpu Filter", tps_filtre)


#transfert du terme de forcage calculé sur grille capteur vers la grille ns
def _interpForcage2grilleNS(t, t_cap, tc_cap2ns, nstep=1):

  zonesD = Internal.getZones(tc_cap2ns)
  zonesR = []
  for z in zonesD: 
     zonesR.append( Internal.getNodeFromName(t_cap,z[0]) )
  vars=['MomentumX_src','MomentumY_src']
  for v in vars: C._cpVars(zonesR, v, zonesD, v)
  #for v in vars: C._cpVars(zonesR, 'centers:'+v, zonesD, v)

  #C.convertPyTree2File(zonesD, 'Kaptor'+str(nstep)+'.cgns')
  X.setInterpTransfers(t , zonesD, variables=vars)
  #C.convertPyTree2File(t, 'forcage'+str(nstep)+'.cgns')

#transfert des donnees NS vers la grille capteur 
# t        : NS tree
# t_cap    : sensor tree
# tc_ns2cap: arbre connectivite pour transfert ns vers sensor
def _interpDataNS2SensorGrid(t, t_cap, tc_ns2cap, vars):

    zonesD_ns2cap = Internal.getZones(tc_ns2cap)
    zonesR_ns2cap=[]
    for z in zonesD_ns2cap: zonesR_ns2cap.append( Internal.getNodeFromName(t,z[0]) )

    for v in vars: C._cpVars(zonesR_ns2cap, 'centers:'+v, zonesD_ns2cap, v)

    X.setInterpTransfers(t_cap , zonesD_ns2cap, variables=vars)

    for z in Internal.getZones(t_cap):
      sol =  Internal.getNodeFromName1(z,'FlowSolution')
      if vars[0]== 'VelocityX': 
        vx  = Internal.getNodeFromName(sol,'VelocityX')[1]
        vy  = Internal.getNodeFromName(sol,'VelocityY')[1]
      else:
        vx  = Internal.getNodeFromName(sol,'VelocityX_P1')[1]
        vy  = Internal.getNodeFromName(sol,'VelocityY_P1')[1]
      if z[0]=='cylindreH':
         vx[:,0]=0
         vy[:,0]=0
      elif z[0]=='cylindreB':
         vx[:,1]=0
         vy[:,1]=0
      elif z[0]=='sillage':
         vx[0,1]=0
         vy[0,1]=0
     
#determination valeur au point de forcage par direc ou gauss
def _nudgingInterpTargetValues(t_cap, data_caps, it, dt_scale, vars=['VelocityX'], verbose=0):

   #boucle sur zones capteur
   for z in Internal.getZones(t_cap):
      sol =  Internal.getNodeFromName1(z,'FlowSolution')
      mx  = Internal.getNodeFromName(sol,'MomentumX_src')[1]
      my  = Internal.getNodeFromName(sol,'MomentumY_src')[1]
      if vars[0]== 'VelocityX': 
        vx  = Internal.getNodeFromName(sol,'VelocityX')[1]
        vy  = Internal.getNodeFromName(sol,'VelocityY')[1]
      else:
        vx  = Internal.getNodeFromName(sol,'VelocityX_P1')[1]
        vy  = Internal.getNodeFromName(sol,'VelocityY_P1')[1]

      data_cap = data_caps[ z[0] ]

      coef = float((it%dt_scale))/float(dt_scale)
      t1 = it//dt_scale
      t2 = t1+1
      #mx = (Vitesse_ns- vitesse_reference) au capteur (i,j) et a l'instant entre t2 et t1
      mx[:,:] =  (data_cap[:,:, t1,0 ]  + (data_cap[:,:, t2 ,0 ]-data_cap[:,:, t1,0 ])*coef) - vx[:,:]
      my[:,:] =  (data_cap[:,:, t1,1 ]  + (data_cap[:,:, t2 ,1 ]-data_cap[:,:, t1,1 ])*coef) - vy[:,:]

      #if z[0]=='cylindre':
      #  tg = (data_cap[0,1, t1,0 ]  + (data_cap[0,1, t2 ,0 ]-data_cap[0,1, t1,0 ])*coef)
      #  print( "forcage", tg, vx[0,1], mx[0,1])
      sh = numpy.shape(data_cap)
      if verbose==1:
        for j in range(sh[1]):
          for i in range(sh[0]):
             tgx = (data_cap[i,j, t1,0 ]  + (data_cap[i,j, t2 ,0 ]-data_cap[i,j, t1,0 ])*coef)
             tgy = (data_cap[i,j, t1,1 ]  + (data_cap[i,j, t2 ,1 ]-data_cap[i,j, t1,1 ])*coef)
             #print( "refxy: ", tgx, tgy, 'ns_xy: ',  vx[2,2], vy[2,2], 'src_xy: ', mx[2,2], my[2,2] )
             print( "refxy: ", tgx, tgy, 'ns_xy: ',  vx[i,j], vy[i,j], 'src_xy: ', mx[i,j], my[i,j], z[0],i+j*sh[0],'F'  )


#determination valeur au point de forcage par direc ou gauss
def _nudgingTargetValues(t, capteurs, data_cap, it, dt_scale, nstep=1, mode='dirac', Dim=2, init=False):

  level=''
  if nstep!=1:level='_P1'

  if Dim ==2:
    for key in capteurs:
      z   = Internal.getNodeFromName(t,key)
      sol =  Internal.getNodeFromName1(z,'FlowSolution#Centers')
      mx  = Internal.getNodeFromName(sol,'MomentumX_src')[1]
      my  = Internal.getNodeFromName(sol,'MomentumY_src')[1]
      vx  = Internal.getNodeFromName(sol,'VelocityX'+level)[1]
      vy  = Internal.getNodeFromName(sol,'VelocityY'+level)[1]

      sh  = numpy.shape(mx)
      capts= capteurs[key]
      coef = float((it%dt_scale))/float(dt_scale)
      t1 = it//dt_scale
      t2 = t1+1

      if mode=='dirac':
        for cap in capts:
          ind = cap[0]; Nocapt = cap[1]
          #if Nocapt==10: print( key, "coef:", coef, it,'t12', t1,t2 )
          j= ind//sh[0]
          i= ind-j* sh[0]
          mx[i,j,0]= data_cap[Nocapt, t1,0 ] + (data_cap[Nocapt, t2 ,0 ]-data_cap[Nocapt, t1,0 ])*coef
          my[i,j,0]= data_cap[Nocapt, t1,1 ] + (data_cap[Nocapt, t2 ,1 ]-data_cap[Nocapt, t1,1 ])*coef
      elif mode=='gauss1':
        vx  = Internal.getNodeFromName(sol,'VelocityX')[1]
        vy  = Internal.getNodeFromName(sol,'VelocityY')[1]
        for cap in capts:
          ind = cap[0]; Nocapt = cap[1]
          j= ind//sh[0]
          i= ind-j* sh[0]
          mx[i,j,0]= data_cap[Nocapt, t1,0 ] + (data_cap[Nocapt, t2 ,0 ]-data_cap[Nocapt, t1,0 ])*coef
          my[i,j,0]= data_cap[Nocapt, t1,1 ] + (data_cap[Nocapt, t2 ,1 ]-data_cap[Nocapt, t1,1 ])*coef

          mx[i-1,j,0]= 0.5*(vx[i-2,j,0]+mx[i,j,0])
          mx[i+1,j,0]= 0.5*(vx[i+2,j,0]+mx[i,j,0])
          my[i-1,j,0]= 0.5*(vy[i-2,j,0]+my[i,j,0])
          my[i+1,j,0]= 0.5*(vy[i+2,j,0]+my[i,j,0])
 
          mx[i,j-1,0]= 0.5*(vx[i,j-2,0]+mx[i,j,0])
          mx[i,j+1,0]= 0.5*(vx[i,j+2,0]+mx[i,j,0])
          my[i,j-1,0]= 0.5*(vy[i,j-2,0]+my[i,j,0])
          my[i,j+1,0]= 0.5*(vy[i,j+2,0]+my[i,j,0])

      elif mode=='gauss3':

        amor =0.8
        for cap in capts:
          ind = cap[0]; Nocapt = cap[1]
          j= ind//sh[0]
          i= ind-j* sh[0]

          mx[i,j,0]= data_cap[Nocapt, t1,0 ] + (data_cap[Nocapt, t2 ,0 ]-data_cap[Nocapt, t1,0 ])*coef
          my[i,j,0]= data_cap[Nocapt, t1,1 ] + (data_cap[Nocapt, t2 ,1 ]-data_cap[Nocapt, t1,1 ])*coef

          delta = mx[i,j,0] - vx[i,j,0]

          mx[i-1,j  ,0]= vx[i-1,j  ,0]+ delta*amor
          mx[i-1,j+1,0]= vx[i-1,j+1,0]+ delta*amor*0.9 
          mx[i-1,j-1,0]= vx[i-1,j-1,0]+ delta*amor*0.9
          mx[i+1,j  ,0]= vx[i+1,j  ,0]+ delta*amor
          mx[i+1,j+1,0]= vx[i+1,j+1,0]+ delta*amor*0.9
          mx[i+1,j-1,0]= vx[i+1,j-1,0]+ delta*amor*0.9
          mx[i  ,j-1,0]= vx[i  ,j-1,0]+ delta*amor
          mx[i  ,j+1,0]= vx[i  ,j+1,0]+ delta*amor

          mx[i-2,j  ,0]= vx[i-2,j  ,0]+ delta*amor**2
          mx[i+2,j  ,0]= vx[i+2,j  ,0]+ delta*amor**2 
          mx[i-2,j+1,0]= vx[i-2,j+1,0]+ delta*amor**2*0.9
          mx[i-2,j-1,0]= vx[i-2,j-1,0]+ delta*amor**2*0.9
          mx[i+2,j+1,0]= vx[i+2,j+1,0]+ delta*amor**2*0.9
          mx[i+2,j-1,0]= vx[i+2,j-1,0]+ delta*amor**2*0.9
          mx[i-2,j+2,0]= vx[i-2,j+2,0]+ delta*amor**2*0.6
          mx[i-2,j-2,0]= vx[i-2,j-2,0]+ delta*amor**2*0.6
          mx[i+2,j+2,0]= vx[i+2,j+2,0]+ delta*amor**2*0.6
          mx[i+2,j-2,0]= vx[i+2,j-2,0]+ delta*amor**2*0.6

          mx[i-1,j-2,0]= vx[i-1,j-2,0]+ delta*amor**2*0.9
          mx[i  ,j-2,0]= vx[i  ,j-2,0]+ delta*amor**2
          mx[i+1,j-2,0]= vx[i+1,j-2,0]+ delta*amor**2*0.9
          mx[i-1,j+2,0]= vx[i-1,j+2,0]+ delta*amor**2*0.9
          mx[i  ,j+2,0]= vx[i  ,j+2,0]+ delta*amor**2
          mx[i+1,j+2,0]= vx[i+1,j+2,0]+ delta*amor**2*0.9

          delta = my[i,j,0] - vy[i,j,0]

          my[i-1,j  ,0]= vy[i-1,j  ,0]+ delta*amor
          my[i-1,j+1,0]= vy[i-1,j+1,0]+ delta*amor*0.9 
          my[i-1,j-1,0]= vy[i-1,j-1,0]+ delta*amor*0.9
          my[i+1,j  ,0]= vy[i+1,j  ,0]+ delta*amor
          my[i+1,j+1,0]= vy[i+1,j+1,0]+ delta*amor*0.9
          my[i+1,j-1,0]= vy[i+1,j-1,0]+ delta*amor*0.9
          my[i  ,j-1,0]= vy[i  ,j-1,0]+ delta*amor
          my[i  ,j+1,0]= vy[i  ,j+1,0]+ delta*amor

          my[i-2,j  ,0]= vy[i-2,j  ,0]+ delta*amor**2
          my[i+2,j  ,0]= vy[i+2,j  ,0]+ delta*amor**2 
          my[i-2,j+1,0]= vy[i-2,j+1,0]+ delta*amor**2*0.9
          my[i-2,j-1,0]= vy[i-2,j-1,0]+ delta*amor**2*0.9
          my[i+2,j+1,0]= vy[i+2,j+1,0]+ delta*amor**2*0.9
          my[i+2,j-1,0]= vy[i+2,j-1,0]+ delta*amor**2*0.9
          my[i-2,j+2,0]= vy[i-2,j+2,0]+ delta*amor**2*0.6
          my[i-2,j-2,0]= vy[i-2,j-2,0]+ delta*amor**2*0.6
          my[i+2,j+2,0]= vy[i+2,j+2,0]+ delta*amor**2*0.6
          my[i+2,j-2,0]= vy[i+2,j-2,0]+ delta*amor**2*0.6

          my[i-1,j-2,0]= vy[i-1,j-2,0]+ delta*amor**2*0.9
          my[i  ,j-2,0]= vy[i  ,j-2,0]+ delta*amor**2
          my[i+1,j-2,0]= vy[i+1,j-2,0]+ delta*amor**2*0.9
          my[i-1,j+2,0]= vy[i-1,j+2,0]+ delta*amor**2*0.9
          my[i  ,j+2,0]= vy[i  ,j+2,0]+ delta*amor**2
          my[i+1,j+2,0]= vy[i+1,j+2,0]+ delta*amor**2*0.9

      else:

        for cap in capts:
          ind = cap[0]; Nocapt = cap[1]
          j= ind//sh[0]
          i= ind-j* sh[0]
          mx[i,j,0]= data_cap[Nocapt, t1,0 ] + (data_cap[Nocapt, t2 ,0 ]-data_cap[Nocapt, t1,0 ])*coef
          my[i,j,0]= data_cap[Nocapt, t1,1 ] + (data_cap[Nocapt, t2 ,1 ]-data_cap[Nocapt, t1,1 ])*coef

          mx[i-1,j  ,0]= 0.5*(vx[i-2,j  ,0]+mx[i,j,0])
          mx[i-1,j+1,0]= 0.5*(vx[i-2,j+2,0]+mx[i,j,0])
          mx[i-1,j-1,0]= 0.5*(vx[i-2,j-2,0]+mx[i,j,0])
          mx[i+1,j  ,0]= 0.5*(vx[i+2,j  ,0]+mx[i,j,0])
          mx[i+1,j+1,0]= 0.5*(vx[i+2,j+2,0]+mx[i,j,0])
          mx[i+1,j-1,0]= 0.5*(vx[i+2,j-2,0]+mx[i,j,0])
          mx[i  ,j-1,0]= 0.5*(vx[i  ,j-2,0]+mx[i,j,0])
          mx[i  ,j+1,0]= 0.5*(vx[i  ,j+2,0]+mx[i,j,0])

          my[i-1,j  ,0]= 0.5*(vy[i-2,j  ,0]+my[i,j,0])
          my[i-1,j+1,0]= 0.5*(vy[i-2,j+2,0]+my[i,j,0])
          my[i-1,j-1,0]= 0.5*(vy[i-2,j-2,0]+my[i,j,0])
          my[i+1,j  ,0]= 0.5*(vy[i+2,j  ,0]+my[i,j,0])
          my[i+1,j+1,0]= 0.5*(vy[i+2,j+2,0]+my[i,j,0])
          my[i+1,j-1,0]= 0.5*(vy[i+2,j-2,0]+my[i,j,0])
          my[i  ,j-1,0]= 0.5*(vy[i  ,j-2,0]+my[i,j,0])
          my[i  ,j+1,0]= 0.5*(vy[i  ,j+2,0]+my[i,j,0])

  #3D: pas d'interp temporelle et pas de nom de zone !!! (a revoir)
  else:

    z= t[0]
    sol =  Internal.getNodeFromName2(z,'FlowSolution#Centers')
    mx  = Internal.getNodeFromName(sol,'MomentumX_src')[1]
    my  = Internal.getNodeFromName(sol,'MomentumY_src')[1]
    mz  = Internal.getNodeFromName(sol,'MomentumZ_src')[1]
    vx  = Internal.getNodeFromName(sol,'VelocityX'+level)[1]
    vy  = Internal.getNodeFromName(sol,'VelocityY'+level)[1]
    vz  = Internal.getNodeFromName(sol,'VelocityZ'+level)[1]

    coef=1.
    if init : coef=0
    Ncapteur = numpy.shape(data_cap)[0]
    if mode =='dirac':
      for l in range(Ncapteur):
        i = int(data_cap[l,6])
        j = int(data_cap[l,7])
        k = int(data_cap[l,8])

        mx[i,j,k]= data_cap[l,3 ]*coef 
        my[i,j,k]= data_cap[l,4 ]*coef 
        mz[i,j,k]= data_cap[l,5 ]*coef 

        if l==2: 
           #print( key, "coef:", coef, it,'t12', t1,t2 )
            print( "refxy: ",  mx[i,j,0], my[i,j,0])

    elif mode=='gauss1':
      for l in range(Ncapteur):
        i = int(data_cap[l,6])
        j = int(data_cap[l,7])
        k = int(data_cap[l,8])
        
        #print(l,i,j,k, numpy.shape(vx) )

        mx[i,j,k]= data_cap[l,3 ]*coef 
        my[i,j,k]= data_cap[l,4 ]*coef 
        mz[i,j,k]= data_cap[l,5 ]*coef 

        mx[i-1,j,k]= coef*0.5*(vx[i-2,j,k]+mx[i,j,k])
        mx[i+1,j,k]= coef*0.5*(vx[i+2,j,k]+mx[i,j,k])
        my[i-1,j,k]= coef*0.5*(vy[i-2,j,k]+my[i,j,k])
        my[i+1,j,k]= coef*0.5*(vy[i+2,j,k]+my[i,j,k])
        mz[i-1,j,k]= coef*0.5*(vz[i-2,j,k]+mz[i,j,k])
        mz[i+1,j,k]= coef*0.5*(vz[i+2,j,k]+mz[i,j,k])

        mx[i,j-1,k]= coef*0.5*(vx[i,j-2,k]+mx[i,j,k])
        mx[i,j+1,k]= coef*0.5*(vx[i,j+2,k]+mx[i,j,k])
        my[i,j-1,k]= coef*0.5*(vy[i,j-2,k]+my[i,j,k])
        my[i,j+1,k]= coef*0.5*(vy[i,j+2,k]+my[i,j,k])
        mz[i,j-1,k]= coef*0.5*(vz[i,j-2,k]+mz[i,j,k])
        mz[i,j+1,k]= coef*0.5*(vz[i,j+2,k]+mz[i,j,k])

        mx[i,j,k-1]= coef*0.5*(vx[i,j,k-2]+mx[i,j,k])
        mx[i,j,k+1]= coef*0.5*(vx[i,j,k+2]+mx[i,j,k])
        my[i,j,k-1]= coef*0.5*(vy[i,j,k-2]+my[i,j,k])
        my[i,j,k+1]= coef*0.5*(vy[i,j,k+2]+my[i,j,k])
        mz[i,j,k-1]= coef*0.5*(vz[i,j,k-2]+mz[i,j,k])
        mz[i,j,k+1]= coef*0.5*(vz[i,j,k+2]+mz[i,j,k])
    else:
      for l in range(Ncapteur):
        i = int(data_cap[l,6])
        j = int(data_cap[l,7])
        k = int(data_cap[l,8])
        
        delta = mx[i,j,0] - vx[i,j,0]

        mx[i   ,j ,k]= data_cap[l,3 ]*coef 
        mx[i-1,j  ,k]= coef*0.5*(vx[i-2,j  ,k]+mx[i,j,k])
        mx[i-1,j+1,k]= coef*0.5*(vx[i-2,j+2,k]+mx[i,j,k])
        mx[i-1,j-1,k]= coef*0.5*(vx[i-2,j-2,k]+mx[i,j,k])
        mx[i+1,j  ,k]= coef*0.5*(vx[i+2,j  ,k]+mx[i,j,k])
        mx[i+1,j+1,k]= coef*0.5*(vx[i+2,j+2,k]+mx[i,j,k])
        mx[i+1,j-1,k]= coef*0.5*(vx[i+2,j-2,k]+mx[i,j,k])
        mx[i  ,j-1,k]= coef*0.5*(vx[i  ,j-2,k]+mx[i,j,k])
        mx[i  ,j+1,k]= coef*0.5*(vx[i  ,j+2,k]+mx[i,j,k])

        mx[i   ,j ,k-1]= coef*0.5*(vx[i  ,j  ,k-2]+mx[i,j,k])
        mx[i-1,j  ,k-1]= coef*0.5*(vx[i-2,j  ,k-2]+mx[i,j,k])
        mx[i-1,j+1,k-1]= coef*0.5*(vx[i-2,j+2,k-2]+mx[i,j,k])
        mx[i-1,j-1,k-1]= coef*0.5*(vx[i-2,j-2,k-2]+mx[i,j,k])
        mx[i+1,j  ,k-1]= coef*0.5*(vx[i+2,j  ,k-2]+mx[i,j,k])
        mx[i+1,j+1,k-1]= coef*0.5*(vx[i+2,j+2,k-2]+mx[i,j,k])
        mx[i+1,j-1,k-1]= coef*0.5*(vx[i+2,j-2,k-2]+mx[i,j,k])
        mx[i  ,j-1,k-1]= coef*0.5*(vx[i  ,j-2,k-2]+mx[i,j,k])
        mx[i  ,j+1,k-1]= coef*0.5*(vx[i  ,j+2,k-2]+mx[i,j,k])

        mx[i   ,j ,k+1]= coef*0.5*(vx[i  ,j  ,k+2]+mx[i,j,k])
        mx[i-1,j  ,k+1]= coef*0.5*(vx[i-2,j  ,k+2]+mx[i,j,k])
        mx[i-1,j+1,k+1]= coef*0.5*(vx[i-2,j+2,k+2]+mx[i,j,k])
        mx[i-1,j-1,k+1]= coef*0.5*(vx[i-2,j-2,k+2]+mx[i,j,k])
        mx[i+1,j  ,k+1]= coef*0.5*(vx[i+2,j  ,k+2]+mx[i,j,k])
        mx[i+1,j+1,k+1]= coef*0.5*(vx[i+2,j+2,k+2]+mx[i,j,k])
        mx[i+1,j-1,k+1]= coef*0.5*(vx[i+2,j-2,k+2]+mx[i,j,k])
        mx[i  ,j-1,k+1]= coef*0.5*(vx[i  ,j-2,k+2]+mx[i,j,k])
        mx[i  ,j+1,k+1]= coef*0.5*(vx[i  ,j+2,k+2]+mx[i,j,k])

        my[i,j,k]= data_cap[l,4 ]*coef 
        my[i-1,j  ,k]= coef*0.5*(vy[i-2,j  ,k]+my[i,j,k])
        my[i-1,j+1,k]= coef*0.5*(vy[i-2,j+2,k]+my[i,j,k])
        my[i-1,j-1,k]= coef*0.5*(vy[i-2,j-2,k]+my[i,j,k])
        my[i+1,j  ,k]= coef*0.5*(vy[i+2,j  ,k]+my[i,j,k])
        my[i+1,j+1,k]= coef*0.5*(vy[i+2,j+2,k]+my[i,j,k])
        my[i+1,j-1,k]= coef*0.5*(vy[i+2,j-2,k]+my[i,j,k])
        my[i  ,j-1,k]= coef*0.5*(vy[i  ,j-2,k]+my[i,j,k])
        my[i  ,j+1,k]= coef*0.5*(vy[i  ,j+2,k]+my[i,j,k])

        my[i   ,j ,k-1]= coef*0.5*(vy[i  ,j  ,k-2]+my[i,j,k])
        my[i-1,j  ,k-1]= coef*0.5*(vy[i-2,j  ,k-2]+my[i,j,k])
        my[i-1,j+1,k-1]= coef*0.5*(vy[i-2,j+2,k-2]+my[i,j,k])
        my[i-1,j-1,k-1]= coef*0.5*(vy[i-2,j-2,k-2]+my[i,j,k])
        my[i+1,j  ,k-1]= coef*0.5*(vy[i+2,j  ,k-2]+my[i,j,k])
        my[i+1,j+1,k-1]= coef*0.5*(vy[i+2,j+2,k-2]+my[i,j,k])
        my[i+1,j-1,k-1]= coef*0.5*(vy[i+2,j-2,k-2]+my[i,j,k])
        my[i  ,j-1,k-1]= coef*0.5*(vy[i  ,j-2,k-2]+my[i,j,k])
        my[i  ,j+1,k-1]= coef*0.5*(vy[i  ,j+2,k-2]+my[i,j,k])

        my[i   ,j ,k+1]= coef*0.5*(vy[i  ,j  ,k+2]+my[i,j,k])
        my[i-1,j  ,k+1]= coef*0.5*(vy[i-2,j  ,k+2]+my[i,j,k])
        my[i-1,j+1,k+1]= coef*0.5*(vy[i-2,j+2,k+2]+my[i,j,k])
        my[i-1,j-1,k+1]= coef*0.5*(vy[i-2,j-2,k+2]+my[i,j,k])
        my[i+1,j  ,k+1]= coef*0.5*(vy[i+2,j  ,k+2]+my[i,j,k])
        my[i+1,j+1,k+1]= coef*0.5*(vy[i+2,j+2,k+2]+my[i,j,k])
        my[i+1,j-1,k+1]= coef*0.5*(vy[i+2,j-2,k+2]+my[i,j,k])
        my[i  ,j-1,k+1]= coef*0.5*(vy[i  ,j-2,k+2]+my[i,j,k])
        my[i  ,j+1,k+1]= coef*0.5*(vy[i  ,j+2,k+2]+my[i,j,k])

        mz[i,j,k]= data_cap[l,5 ]*coef 
        mz[i-1,j  ,k]= coef*0.5*(vz[i-2,j  ,k]+mz[i,j,k])
        mz[i-1,j+1,k]= coef*0.5*(vz[i-2,j+2,k]+mz[i,j,k])
        mz[i-1,j-1,k]= coef*0.5*(vz[i-2,j-2,k]+mz[i,j,k])
        mz[i+1,j  ,k]= coef*0.5*(vz[i+2,j  ,k]+mz[i,j,k])
        mz[i+1,j+1,k]= coef*0.5*(vz[i+2,j+2,k]+mz[i,j,k])
        mz[i+1,j-1,k]= coef*0.5*(vz[i+2,j-2,k]+mz[i,j,k])
        mz[i  ,j-1,k]= coef*0.5*(vz[i  ,j-2,k]+mz[i,j,k])
        mz[i  ,j+1,k]= coef*0.5*(vz[i  ,j+2,k]+mz[i,j,k])

        mz[i   ,j ,k-1]= coef*0.5*(vz[i  ,j  ,k-2]+mz[i,j,k])
        mz[i-1,j  ,k-1]= coef*0.5*(vz[i-2,j  ,k-2]+mz[i,j,k])
        mz[i-1,j+1,k-1]= coef*0.5*(vz[i-2,j+2,k-2]+mz[i,j,k])
        mz[i-1,j-1,k-1]= coef*0.5*(vz[i-2,j-2,k-2]+mz[i,j,k])
        mz[i+1,j  ,k-1]= coef*0.5*(vz[i+2,j  ,k-2]+mz[i,j,k])
        mz[i+1,j+1,k-1]= coef*0.5*(vz[i+2,j+2,k-2]+mz[i,j,k])
        mz[i+1,j-1,k-1]= coef*0.5*(vz[i+2,j-2,k-2]+mz[i,j,k])
        mz[i  ,j-1,k-1]= coef*0.5*(vz[i  ,j-2,k-2]+mz[i,j,k])
        mz[i  ,j+1,k-1]= coef*0.5*(vz[i  ,j+2,k-2]+mz[i,j,k])

        mz[i   ,j ,k+1]= coef*0.5*(vz[i  ,j  ,k+2]+mz[i,j,k])
        mz[i-1,j  ,k+1]= coef*0.5*(vz[i-2,j  ,k+2]+mz[i,j,k])
        mz[i-1,j+1,k+1]= coef*0.5*(vz[i-2,j+2,k+2]+mz[i,j,k])
        mz[i-1,j-1,k+1]= coef*0.5*(vz[i-2,j-2,k+2]+mz[i,j,k])
        mz[i+1,j  ,k+1]= coef*0.5*(vz[i+2,j  ,k+2]+mz[i,j,k])
        mz[i+1,j+1,k+1]= coef*0.5*(vz[i+2,j+2,k+2]+mz[i,j,k])
        mz[i+1,j-1,k+1]= coef*0.5*(vz[i+2,j-2,k+2]+mz[i,j,k])
        mz[i  ,j-1,k+1]= coef*0.5*(vz[i  ,j-2,k+2]+mz[i,j,k])
        mz[i  ,j+1,k+1]= coef*0.5*(vz[i  ,j+2,k+2]+mz[i,j,k])

#activation spatiale terme source à l'emplacememnt des capteurs
def _nudgingLocalisation(t, capteurs, mode='dirac', tc=None):

  if tc == None:
    amor=0.66
    for key in capteurs:
       z = Internal.getNodeFromName(t,key)
       cellN = Internal.getNodeFromName(z,'cellN_src')[1]
       sh    = numpy.shape(cellN)
       capts= capteurs[key]
       if mode =='dirac':
         for cap in capts:
            l = cap[0]
            print(key, l, cap[1])
            j= l//sh[0]
            i= l-j* sh[0]
            cellN[i,j,0]=1.
       elif mode=='gauss1':
         for cap in capts:
            l = cap[0]
            j= l//sh[0]
            i= l-j* sh[0]
            cellN[i  ,j  ,0]=1.
            cellN[i-1,j  ,0]=amor
            cellN[i+1,j  ,0]=amor
            cellN[i  ,j+1,0]=amor
            cellN[i  ,j-1,0]=amor

       elif mode=='gauss3':
         alpha=0.250
         for cap in capts:
            l = cap[0]
            j= l//sh[0]
            i= l-j* sh[0]

            cellN[i-2:i+3  ,j-2:j+3  ,0]=1.

       else:
         for cap in capts:
            l = cap[0]
            j= l//sh[0]
            i= l-j* sh[0]
            cellN[i  ,j  ,0]=1.
            cellN[i-1,j  ,0]=amor
            cellN[i-1,j+1,0]=amor
            cellN[i-1,j-1,0]=amor
            cellN[i+1,j  ,0]=amor
            cellN[i+1,j+1,0]=amor
            cellN[i+1,j-1,0]=amor
            cellN[i  ,j+1,0]=amor
            cellN[i  ,j-1,0]=amor
  else:

    for z in Internal.getZones(tc):

      subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
      for s in subRegions:
         zRname = Internal.getValue(s)
         zt = Internal.getNodeFromName(t,zRname)
         sol = Internal.getNodeFromName(zt,"FlowSolution#Centers")
         print(zRname)
         cellN_src = Internal.getNodeFromName(sol,"cellN_src")[1]
         pointlistD    =  Internal.getNodeFromName1(s, 'PointListDonor')[1]
         sh  = numpy.shape(cellN_src)
         for l in range(numpy.size(pointlistD)):
            ind = pointlistD[l]
            j= ind//sh[0]
            i= ind-j* sh[0]
            cellN_src[i,j,0 ]=1

def SensorLocalisation(t, tc, xyz_cap):

  bbox={}
  for z in Internal.getZones(tc):
    bbox[ z[0] ] = G.bbox(z)
    #print(  bbox[ z[0] ])

  sh_capt = numpy.shape(xyz_cap)
  NbCapt  = sh_capt[0]

  zone_tg={}
  for c in range(NbCapt):
    for z in Internal.getZones(tc):
      z_t = Internal.getNodeFromName(t,z[0])
      sol = Internal.getNodeFromName1(z_t,'FlowSolution#Centers')
      #cellN = Internal.getNodeFromName1(sol,'cellN')[1]
      dens = Internal.getNodeFromName1(sol,'Density')[1]
      sh = numpy.shape(dens)
      #print(sh)
  
      if xyz_cap[c,0] > bbox[ z[0] ][0] and xyz_cap[c,0] < bbox[ z[0] ][3] and xyz_cap[c,1] > bbox[ z[0] ][1] and xyz_cap[c,1] < bbox[ z[0] ][4]:
        inds = D.getNearestPointIndex(z, (xyz_cap[c,0] ,xyz_cap[c,1],0.005) )
        j= inds[0]//sh[0]
        i= inds[0]-j* sh[0]
        #if cellN[i,j,0] ==1 :
        if i > 1 and i < sh[0]-2 and  j > 1 and j < sh[1]-2:
          if z[0] in zone_tg:
            zone_tg[z[0]].append( [inds[0], c ])
          else:
            zone_tg[z[0]]=[ [inds[0], c] ]

          #print('capteur:', c, 'Zone:', z[0], 'ij',  i,j, 'cellN:', cellN[i,j,0] )
          print('capteur:', c, 'Zone:', z[0], 'ij',  i,j  )

  return zone_tg

def filter(t,capteurs , nstep):

   #tmp_i = numpy.empty( (6,3), dtype=float)
   #tmp_j = numpy.empty( (3,6), dtype=float)
   tmp_i = numpy.empty( (8,5), dtype=float)
   tmp_j = numpy.empty( (5,8), dtype=float)
   for key in capteurs:
      z   = Internal.getNodeFromName(t,key)
      sol =  Internal.getNodeFromName1(z,'FlowSolution#Centers')

      level=''
      if nstep !=1: level='_P1'
      for v in ['VelocityX','VelocityY','Density','Temperature']:
        var= Internal.getNodeFromName1(sol,v+level)[1]
        sh  = numpy.shape(var)
        capts= capteurs[key]

        for cap in capts:
          ind = cap[0]; Nocapt = cap[1]
          j= ind//sh[0]
          i= ind-j* sh[0]
          #tmp_i[:,:]= 0.5*(var[i-3:i+3,j-1:j+2,0]+var[i-2:i+4,j-1:j+2,0])
          #var[i-2:i+3,j-1:j+2,0] = 0.5*(tmp_i[0:5,:]+tmp_i[1:6,:])
          #tmp_i[:,:]= 0.5*(var[i-4:i+4,j-2:j+3,0]+var[i-3:i+5,j-2:j+3,0])
          #print(Nocapt, i, j, numpy.shape(tmp_i), sh , len(capts) , flush=True)

          #if i >=4 and i+3 <= sh[0] and j >= 2 and j+2 < sh[1]:
          #  #tmp_i[:,:]= var[i-4:i+4,j-2:j+3,0]
          #  #var[i-3:i+4,j-2:j+3,0] = 0.5*(tmp_i[0:7,:]+tmp_i[1:8,:])
          #  #tmp_i[:,0:3]= var[i-4:i+4,j-1:j+2,0]
          #  #var[i-3:i+4,j-1:j+2,0] = 0.5*(tmp_i[0:7,0:3]+tmp_i[1:8,0:3])
          #  tmp_i[0:6,0:5]= var[i-3:i+3,j-2:j+3,0]
          #  var[i-2:i+3,j-2:j+3,0] = 0.5*(tmp_i[0:5,0:5]+tmp_i[1:6,0:5])
          #else:
          tmp_i[0:6,0:3]= 0.5*(var[i-3:i+3,j-1:j+2,0]+var[i-2:i+4,j-1:j+2,0])
          #var[i-2:i+3,j-1:j+2,0] = 0.5*(tmp_i[0:5,0:3]+tmp_i[1:6,0:3])
          var[i-2:i  ,j-1:j+2,0] = 0.5*(tmp_i[0:2,0:3]+tmp_i[1:3,0:3])
          var[i+1:i+3,j-1:j+2,0] = 0.5*(tmp_i[3:5,0:3]+tmp_i[4:6,0:3])

          #if j >=4 and j+3 <= sh[1] and i >= 2 and i+2 < sh[0]:
          #   #tmp_j[:,:]= 0.5*(var[i-2:i+3, j-4:j+4,0]+var[i-2:i+3, j-3:j+5,0])
          #   #var[i-2:i+3, j-3:j+4,0] = 0.5*(tmp_j[:, 0:7]+tmp_j[:,1:8])
          #   #tmp_j[0:3,:]= 0.5*(var[i-1:i+2, j-4:j+4,0]+var[i-1:i+2, j-3:j+5,0])
          #   #var[i-1:i+2, j-3:j+4,0] = 0.5*(tmp_j[0:3, 0:7]+tmp_j[0:3,1:8])
          #   tmp_j[0:5,0:6]= 0.5*(var[i-2:i+3, j-3:j+3,0]+var[i-2:i+3, j-2:j+4,0])
          #   var[i-2:i+3, j-2:j+3,0] = 0.5*(tmp_j[0:5, 0:5]+tmp_j[0:5,1:6])
          #else:
          tmp_j[0:3,0:6]= 0.5*(var[i-1:i+2, j-3:j+3,0]+var[i-1:i+2, j-2:j+4,0])
          #var[i-1:i+2, j-2:j+3,0] = 0.5*(tmp_j[0:3, 0:5]+tmp_j[0:3,1:6])
          var[i-1:i+2, j-2:j  ,0] = 0.5*(tmp_j[0:3, 0:2]+tmp_j[0:3,1:3])
          var[i-1:i+2, j+1:j+3,0] = 0.5*(tmp_j[0:3, 3:5]+tmp_j[0:3,4:6])

def filterSquareDuct(t,capteurs , nstep):

   Ncapteur = numpy.shape(capteurs)[0]
   sol =  Internal.getNodeFromName(t,'FlowSolution#Centers')

   tmp_i = numpy.empty( (6,3,3), dtype=float)
   tmp_j = numpy.empty( (3,6,3), dtype=float)
   tmp_k = numpy.empty( (3,3,6), dtype=float)
   level=''
   if nstep !=1: level='_P1'
   for v in ['VelocityX','VelocityY','VelocityZ','Density','Temperature']:
     var= Internal.getNodeFromName1(sol,v+level)[1]
     sh  = numpy.shape(var)

     for l in range(Ncapteur):
        i = int( capteurs[l,6] )
        j = int( capteurs[l,7] )
        k = int( capteurs[l,8] )

        tmp_i[0:6,0:3, 0:3]= 0.5*(var[i-3:i+3, j-1:j+2, k-1:k+2]+var[i-2:i+4, j-1:j+2, k-1:k+2])
        var[i-2:i  ,j-1:j+2,k-1:k+2] = 0.5*(tmp_i[0:2,0:3, 0:3]+tmp_i[1:3,0:3, 0:3])
        var[i+1:i+3,j-1:j+2,k-1:k+2] = 0.5*(tmp_i[3:5,0:3, 0:3]+tmp_i[4:6,0:3, 0:3])

        tmp_j[0:3,0:6, 0:3]= 0.5*(var[i-1:i+2, j-3:j+3, k-1:k+2]+var[i-1:i+2, j-2:j+4, k-1:k+2])
        var[i-1:i+2, j-2:j  ,k-1:k+2] = 0.5*(tmp_j[0:3, 0:2, 0:3]+tmp_j[0:3,1:3, 0:3])
        var[i-1:i+2, j+1:j+3,k-1:k+2] = 0.5*(tmp_j[0:3, 3:5, 0:3]+tmp_j[0:3,4:6, 0:3])

        if k> 2:
          #print(i,j,k, sh)
          tmp_k[0:3, 0:3, 0:6]= 0.5*(var[i-1:i+2, j-1:j+2, k-3:k+3]+var[i-1:i+2, j-1:j+2, k-2:k+4])
          var[i-1:i+2, j-1:j+2, k-2:k  ] = 0.5*(tmp_k[0:3,  0:3, 0:2]+tmp_k[0:3, 0:3, 1:3])
          var[i-1:i+2, j-1:j+2, k+1:k+3] = 0.5*(tmp_k[0:3,  0:3, 3:5]+tmp_k[0:3, 0:3, 4:6])


def errorEstimate(t, zone_tg, data, it_loc, dt_scale):
   #calcul erreur
   ex= 0.
   ey= 0.
   Nbcapt=0
   for key in zone_tg:
      #print(key)
      z = Internal.getNodeFromName(t,key)
      sol =  Internal.getNodeFromName1(z,'FlowSolution#Centers')
      Vx = Internal.getNodeFromName(sol,'VelocityX')[1]
      Vy = Internal.getNodeFromName(sol,'VelocityY')[1]

      sh    = numpy.shape(Vx)
      capts= zone_tg[key]
      Nbcapt +=len(capts)
      coef = float((it_loc%dt_scale))/float(dt_scale)
      t1 = it_loc//dt_scale
      t2 = t1+1
      for cap in capts:
        ind = cap[0]; Nocapt = cap[1]
        j= ind//sh[0]
        i= ind-j* sh[0]
        mx= data[Nocapt, t1,0 ] + (data[Nocapt, t2 ,0 ]-data[Nocapt, t1,0 ])*coef
        my= data[Nocapt, t1,1 ] + (data[Nocapt, t2 ,1 ]-data[Nocapt, t1,1 ])*coef

        ex += (mx-Vx[i,j,0])**2
        ey += (my-Vy[i,j,0])**2

   ex =  math.sqrt(ex/float(Nbcapt)); ey =  math.sqrt(ey/float(Nbcapt))
   #print('erreur:', ex, ey)
   return ex, ey


#modif terme source pour qu'il soit a divergence nulle
def divFree(t, nstep=1):

  level=''
  if nstep!=1:level='_P1'
