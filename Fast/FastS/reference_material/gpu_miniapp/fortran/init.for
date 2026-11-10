c***********************************************************************
c     Donnees synthetiques (identiques a cpp/init.hpp): champ primitif
c     regulier, metrique 3dfull legerement deformee, coefficients
c     implicites coe (dt/vol, rayons spectraux, diagonale).
c***********************************************************************
      subroutine init_fields(param_int, param_real, rop, ti,tj,tk, vol,
     &                       coe)
      implicit none
#include "mini_param.h"
      INTEGER_E param_int(0:*)
      REAL_E    param_real(0:*)
      REAL_E rop(param_int(NDIMDX),5), coe(param_int(NDIMDX),5)
      REAL_E ti(param_int(NDIMDX_MTR),3), tj(param_int(NDIMDX_MTR),3)
      REAL_E tk(param_int(NDIMDX_MTR),3), vol(param_int(NDIMDX_MTR))
      INTEGER_E i,j,k,l,lt,ni,nj,nk,inci,incj,inck
      REAL_E dx,dy,dz,xc,yc,zc,xf,yf,zf,twopi,gam,rgp,c,dt,
     & sx,sy,sz,qn,u,v,w,t
#include "formule_param.h"
#include "formule_mtr_param.h"

      ni = param_int(IJKV)
      nj = param_int(IJKV+1)
      nk = param_int(IJKV+2)
      inci = 1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)
      dx = 1.d0/ni
      dy = 1.d0/nj
      dz = 1.d0/nk
      twopi = 8.d0*atan(1.d0)
      gam = param_real(GAMMA)
      rgp = param_real(CVINF)*(gam-1.d0)
      dt  = param_real(DTC)

      do k = -1, nk+2
      do j = -1, nj+2
      do i = -1, ni+2
        l  = inddm(i,j,k)
        lt = indmtr(i,j,k)
        xc = (i-0.5d0)*dx
        yc = (j-0.5d0)*dy
        zc = (k-0.5d0)*dz
        xf = (i-1)*dx
        yf = (j-1)*dy
        zf = (k-1)*dz
        rop(l,1) = 1.d0 + 0.1d0*sin(twopi*xc)*cos(twopi*yc)*cos(twopi*zc)
        rop(l,2) = 0.5d0 + 0.1d0*cos(twopi*xc)*sin(twopi*yc)
        rop(l,3) = 0.1d0*sin(twopi*zc)*cos(twopi*xc)
        rop(l,4) = 0.05d0*cos(twopi*yc)*sin(twopi*zc)
        rop(l,5) = 1.d0 + 0.05d0*cos(twopi*xc)*cos(twopi*yc)*sin(twopi*zc)

        ti(lt,1) = dy*dz*(1.d0+0.1d0*sin(twopi*yc)*cos(twopi*zc))
        ti(lt,2) = dy*dz*0.05d0*cos(twopi*xf)
        ti(lt,3) = dy*dz*0.02d0*sin(twopi*zc)
        tj(lt,1) = dx*dz*0.05d0*sin(twopi*yf)
        tj(lt,2) = dx*dz*(1.d0+0.1d0*cos(twopi*xc)*sin(twopi*zc))
        tj(lt,3) = dx*dz*0.02d0*cos(twopi*xc)
        tk(lt,1) = dx*dy*0.02d0*cos(twopi*zf)
        tk(lt,2) = dx*dy*0.05d0*sin(twopi*xc)
        tk(lt,3) = dx*dy*(1.d0+0.1d0*sin(twopi*xc)*cos(twopi*yc))
        vol(lt)  = dx*dy*dz*(1.d0
     &           + 0.1d0*sin(twopi*xc)*sin(twopi*yc)*cos(twopi*zc))
        coe(l,1) = 0.d0
        coe(l,2) = 0.d0
        coe(l,3) = 0.d0
        coe(l,4) = 0.d0
        coe(l,5) = 1.d0
      enddo
      enddo
      enddo

c     coefficients implicites sur les cellules interieures
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
        l  = inddm(i,j,k)
        lt = indmtr(i,j,k)
        u = rop(l,2)
        v = rop(l,3)
        w = rop(l,4)
        t = rop(l,5)
        c = sqrt(gam*rgp*t)
        coe(l,1) = dt/vol(lt)
        sx = 0.5d0*(ti(lt,1)+ti(lt+inci,1))
        sy = 0.5d0*(ti(lt,2)+ti(lt+inci,2))
        sz = 0.5d0*(ti(lt,3)+ti(lt+inci,3))
        qn = u*sx+v*sy+w*sz
        coe(l,2) = 2.d0*(abs(qn) + c*sqrt(sx*sx+sy*sy+sz*sz))
        sx = 0.5d0*(tj(lt,1)+tj(lt+incj,1))
        sy = 0.5d0*(tj(lt,2)+tj(lt+incj,2))
        sz = 0.5d0*(tj(lt,3)+tj(lt+incj,3))
        qn = u*sx+v*sy+w*sz
        coe(l,3) = 2.d0*(abs(qn) + c*sqrt(sx*sx+sy*sy+sz*sz))
        sx = 0.5d0*(tk(lt,1)+tk(lt+inck,1))
        sy = 0.5d0*(tk(lt,2)+tk(lt+inck,2))
        sz = 0.5d0*(tk(lt,3)+tk(lt+inck,3))
        qn = u*sx+v*sy+w*sz
        coe(l,4) = 2.d0*(abs(qn) + c*sqrt(sx*sx+sy*sy+sz*sz))
        coe(l,5) = 1.d0 + coe(l,1)*(coe(l,2)+coe(l,3)+coe(l,4))*0.5d0
      enddo
      enddo
      enddo
      end

c***********************************************************************
c     Sommes de controle sur les cellules interieures, par variable
c***********************************************************************
      subroutine checksum(param_int, name, field)
      implicit none
#include "mini_param.h"
      INTEGER_E param_int(0:*)
      character*(*) name
      REAL_E field(param_int(NDIMDX),5)
      INTEGER_E i,j,k,l,ni,nj,nk,nv
      REAL_E s,sa,smax
#include "formule_param.h"
      ni = param_int(IJKV)
      nj = param_int(IJKV+1)
      nk = param_int(IJKV+2)
      do nv = 1, 5
        s    = 0.d0
        sa   = 0.d0
        smax = 0.d0
        do k = 1, nk
        do j = 1, nj
        do i = 1, ni
          l = inddm(i,j,k)
          s    = s  + field(l,nv)
          sa   = sa + abs(field(l,nv))
          smax = max(smax, abs(field(l,nv)))
        enddo
        enddo
        enddo
        write(*,'(a,1x,a,1x,a,i1,3(1x,a,es24.16))') 'CHECK', name,
     &    'v',nv, 'sum=',s, 'abssum=',sa, 'max=',smax
      enddo
      end

c***********************************************************************
c     Ecart max relatif entre deux champs (cellules interieures)
c***********************************************************************
      subroutine maxdiff(param_int, name, f1, f2)
      implicit none
#include "mini_param.h"
      INTEGER_E param_int(0:*)
      character*(*) name
      REAL_E f1(param_int(NDIMDX),5), f2(param_int(NDIMDX),5)
      INTEGER_E i,j,k,l,ni,nj,nk,nv
      REAL_E d,dmax,fmax
#include "formule_param.h"
      ni = param_int(IJKV)
      nj = param_int(IJKV+1)
      nk = param_int(IJKV+2)
      dmax = 0.d0
      fmax = 0.d0
      do nv = 1, 5
        do k = 1, nk
        do j = 1, nj
        do i = 1, ni
          l = inddm(i,j,k)
          d = abs(f1(l,nv)-f2(l,nv))
          dmax = max(dmax, d)
          fmax = max(fmax, abs(f1(l,nv)))
        enddo
        enddo
        enddo
      enddo
      write(*,'(a,1x,a,1x,a,es12.4,1x,a,es12.4)') 'DIFF', name,
     &  'maxabs=', dmax, 'rel=', dmax/max(fmax,1.d-300)
      end
