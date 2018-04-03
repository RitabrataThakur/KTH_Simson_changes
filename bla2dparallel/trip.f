c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/trip.f $
c $LastChangedDate: 2007-10-29 14:23:13 +0100 (Mon, 29 Oct 2007) $
c $LastChangedBy: qiang@MECH.KTH.SE $
c $LastChangedRevision: 786 $
c
c ***********************************************************************
      subroutine trip(om2r,om2i,yb,xl,xsc,eta,txsc,tysc,ttzc,tx0,
     &     my_node_world)
      
      implicit none
      include 'par.f'
      integer yb
      real om2r(nxp/2+1,mby,nzd_new/nprocz,3)
      real om2i(nxp/2+1,mby,nzd_new/nprocz,3)
      real txsc,tysc,ttzc(nzp+2,4),tx0
      real xl,xsc
      real eta(nyp)

      integer x,y,z,zp,zb,my_node_world
      real yc,xc11,xc22,ex1(nxp/2),ex2(nxp/2),expy
c
c     Adding a volume force of the form
c     Fv=exp(-((x-tx0)/txsc)**2-(y/tysc)**2)*f(z,t)
c
c     Note that we need a coordinate that runs between -xl/2 and xl/2
c     regardless of the shift xsc, otherwise the force would be turned off
c     abruptly when shifted out of the box
c
c     f(z,t) is computed in gtrip.f
c
c     Compute f(x)=exp(-((x-tx0)/txsc)**2)
c
      zb=mod(my_node_world,nprocx)*nzd_new/nprocz
      do x=1,nxp/2
         xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
         xc11=xc11-int((xc11+xl/2.)/xl)*xl
         ex1(x)=exp(-((xc11-tx0)/txsc)**2)
         xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
         xc22=xc22-int((xc22+xl/2.)/xl)*xl
         ex2(x)=exp(-((xc22-tx0)/txsc)**2)
      end do
      do x=1,nxp/2
         if(yb.eq.1) then
c            write(4567,*)x,z,yb,expy,ex1(x),ex2(x),ttzc(z,1)
c            write(400+my_node_world,*)x,ex1(x),ex2(x)
         end if
         if(yb.eq.2) then
c            write(4567,*)x,z,yb,expy,ex1(x),ex2(x),ttzc(z,1)
c            write(500+my_node_world,*)x,ex1(x),ex2(x)
         end if
         if(yb.eq.5) then
c            write(4567,*)x,z,yb,expy,ex1(x),ex2(x),ttzc(z,1)
c            write(600+my_node_world,*)x,ex1(x),ex2(x)
         end if
      end do
      do zp=1,nzd_new/nprocz
         z=zp+zb
         do y=1,min(mby,nyp-yb+1)
            yc=1.+eta(y+yb-1)
            expy=exp(-(yc/tysc)**2)
            do x=1,nxp/2
c
c     Volume force, note that it is nonzero on the wall
c
               om2r(x,y,zp,2)=om2r(x,y,zp,2)+expy*ex1(x)*ttzc(z,1)
               om2i(x,y,zp,2)=om2i(x,y,zp,2)+expy*ex2(x)*ttzc(z,1)
            end do
         end do
      end do
      
      end subroutine trip
