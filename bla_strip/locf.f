c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/bla/locf.f $
c $LastChangedDate: 2013-01-28 21:56:18 +0100 (Mon, 28 Jan 2013) $
c $LastChangedBy: pschlatt@MECH.KTH.SE $
c $LastChangedRevision: 1813 $
c
c ***********************************************************************
      subroutine locf(om2r,om2i,yb,xl,zl,xsc,zsc,eta,tc,
     &     loctyp,fp1,fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7,fpds8,
     &     fpdds1,fpdds2,fpdds3,fpdds4,fpdds5,g1,g2,
     &     th2r,th2i,u2r,u2i)

c
c     Localized forcing
c
c==== loctyp=1:
c
c     Adding a volume force of the form
c     F=(ampx,ampy,ampz)*exp(-(y/yscale)**2)*g(x,z)*f(t)
c
c     zscale>0   g(x,z)=exp(-(x-xloc0)/xscale**2-(z/zscale)**2)
c     zscale<0   g(x,z)=exp(-(x-xloc0)/xscale**2)*cos((z-x*lskew)/zscale*2*pi)
c
c     tscale>0 f(t) is a smooth turn on   : f(t)=exp(-(t/tscale)**2)
c     tscale<0 f(t) is a smooth turn off  : f(t)=step(-t/tscale))
c     tscale=0 f(t)=1.
c
c     where step is defined in step.f
c
c     the volume force is only calculated if locfor is true
c     and the time is in the interval [0-5 tscale] or tscale<0
c
c
c==== loctyp=2:
c
c     Adding a volume force of the form
c     F=(ampx,ampy,ampz).*(g(1,z't),g(2,z't),g(3,z't))**fy(y),fx(x)
c
c     g(1,z't)=cos(zbet*z)*cos(tomeg*t)/(2*tomeg)
c     g(2,z't)=cos(zbet*z)*sin(tomeg*t)
c     g(3,z't)=-sin(zbet*z)*sin(tomeg*t)/(2*zbet)
c
c     fx(x)=step((x-xstart)/xrise)-step((x-xend)/xfall+1)
c
c     fy(y)=step((y-ystart)/yrise)-step((y-yend)/yfall+1)
c
c     where step is defined in step.f
c
c==== loctyp=3:
c
c     Adding a volume force of the form
c     F=(ampx,ampy,ampz)*exp(-(y/yscale)**2)*fx(x)*f(t)
c
c     xtype=0  f(x)=                 exp(-((x-xloc0)/xscale)**2)
c     xtype=1  f(x)=(x-xloc0)/xscale*exp(-((x-xloc0)/xscale)**2)
c
c     f(t)=sin(tomeg*t)
c
c==== loctyp=4:
c
c     Adding a volume force of the form
c     F=(ampx,ampy,ampz)*exp(-y/yscale)*ft(t)
c
c     f(t)=(step((t-tstart)/tscale))-step((t-tend)/tscale+1))*cos(tomeg*t)
c
c==== loctyp=5:
c
c     Related to loctyp=1
c
c     adding a volume force of the form
c     F=(ampx,ampy,ampz)*exp(-(y/yscale)**2)*g(x,z)*f(t)
c
c
c     zscale>0   g(x,z)=exp(-(x-xloc0)/xscale**2-(z/zscale)**2)
c     zscale<0   g(x,z)=cos((z-x*lskew)/zscale*2*pi)*exp(-(x-xloc0)/xscale**2)
c
c     tscale>0 f(t) is a smooth turn on   : f(t)=exp(-(t/tscale)**2)
c     tscale<0 f(t) is a smooth turn off  : f(t)=step(-t/tscale))
c     tscale=0 f(t)=1.
c
c     h1(t)=aomeg*sin(tomeg*tc) (note: aomeg is relative amplitude
c     between oscillation and stationary force)
c
c     where step is defined in step.f
c
c     the volume force is only calculated if locfor is true
c     and the time is in the interval [0-5 tscale] or tscale<0
c
c==== loctyp=6:
c
c     Adding a volume force of the form
c     F=(ampx,ampy,ampz)*exp(-(y/yscale)**2)*fx(x)*f(t)
c     useful for TS waves and corresponding secondary instability (K/H-type)
c
c     f(x)=exp(-((x-xloc0)/xscale)**2)
c
c     g(z) = cos(2pi/zl)
c
c     f2d(t)   = sin(tomeg*t)
c     f3d(t) = sin(tomeg3D*t)
c
c     F=(0,1,0)*exp(-(yc/yscale)**2)*(amp2d*f2d(t)*f(x) +
c                                     amp3d*f3d(t)*f(x)*g(z) )
c
c==== loctyp=7:
c
c     Adding a localised forcing of the temperature (line source)
c
c==== loctyp=8:
c
c     Approximated an impulse input
c
c     Adding a volume force of the form
c     F=(ampx,ampy,ampz)*exp(-((y-yloc)/yscale)**2)*fx(x)*f(t)
c
c     xtype=0  f(x)=                 exp(-((x-xloc0)/xscale)**2)
c
c     f(t)=exp(-((t-tstart)/tscale)**2)
c

      implicit none

      include 'par.f'

      integer yb,loctyp,ith
      real om2r(nxp/2+1,mby,nzd,3),om2i(nxp/2+1,mby,nzd,3)
      real th2r(nxp/2+1,mby,nzd,4*scalar),th2i(nxp/2+1,mby,nzd,4*scalar)
      real u2r(nxp/2+1,mby,nzd,3),u2i(nxp/2+1,mby,nzd,3)
      real xl,zl,xsc,zsc,tc
      real eta(nyp)
      real g1(nxp/2,nzd),g2(nxp/2,nzd)
      real fp1,fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7
      real fpds8,fpdds1,fpdds2,fpdds3,fpdds4,fpdds5

      real ampx,ampy,ampz,xscale,yscale,zscale,tscale,xloc0,lskew
      real xtype,xstart,xend,xrise,xfall,ystart,yend,yrise,yfall
      real zbet,tomeg,tstart,tend,xo,aomeg,y0,yscale1
      real amp,yloc0

      real xc1(nxp/2),xc2(nxp/2)
      integer x,y,z
      real xc11,xc22,fx1(nxp/2),fx2(nxp/2),f2x1(nxp/2),f2x2(nxp/2)
      real fy,dfy,f2y,ft,f2t,yc,zc,k2
      real pi
      real amp2d,amp3d,tomeg3d,ft3d
      real rad,lam,lam0
      parameter (pi = 3.1415926535897932385)
c
c     Functions
c
      real,external :: step,dstep

c
c Note that we need a coordinate that runs between -xl/2 and xl/2
c regardless of the shift xsc, otherwise the force would be turned off
c abruptly when shifted out of the box
c
      if (loctyp.eq.1) then
c
c     Rename variables to simplify understanding
c
         ampx=fpdds1
         ampy=fpdds2
         ampz=fpdds3
         xloc0=fpds1
         xscale=fpds2
         yscale=fpds3
         zscale=fpds4
         lskew=fp1
         tscale=fpds5

         if (tscale.gt.0..and.tc.gt.5.*tscale) return

         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x)=exp(-((xc11-xloc0)/xscale)**2)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x)=exp(-((xc22-xloc0)/xscale)**2)
         end do
c
c     Construct g(x,z)
c
         if (zscale.gt.0) then
            do z=1,nzpc
               do x=1,nxp/2
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  g1(x,z)=exp(-(zc/zscale)**2)*fx1(x)
                  g2(x,z)=exp(-(zc/zscale)**2)*fx2(x)
               end do
            end do
         else
            do z=1,nzpc
               do x=1,nxp/2
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  g1(x,z)=cos((zc-(xc1(x)-xloc0)*lskew)/zscale*2.*pi)
     &                 *fx1(x)
                  g2(x,z)=cos((zc-(xc2(x)-xloc0)*lskew)/zscale*2.*pi)
     &                 *fx2(x)
               end do
            end do
        end if
c
c     Construct f(t)
c
        if (tscale.gt.0.) ft = exp(-(tc/tscale)**2)
        if (tscale.lt.0.) ft = step(-(tc/tscale))
        if (tscale.eq.0.) ft=1.

        do z=1,nzpc
           do y=1,min(mby,nyp-yb+1)
              yc=1.+eta(y+yb-1)
              fy=exp(-(yc/yscale)**2)
              do x=1,nxp/2
c
c     Volume force, note that it is nonzero on the wall
c
                 om2r(x,y,z,1)=om2r(x,y,z,1)+ampx*fy*ft*g1(x,z)
                 om2i(x,y,z,1)=om2i(x,y,z,1)+ampx*fy*ft*g2(x,z)
                 om2r(x,y,z,2)=om2r(x,y,z,2)+ampy*fy*ft*g1(x,z)
                 om2i(x,y,z,2)=om2i(x,y,z,2)+ampy*fy*ft*g2(x,z)
                 om2r(x,y,z,3)=om2r(x,y,z,3)+ampz*fy*ft*g1(x,z)
                 om2i(x,y,z,3)=om2i(x,y,z,3)+ampz*fy*ft*g2(x,z)
              end do
           end do
        end do

      end if

      if (loctyp.eq.2) then
c
c     Stellans note: disturbance slightly changed since wiegel and
c     turb simulations.
c
c     Rename variables to simplify understanding
c
         ampx=fpdds1
         ampy=fpdds2
         ampz=fpdds3
         xstart=fpds1
         xend=fpds2
         xrise=fpds3
         xfall=fpds4
         ystart=fpds5
         yend=fpds6
         yrise=fpds7
         yfall=fpds8
         zbet=fpdds4
         tomeg=fpdds5
         k2=sqrt(tomeg*tomeg+zbet*zbet)
c
c     Construct fx(x)
c
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x) = step((xc1(x)-xstart)/xrise)-
     &           step((xc1(x)-xend)/xfall+1)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x) = step((xc2(x)-xstart)/xrise)-
     &           step((xc2(x)-xend)/xfall+1)
         end do
c
c     Construct g(i,z)
c
         do z=1,nzpc
            zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
            g1(1,z)=-tomeg*cos(zbet*zc)*sin(tomeg*tc)/k2
            g1(2,z)=cos(zbet*zc)*cos(tomeg*tc)
            g1(3,z)=-zbet*sin(zbet*zc)*cos(tomeg*tc)/k2
         end do

         do y=1,min(mby,nyp-yb+1)
            yc=1.+eta(y+yb-1)
            fy=step((yc-ystart)/yrise)-step((yc-yend)/yfall+1)
            dfy=dstep((yc-ystart)/yrise)/yrise
     &           -dstep((yc-yend)/yfall+1)/yfall
            do z=1,nzpc
               do x=1,nxp/2
                  om2r(x,y,z,1)=om2r(x,y,z,1)+ampx*g1(1,z)*dfy*fx1(x)
                  om2i(x,y,z,1)=om2i(x,y,z,1)+ampx*g1(1,z)*dfy*fx2(x)
                  om2r(x,y,z,2)=om2r(x,y,z,2)+ampy*g1(2,z)*fy *fx1(x)
                  om2i(x,y,z,2)=om2i(x,y,z,2)+ampy*g1(2,z)*fy *fx2(x)
                  om2r(x,y,z,3)=om2r(x,y,z,3)+ampz*g1(3,z)*dfy*fx1(x)
                  om2i(x,y,z,3)=om2i(x,y,z,3)+ampz*g1(3,z)*dfy*fx2(x)
               end do
            end do
         end do
      end if

      if (loctyp.eq.3) then
c
c     Rename variables to simplify understanding
c
         ampx=fpdds1
         ampy=fpdds2
         ampz=fpdds3
         xscale=fpds1
         xloc0=fpds2
         yscale=fpds3
         xtype=fp1
         tomeg=fpdds4

c
c     Construct f(x)
c
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x)=((xc11-xloc0)/xscale)**xtype
     &           *exp(-((xc11-xloc0)/xscale)**2)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x)=((xc22-xloc0)/xscale)**xtype
     &           *exp(-((xc22-xloc0)/xscale)**2)
         end do
c
c     Construct f(t)
c
         ft=sin(tc*tomeg)

         do y=1,min(mby,nyp-yb+1)
            yc = 1.+eta(y+yb-1)
            fy = exp(-(yc/yscale)**2)
c            write(*,*) yb,fy

            do z=1,nzpc
               do x=1,nxp/2
c
c     Volume force, note that it is nonzero on the wall
c
                  om2r(x,y,z,1)=om2r(x,y,z,1)+ampx*fy*ft*fx1(x)
                  om2i(x,y,z,1)=om2i(x,y,z,1)+ampx*fy*ft*fx2(x)
                  om2r(x,y,z,2)=om2r(x,y,z,2)+ampy*fy*ft*fx1(x)
                  om2i(x,y,z,2)=om2i(x,y,z,2)+ampy*fy*ft*fx2(x)
                  om2r(x,y,z,3)=om2r(x,y,z,3)+ampz*fy*ft*fx1(x)
                  om2i(x,y,z,3)=om2i(x,y,z,3)+ampz*fy*ft*fx2(x)
               end do
            end do
         end do
      end if

      if (loctyp.eq.4) then
         ampz=fpdds1
         yscale=fpds1
         tstart=fpds2
         tend=fpds3
         tscale=fpds4
         tomeg=fpdds2
c
c     Construct f(t)
c
         ft=(step((tc-tstart)/tscale)-step((tc-tend)/tscale+1))*
     &        cos(tc*tomeg)

         do y=1,min(mby,nyp-yb+1)
            yc=1.+eta(y+yb-1)
            fy=exp(-yc/yscale)
            do z=1,nzpc
               do x=1,nxp/2
c
c     Volume force, note that it is nonzero on the wall
c
                  om2r(x,y,z,3)=om2r(x,y,z,3)+ampz*fy*ft
                  om2i(x,y,z,3)=om2i(x,y,z,3)+ampz*fy*ft
               end do
            end do
         end do
      end if

      if (loctyp.eq.5) then
c
c     Rename variables to simplify understanding
c
         ampx=fpdds1
         ampy=fpdds2
         ampz=fpdds3
         xloc0=fpds1
         xscale=fpds2
         yscale=fpds3
         zscale=fpds4
         lskew=fp1
         tscale=fpds5
         xo=fpds6
         aomeg=fpdds4
         tomeg=fpdds5
         y0=fpds7
         yscale1=fpds8

         if (tscale.gt.0..and.tc.gt.5.*tscale) return
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x)=exp(-((xc11-xloc0)/xscale)**2)
            f2x1(x)=exp(-((xc11-xo)/xscale)**2)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x)=exp(-((xc22-xloc0)/xscale)**2)
            f2x2(x)=exp(-((xc22-xo)/xscale)**2)
         end do
c
c     Construct g(x,z)
c
         if (zscale.gt.0) then
            do z=1,nzpc
               do x=1,nxp/2
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  g1(x,z)=exp(-(zc/zscale)**2)*fx1(x)
                  g2(x,z)=exp(-(zc/zscale)**2)*fx2(x)
               end do
            end do
         else
            do z=1,nzpc
               do x=1,nxp/2
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  g1(x,z)=cos((zc-(xc1(x)-xloc0)*lskew)/zscale*2.*pi)
     &                 *fx1(x)
                  g2(x,z)=cos((zc-(xc2(x)-xloc0)*lskew)/zscale*2.*pi)
     &                 *fx2(x)
               end do
            end do
         end if
c
c     Construct f(t)
c
         if (tscale.gt.0.) ft=exp(-(tc/tscale)**2)
         if (tscale.lt.0.) ft=step(-(tc/tscale))
         if (tscale.eq.0.) ft=1.
c
c     Construct f2(t)
c
         f2t=aomeg*sin(tomeg*tc)

         do z=1,nzpc
            do y=1,min(mby,nyp-yb+1)
               yc=1.+eta(y+yb-1)
               fy=exp(-(yc/yscale)**2)
               f2y=exp(-((yc-y0)/yscale1)**2)
               do x=1,nxp/2
c
c     Volume force, note that it is nonzero on the wall
c
                  om2r(x,y,z,1)=om2r(x,y,z,1)+ampx*(fy*ft*g1(x,z)
     &                 +f2y*f2t*f2x1(x))
                  om2i(x,y,z,1)=om2i(x,y,z,1)+ampx*(fy*ft*g2(x,z)
     &                 +f2y*f2t*f2x2(x))
                  om2r(x,y,z,2)=om2r(x,y,z,2)+ampy*(fy*ft*g1(x,z)
     &                 +f2y*f2t*f2x1(x))
                  om2i(x,y,z,2)=om2i(x,y,z,2)+ampy*(fy*ft*g2(x,z)
     &                 +f2y*f2t*f2x2(x))
                  om2r(x,y,z,3)=om2r(x,y,z,3)+ampz*(fy*ft*g1(x,z)
     &                 +f2y*f2t*f2x1(x))
                  om2i(x,y,z,3)=om2i(x,y,z,3)+ampz*(fy*ft*g2(x,z)
     &                 +f2y*f2t*f2x2(x))
               end do
            end do
         end do
      end if

      if (loctyp.eq.6) then
c
c     Rename variables to simplify understanding
c
         amp2d=fpdds1
         xscale=fpds1
         xloc0=fpds2
         yscale=fpds3
         tomeg=fpdds4
         amp3d=fpdds2
         tomeg3D=fpdds3
c
c     Construct f(x)
c
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x)=exp(-((xc11-xloc0)/xscale)**2)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x)=exp(-((xc22-xloc0)/xscale)**2)
         end do
c
c     Construct g(1,z)
c
         do z=1,nzpc
            zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
            g1(1,z) = cos(2*pi/zl*zc)
         end do
c
c     Construct f(t)
c
         ft  = sin(tc*tomeg)
         ft3d = sin(tc*tomeg3D)

         do y=1,min(mby,nyp-yb+1)
            yc = 1.+eta(y+yb-1)
            fy = exp(-(yc/yscale)**2)

            do z=1,nzpc
               do x=1,nxp/2
c
c     Volume force, note that it is nonzero on the wall
c
                  om2r(x,y,z,2)=om2r(x,y,z,2)+amp2d*fy*ft*fx1(x)+
     &                 amp3d*fy*ft3d*fx1(x)*g1(1,z)
                  om2i(x,y,z,2)=om2i(x,y,z,2)+amp2d*fy*ft*fx2(x)+
     &                 amp3d*fy*ft3d*fx2(x)*g1(1,z)
               end do
            end do
         end do
      end if

      if (loctyp.eq.7) then
c
c     Rename variables to simplify understanding
c
         amp=fpdds1
         xloc0=fpds1
         yloc0=fpds2
         xscale=fpds3
         yscale=fpds4
         tomeg=fpdds2
c
c     Construct f(x)
c
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x)=exp(-((xc11-xloc0)/xscale)**2)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x)=exp(-((xc22-xloc0)/xscale)**2)
         end do
c
c     Construct f(t)
c
         ft  = cos(tc*tomeg)

         do ith=1,scalar
            do y=1,min(mby,nyp-yb+1)
               yc = 1.+eta(y+yb-1)
               fy = exp(-((yc-yloc0)/yscale)**2)
               do z=1,nzpc
                  do x=1,nxp/2
                     th2r(x,y,z,4+4*(ith-1))=th2r(x,y,z,4+4*(ith-1))
     &                    +amp*fy*ft*fx1(x)
                     th2i(x,y,z,4+4*(ith-1))=th2i(x,y,z,4+4*(ith-1))
     &                    +amp*fy*ft*fx2(x)
                  end do
               end do
            end do
         end do
      end if

      if (loctyp.eq.8) then
c
c     Rename variables to simplify understanding
c
         ampx=fpdds1
         ampy=fpdds2
         ampz=fpdds3
         xscale=fpds1
         xloc0=fpds2
         yscale=fpds3
         yloc0=fpds6
         tscale=fpds4
         tstart=fpds5
c
c     Construct f(x)
c
         xtype = 0.
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            xc1(x)=xc11
            fx1(x)=((xc11-xloc0)/xscale)**xtype
     &           *exp(-((xc11-xloc0)/xscale)**2)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            xc2(x)=xc22
            fx2(x)=((xc22-xloc0)/xscale)**xtype
     &           *exp(-((xc22-xloc0)/xscale)**2)
         end do
c
c     Construct f(t)
c
         ft=exp(-((tc-tstart)/tscale)**2)

         do y=1,min(mby,nyp-yb+1)
            yc = 1.+eta(y+yb-1)
            fy = exp(-((yc-yloc0)/yscale)**2)

            do z=1,nzpc
               do x=1,nxp/2
c
c     Volume force, note that it is nonzero on the wall
c
                  om2r(x,y,z,1)=om2r(x,y,z,1)+ampx*fy*ft*fx1(x)
                  om2i(x,y,z,1)=om2i(x,y,z,1)+ampx*fy*ft*fx2(x)
                  om2r(x,y,z,2)=om2r(x,y,z,2)+ampy*fy*ft*fx1(x)
                  om2i(x,y,z,2)=om2i(x,y,z,2)+ampy*fy*ft*fx2(x)
                  om2r(x,y,z,3)=om2r(x,y,z,3)+ampz*fy*ft*fx1(x)
                  om2i(x,y,z,3)=om2i(x,y,z,3)+ampz*fy*ft*fx2(x)
               end do
            end do
         end do
      end if


      if (loctyp.eq.9) then

         do y=1,min(mby,nyp-yb+1)
            yc = eta(y+yb-1)
            do z=1,nzpc
               do x=1,nxp/2

                  lam0 = 0.8

                  rad = sqrt((real(2*x-1)/real(nxp)*xl-xl/2.)**2 +
     &                 (real(z)/real(nzpc)*zl-zl/2.)**2)

                  lam = step( (rad-170.)/20. )

                  om2r(x,y,z,1)=om2r(x,y,z,1)-lam*(u2r(x,y,z,1)-yc)
                  om2r(x,y,z,2)=om2r(x,y,z,2)-lam*(u2r(x,y,z,2))
                  om2r(x,y,z,3)=om2r(x,y,z,3)-lam*(u2r(x,y,z,3))

                  rad = sqrt((real(2*x)/real(nxp)*xl-xl/2.)**2 +
     &                 (real(z)/real(nzpc)*zl-zl/2.)**2)

                  lam = step( (rad-170.)/20. )

                  om2i(x,y,z,1)=om2i(x,y,z,1)-lam*(u2i(x,y,z,1)-yc)
                  om2i(x,y,z,2)=om2i(x,y,z,2)-lam*(u2i(x,y,z,2))
                  om2i(x,y,z,3)=om2i(x,y,z,3)-lam*(u2i(x,y,z,3))

               end do
            end do
         end do
      end if
      if (loctyp.eq.10) then
c
c     Uniform heating
c     for the moment only uniform in y
c
c     Rename variables to simplify understanding
c
         amp=fpdds1
         yloc0=fpds1
         yscale=fpds2
c
c     Construct f(t)
c
         do ith=1,scalar
            do y=1,min(mby,nyp-yb+1)
               yc = 1.+eta(y+yb-1)
c              fy = exp(-((yc-yloc0)/yscale)**2)
               fy=1.
               do z=1,nzpc
                  do x=1,nxp/2
                     th2r(x,y,z,4+4*(ith-1))=th2r(x,y,z,4+4*(ith-1))
     &                    +amp*fy
                     th2i(x,y,z,4+4*(ith-1))=th2i(x,y,z,4+4*(ith-1))
     &                    +amp*fy
                  end do
               end do
            end do
         end do
      end if



      if (loctyp.gt.10) then
         write(*,*) 'loctyp = ',loctyp,' not implemented.'
         call stopnow(5433332)
      end if

      end subroutine locf
