c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/cwallbc.f $
c $LastChangedDate: 2007-10-29 14:23:13 +0100 (Mon, 29 Oct 2007) $
c $LastChangedBy: qiang@MECH.KTH.SE $
c $LastChangedRevision: 786 $
c
c ***********************************************************************
      subroutine cwallbc(wallvr,wallvi,wp1,wp2,wpds1,
     &     wpds2,wpds3,wpds4,wpds5,wpds6,wpds7,wpds8,wpds9,wpds10,
     &     wpds11,wpds12,wpdds1,wpdds2,
     &     prex,prez,pres,wr,wi,xl,xsc,zl,zsc,tc,wbci)
c      
c     wbci=1:
c     Setting the boundary condition at the wall of v to
c     v(y=0) = amp*f(x)*cos(zbet*zc)*sin(tomeg*tc)
c     
c     f(x) = step((x-xstart)/xrise)-step((x-xend)/xfall+1)
c
      implicit none

      include 'par.f'

      real wallvr(nxp/2+1,nzd),wallvi(nxp/2+1,nzd)
      real wp1,wp2,wpds1,wpds2,wpds3,wpds4,wpds5
      real wpds6,wpds7,wpds8,wpds9
      real wpds10,wpds11,wpds12
      real wpdds1,wpdds2
      real prex(nxp+15),prez(nzp*2+15),pres(nzst*2+15)
      real wr(nxp/2+1,nzd),wi(nxp/2+1,nzd)
      real xl,xsc,zl,zsc,tc
      integer wbci

      integer x,z
      real zc,tmpr,tmpi
      real xc11,xc22,fx1(nxp/2),fx2(nxp/2),ft
      real amp,damp,xstart,xend,xrise,xfall,zbet,tomeg
      real zstart,zend,zrise,zfall,tstart,tend,trise,tfall
      real step
      real pi
      real dist,distom
      parameter (pi = 3.1415926535897932385)

      if (wbci.eq.1) then
         amp=wp1
         damp=wp2
         xstart=wpds1
         xend=wpds2
         xrise=wpds3
         xfall=wpds4
         zbet=wpdds1
         tomeg=wpdds2
        
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            fx1(x)=step((xc11-xstart)/xrise)-step((xc11-xend)/xfall+1)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            fx2(x)=step((xc22-xstart)/xrise)-step((xc22-xend)/xfall+1)
            do z=1,nzpc
               zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
               tmpr=amp*fx1(x)*cos(zbet*zc)*sin(tomeg*tc)
               tmpi=amp*fx2(x)*cos(zbet*zc)*sin(tomeg*tc)
               wallvr(x,z)=max(tmpr,tmpr*damp)
               wallvi(x,z)=max(tmpi,tmpi*damp)
            end do
         end do 
      end if
      if (wbci.eq.2) then
         amp=wp1
         xstart=wpds1
         xend=wpds2
         xrise=wpds3
         xfall=wpds4
         zstart=wpds5
         zend=wpds6
         zrise=wpds7
         zfall=wpds8
         tstart=wpds9
         tend=wpds10
         trise=wpds11
         tfall=wpds12
         distom=wpdds1
         
         ft=step((tc-tstart)/trise)-step((tc-tend)/tfall+1)
         dist=amp/100.*sin(distom*tc)
         do x=1,nxp/2
            xc11=real(2*x-1-nxp/2-1)/real(nxp)*xl+xsc
            xc11=xc11-int((xc11+xl/2.)/xl)*xl
            fx1(x)=step((xc11-xstart)/xrise)-step((xc11-xend)/xfall+1)
            xc22=real(2*x-nxp/2-1)/real(nxp)*xl+xsc
            xc22=xc22-int((xc22+xl/2.)/xl)*xl
            fx2(x)=step((xc22-xstart)/xrise)-step((xc22-xend)/xfall+1)
            do z=1,nzpc
               zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
               tmpr=step((zc-zstart)/zrise)-step((zc-zend)/zfall+1)
               wallvr(x,z)=(amp*fx1(x)*tmpr+dist)*ft
               wallvi(x,z)=(amp*fx2(x)*tmpr+dist)*ft
            end do
         end do
      end if        
c
c     Real to half complex transform in x-direction first
c
      call vrfftf(wallvr,wallvi,wr,wi,nxp,nzpc,1,nxp/2+1,prex)
c
c     Then complex transform in z direction
c
      if (nfzsym.eq.0) then
         call vcfftf(wallvr,wallvi,wr,wi,nzp,nx/2,nxp/2+1,1,prez)
      else
         call vcffts(wallvr,wallvi,wr,nzst,nx/2,nxp/2+1,1,pres)
      end if
c
c     Normalize
c
      do x=1,nxp/2+1 
         do z=1,nz/2
            wallvr(x,z)=wallvr(x,z)/(nxp*nzp)
            wallvi(x,z)=wallvi(x,z)/(nxp*nzp)
         end do
      end do
      if (nfzsym.eq.0) then
         do x=1,nxp/2+1
            do z=nz/2+1,nz
               wallvr(x,z)=wallvr(x,z+nzp-nz)/real(nxp*nzp)
               wallvi(x,z)=wallvi(x,z+nzp-nz)/real(nxp*nzp)
            end do
         end do
      end if

      end subroutine cwallbc
      
        
