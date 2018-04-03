c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/pxyst/mxpl.f $
c $LastChangedDate: 2006-09-29 14:54:06 +0200 (Fri, 29 Sep 2006) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 216 $
c
c ***********************************************************************
      subroutine mxpl(xgraph,xgrid,plane,grid,y1,npy,mx,my,mgr)
c
c     Selects data for multiple function of x plots
c
      implicit none

      integer mx,my,mgr,npy,y1(mgr)
      real plane(mx,my),grid(mx,my,2),xgraph(mx,mgr),xgrid(mx,mgr)
c
      integer i,x

      do i=1,npy
         do x=1,mx
            xgraph(x,i)=plane(x,y1(i))
            xgrid(x,i)=grid(x,y1(i),1)
         end do
      end do

      return

      end subroutine mxpl
