c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/rps/fofx.f $
c $LastChangedDate: 2006-11-19 23:43:46 +0100 (Sun, 19 Nov 2006) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 342 $
c
c ***********************************************************************
      subroutine fofx(plxy,plxz,gridxy,gridxz,y,z,tpl,fx,xc)
c
c     Get f of x data for one y or z
c
      implicit none

      include 'par.f'

      integer y,z,tpl
      real fx(nx+1),xc(nx+1)
      real plxy(nx+1,nyp),gridxy(nx+1,nyp,2)
      real plxz(nx+1,nz+1),gridxz(nx+1,nz+1,2)

      integer x

      do x=1,nx+1
         if (tpl.eq.1) then
            fx(x)=plxy(x,y)
            xc(x)=gridxy(x,y,1)
         else
            fx(x)=plxz(x,z)
            xc(x)=gridxz(x,z,1)
         end if
      end do

      return

      end subroutine fofx
