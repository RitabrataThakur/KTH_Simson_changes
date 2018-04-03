c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/rit/zwall.f $
c $LastChangedDate: 2006-09-28 17:40:52 +0200 (Thu, 28 Sep 2006) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 205 $
c
c ***********************************************************************
      subroutine zwall(ur,pxy)

      implicit none

      include 'par.f'

      complex ur(memnx,memny,memnz,3)
      complex pxy(nx/2,nyp)
c
      integer y

      call getxyp(pxy,1,1,ur)
      do y=1,nyp
         pxy(1,y)=pxy(1,y)-pxy(1,nyp)
      end do
      call putxyp(pxy,1,1,ur)

      return

      end subroutine zwall
