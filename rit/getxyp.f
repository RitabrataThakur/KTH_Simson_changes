c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/rit/getxyp.f $
c $LastChangedDate: 2007-11-12 13:19:01 +0100 (Mon, 12 Nov 2007) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 850 $
c
c ***********************************************************************
      subroutine getxyp(pln,z,i,ur)
c
c     Get an xy box from ur
c
      implicit none

      include 'par.f'

      integer z,i
      complex pln(nx/2,nyp)
      complex ur(memnx,memny,memnz,7)
c
      integer x,y
      do y=1,nyp
         do x=1,nx/2
            pln(x,y)=ur(x,y,z,i)
         end do
      end do

      end subroutine getxyp
