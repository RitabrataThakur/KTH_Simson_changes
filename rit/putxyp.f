c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/rit/putxyp.f $
c $LastChangedDate: 2007-11-12 13:19:01 +0100 (Mon, 12 Nov 2007) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 850 $
c
c ***********************************************************************
      subroutine putxyp(pln,z,i,ur)
c
c     Put an xy box to ur
c
      implicit none

      include 'par.f'

      integer z,i
      complex pln(nx/2,nyp)
      complex ur(memnx,memny,memnz,7)
c
      integer x,y
c
      do y=1,nyp
         do x=1,nx/2
            ur(x,y,z,i)=pln(x,y)
         end do
      end do

      return

      end subroutine putxyp
