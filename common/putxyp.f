c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/common/putxyp.f $
c $LastChangedDate: 2007-11-19 11:33:23 +0100 (Mon, 19 Nov 2007) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 882 $
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

      integer x,y

      do y=1,nyp
         do x=1,nx/2
            ur(x,y,z,i)=pln(x,y)
         end do
      end do

      return

      end subroutine putxyp
