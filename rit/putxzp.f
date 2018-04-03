c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/rit/putxzp.f $
c $LastChangedDate: 2007-11-12 13:19:01 +0100 (Mon, 12 Nov 2007) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 850 $
c
c ***********************************************************************
      subroutine putxzp(pln,yp,i,ur)
c
c     Put an xy plane into ur
c
      implicit none

      include 'par.f'

      integer yp,i
      complex pln(nx/2+1,nz)
      complex ur(memnx,memny,memnz,7)
c
      integer x,z
c
c     Put plane in core
c
      do z=1,nzc
         do x=1,nx/2
            ur(x,yp,z,i)=pln(x,z)
         end do
      end do

      return

      end subroutine putxzp
