c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/cmp/putxy.f $
c $LastChangedDate: 2007-11-12 13:19:01 +0100 (Mon, 12 Nov 2007) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 850 $
c
c ***********************************************************************
      subroutine putxy(boxr,boxi,zb,i,ur,ui)
c
c     Put an xy box to ur
c
      implicit none

      include 'par.f'

      integer zb,i
      real boxr(nx/2,mbz,nyp),boxi(nx/2,mbz,nyp)
      real ur(memnx,memny,memnz,7),ui(memnx,memny,memnz,7)
c
      integer x,y,z
c
      do z=zb,zb+mbz-1
         do y=1,nyp
            do x=1,nx/2
               ur(x,y,z,i)=boxr(x,z-zb+1,y)
               ui(x,y,z,i)=boxi(x,z-zb+1,y)
            end do
         end do
      end do

      return

      end subroutine putxy
