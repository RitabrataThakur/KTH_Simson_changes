c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/rit/getxzp.f $
c $LastChangedDate: 2006-09-28 17:40:52 +0200 (Thu, 28 Sep 2006) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 205 $
c
c ***********************************************************************
      subroutine getxzp(pln,yp,i,ur,sym)
c
c     Get an xz plane from ur
c     sym 1 for z-symmetric variable,-1 for z-antisymmetric variable
c
      implicit none

      include 'par.f'

      integer yp,i
      real sym
      complex pln(nx/2+1,nz)
      complex ur(memnx,memny,memnz,7)
c
      integer x,z
c
c     Get plane from core
c
      do z=1,nzc
         do x=1,nx/2
            pln(x,z)=ur(x,yp,z,i)
         end do
      end do
      if (nfzsym.eq.1) then
         do z=nz/2+2,nz
            do x=1,nx/2
               pln(x,z)=sym*pln(x,nz+2-z)
            end do
         end do
         do x=1,nx/2
            pln(x,nz/2+1)=0.0
         end do
      end if
      do z=1,nz
         pln(nx/2+1,z)=0.0
      end do

      return

      end subroutine getxzp
