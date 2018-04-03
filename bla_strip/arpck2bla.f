c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/bla/arpck2bla.f $
c $LastChangedDate: 2010-09-03 14:36:38 +0200 (Fri, 03 Sep 2010) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 1503 $
c
c ***********************************************************************

      subroutine arpck2bla(workd,ur,ui,maxn)
c
c     Transform workd(ipntr) to bla-format
c
      implicit none

      include 'par.f'

      integer maxn
      real ur(memnx,memny,memnz,3),ui(memnx,memny,memnz,3)
      real workd(maxn)
      integer i,ll,x,y,z,jj

      ll = 3
      do i=1,ll
         do z=1,nz/nproc
            do y = 1,ny
               do x=1,nx/2
                  jj = (i-1)*nz*ny*nx/nproc +
     &                 (z-1)*ny*nx + (y-1)*nx

                  ur(x,y,z,i)=workd(jj+2*x-1)
                  ui(x,y,z,i)=workd(jj+2*x)
               end do
            end do
         end do
      end do

      end subroutine arpck2bla
