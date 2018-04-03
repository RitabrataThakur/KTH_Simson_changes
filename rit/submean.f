c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/rit/submean.f $
c $LastChangedDate: 2007-11-12 13:19:01 +0100 (Mon, 12 Nov 2007) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 850 $
c
c ***********************************************************************
      subroutine submean(ur,meanfr)
c
c     Version of add from cmp used with rit
c
c     Adds field n to field 1 with coefficient c
c     If n=1 simply multiply by c
c
      implicit none

      include 'par.f'
c
c     Global variables
c
      complex ur(memnx,memny,memnz,5)
      complex meanfr(memnx,memny,memnz,5)
c
c     Local variables
c
      complex boxr(nx/2,nyp)
      complex box2r(nx/2,nyp)
      integer x,y,i,zb

      do i=1,3
         do zb=1,nzc
            call getxyp(boxr,zb,i,ur)
            call getxyp(box2r,zb,i,meanfr)
c
c     Loop over y-points
c
            do y=1,nyp
               if (zb.eq.1.and.i.eq.1) then
                  boxr(1,y)=boxr(1,y)-(box2r(1,y))
                  do x=2,nx/2
                     boxr(x,y)=boxr(x,y)-box2r(x,y)
                  end do
               else
                  do x=1,nx/2
                     boxr(x,y)=boxr(x,y)-box2r(x,y)
                  end do
               end if
            end do
            call putxyp(boxr,zb,i,ur)
         end do
      end do

      return

      end subroutine submean
