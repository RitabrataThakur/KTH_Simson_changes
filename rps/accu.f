c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/rps/accu.f $
c $LastChangedDate: 2007-11-12 13:19:01 +0100 (Mon, 12 Nov 2007) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 850 $
c
c ***********************************************************************
      subroutine accu(graf,grafm,np,tint,nint,tf,tl,ipln)
c
c     Take temporal average of 1-d graf
c     this routine recieves one graf at the time and accumulates
c     average
c
      implicit none

      integer nint,ipln,np
      real tint(nint),tf,tl
      real graf(np),grafm(np)

      real c
      integer i
c
c     Calculate integration weights
c
      call intwgt(c,tint,nint,tf,tl,ipln)
c
c     Reset statistics arrays
c
      if (ipln.eq.1) then
         do i=1,np
            grafm(i)=0.
         end do
      end if
c
c     Add upp statistics
c
      do i=1,np
        grafm(i)=grafm(i)+c*graf(i)
      end do

      return

      end subroutine accu
