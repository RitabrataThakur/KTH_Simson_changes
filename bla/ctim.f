c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/bla/ctim.f $
c $LastChangedDate: 2007-11-12 13:19:01 +0100 (Mon, 12 Nov 2007) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 850 $
c
c ***********************************************************************
      subroutine ctim(ctime,wtime)
c
c     Get wall time and CPU time
c
      implicit none

      real wtime,ctime

      call cpu_time(ctime)
      call wall_time(wtime)

      end subroutine ctim
