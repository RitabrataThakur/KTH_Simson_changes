c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/ctim.f $
c $LastChangedDate: 2007-10-29 14:23:13 +0100 (Mon, 29 Oct 2007) $
c $LastChangedBy: qiang@MECH.KTH.SE $
c $LastChangedRevision: 786 $
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
