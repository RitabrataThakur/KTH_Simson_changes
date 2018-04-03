c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/cflush.f $
c $LastChangedDate: 2007-10-29 14:23:13 +0100 (Mon, 29 Oct 2007) $
c $LastChangedBy: qiang@MECH.KTH.SE $
c $LastChangedRevision: 786 $
c
c ***********************************************************************
      subroutine cflush(iunit)
c
c     Compiler-dependent flush of unit iunit
c
      implicit none
      integer iunit

#ifdef IBM
c
c     IBM xlf version
c
      call flush_(iunit)
      return
#endif
c
c     Sun hp sgi dec linux version
c
      call flush(iunit)

      end subroutine cflush
