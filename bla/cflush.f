c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/bla/cflush.f $
c $LastChangedDate: 2007-04-21 16:01:24 +0200 (Sat, 21 Apr 2007) $
c $LastChangedBy: pschlatt@MECH.KTH.SE $
c $LastChangedRevision: 532 $
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
