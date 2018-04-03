c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/bla/xlim.f $
c $LastChangedDate: 2006-12-02 01:17:16 +0100 (Sat, 02 Dec 2006) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 372 $
c
c ***********************************************************************
      real function xlim(x,xlow)
c
c     Returns the argument if argument is larger than xlow
c     If less than this it returns a value which
c     is always at least xlow/2
c     xlim has two continuous derivatives
c
      implicit none

      real x,xlow
      if (x.ge.xlow) then
         xlim=x
      else
         if (x.le.0.) then
            xlim=xlow*0.5
         else
            xlim=xlow*(0.5+(x/xlow)**3-0.5*(x/xlow)**4)
         end if
      end if

      return

      end function xlim
