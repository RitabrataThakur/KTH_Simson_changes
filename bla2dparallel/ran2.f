c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/ran2.f $
c $LastChangedDate: 2007-10-29 14:23:13 +0100 (Mon, 29 Oct 2007) $
c $LastChangedBy: qiang@MECH.KTH.SE $
c $LastChangedRevision: 786 $
c
c ***********************************************************************
      real function ran2(idum)
c
c     A simple portable random number generator
c
c     Requires 32-bit integer arithmetic
c     Taken from Numerical Recipes, William Press et al.
c     gives correlation free random numbers but does not have a very large
c     dynamic range, i.e only generates 714025 different numbers
c     for other use consult the above
c     Set idum negative for initialization
c
      implicit none

      integer idum,ir(97),m,ia,ic,iff,iy,j
      real rm
      parameter (m=714025,ia=1366,ic=150889,rm=1./m)
      save iff,ir,iy
      data iff /0/
      
      if (idum.lt.0.or.iff.eq.0) then
c
c     Initialize
c
         iff=1
         idum=mod(ic-idum,m)
         do j=1,97
            idum=mod(ia*idum+ic,m)
            ir(j)=idum
         end do
         idum=mod(ia*idum+ic,m)
         iy=idum        
      end if
c
c     Generate random number
c
      j=1+(97*iy)/m
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      
      end function ran2
