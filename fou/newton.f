c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/fou/newton.f $
c $LastChangedDate: 2007-11-12 13:19:01 +0100 (Mon, 12 Nov 2007) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 850 $
c
c ***********************************************************************
      subroutine newton(z,a,deg)
c
c     Finds a zero of the derivative of a deg degre polynomial
c     close to z. The coefficients of the polynomial are given
c     in a
c
      implicit none

      integer deg,i,j
      real z,a(10),zold,dp,ddp,tmp

      tmp=z
      do i=1,100
         dp=0.
         ddp=00.
         do j=1,deg
            dp=dp+(j-1)*a(j)*z**(j-2)
            ddp=ddp+(j-1)*(j-2)*a(j)*z**(j-3)
         end do
         zold=z
         z=zold-dp/ddp
         write(52,*) z,dp,ddp
         if (abs(z-zold).lt.1e-2) goto 5
      end do

 5    if (i.ge.100) then
         write(*,*) 'No convergence in 100 iterations'
         z=tmp
      end if

      return

      end subroutine newton
