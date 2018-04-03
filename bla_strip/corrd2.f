c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/bla/corrd2.f $
c $LastChangedDate: 2007-11-12 13:19:01 +0100 (Mon, 12 Nov 2007) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 850 $
c
c ***********************************************************************
      subroutine corrd2(a,n,m,md)

      implicit none

      integer n,m,md
      real a(md,n)
      integer i,j

      if (mod(n-1,2).eq.0) then
         do j=1,m
            a(j,1)=a(j,1)-.5*a(j,n)
         end do
         do i=3,n,2
            do j=1,m
               a(j,i-1)=a(j,i-1)-
     &              (1.+real((n-2)**2-(i-2)**2)/real(4*n-4))*a(j,n-1)
               a(j,i)=a(j,i)-a(j,n)
            end do
         end do
      else
         do j=1,m
            a(j,1)=a(j,1)-.5*(1.+real((n-2)**2)/real(4*n-4))*a(j,n-1)
            a(j,2)=a(j,2)-a(j,n)
         end do
         do i=4,n,2
            do j=1,m
               a(j,i-1)=a(j,i-1)-
     &              (1.+real((n-2)**2-(i-2)**2)/real(4*n-4))*a(j,n-1)
               a(j,i)=a(j,i)-a(j,n)
            end do
         end do
      end if

      end subroutine corrd2
