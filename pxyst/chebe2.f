c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/pxyst/chebe2.f $
c $LastChangedDate: 2008-07-03 13:15:41 +0200 (Thu, 03 Jul 2008) $
c $LastChangedBy: pschlatt@MECH.KTH.SE $
c $LastChangedRevision: 1248 $
c
c ***********************************************************************
      subroutine chebe2(e,a,b,c,n,m,x)
c
c     Multiple chebyshev interpolation evaluation
c     performs n simultaneous interpolations
c
      integer n,m
      real a,b,x
      real e(n),c(n,m)
      real d(5049),dd(5049),sv(5049),y,y2

      integer i,j

      if ((x-a)*(x-b).gt.0.) stop
      if (n.gt.5049) then
         write(*,*) 'increase n in chebe2.f'
         stop
      end if
      do i=1,n
         d(i)=0.
         dd(i)=0.
      end do
      y=(2.*x-a-b)/(b-a)
      y2=2.*y
      do j=m,2,-1
         do i=1,n
            sv(i)=d(i)
            d(i)=y2*d(i)-dd(i)+c(i,j)
            dd(i)=sv(i)
         end do
      end do
      do i=1,n
         e(i)=y*d(i)-dd(i)+c(i,1)
      end do

      end subroutine chebe2
