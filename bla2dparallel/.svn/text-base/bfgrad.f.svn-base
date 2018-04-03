      subroutine bfgrad(bu1,bu2,bfgrad1,bfgrad2,alfa,prex,prey)
c
c     Computes vorticity associated with base flow
c
      implicit none

      include 'par.f'

      real bu1(nxp/2+1,nyp,3),bu2(nxp/2+1,nyp,3)
      real wbr(nxp/2+1,nyp),wbi(nxp/2+1,nyp)
      real bfgrad1(nxp/2+1,nyp,3,2),bfgrad2(nxp/2+1,nyp,3,2)
      real app1(nyp),app2(nyp)
      real appr(nxp/2+1),appi(nxp/2+1)
      real app2r(nxp/2+1),app2i(nxp/2+1)
      real prey(nyp*2+15),w(nxp/2+1,nyp)
      real alfa(nx/2*mbz),prex(nxp+15)

      integer x,i,y
c
c     y derivatives
c
      do i=1,3
         do x=1,nxp/2+1
            do y=1,nyp
               app1(y)=bu1(x,y,i)
               app2(y)=bu2(x,y,i)
            end do
            call vchbf(app1,w,nyp,1,1,1,prey)
            call vchbf(app2,w,nyp,1,1,1,prey)
            call rdcheb(app1,nyp,1,1)
            call rdcheb(app2,nyp,1,1)
            call vchbb(app1,w,nyp,1,1,1,prey)
            call vchbb(app2,w,nyp,1,1,1,prey)
            do y=1,nyp
               bfgrad1(x,y,i,2)=app1(y)*(2./real(nyp-1))
               bfgrad2(x,y,i,2)=app2(y)*(2./real(nyp-1))
            end do
         end do
      end do
c
c     x-derivatives 
c
      do i=1,3
         do y=1,nyp
            do x=1,nxp/2
               appr(x)=bu1(x,y,i)
               appi(x)=bu2(x,y,i)
            end do
            call vrfftf(appr,appi,wbr,wbi,nxp,1,1,1,prex)
            do x=1,nxp/2
               app2r(x)=-alfa(x)*appi(x)*(1./real(nxp))
               app2i(x)=alfa(x)*appr(x)*(1./real(nxp))
            end do
            call vrfftb(app2r,app2i,wbr,wbi,nxp,1,1,1,prex)
            do x=1,nxp/2
               bfgrad1(x,y,i,1)=app2r(x)
               bfgrad2(x,y,i,1)=app2i(x)
            end do
         end do
      end do

      do i=1,3
         do y=1,nyp
            bfgrad1(nxp/2+1,y,i,1)=bfgrad1(1,y,i,1)
            bfgrad2(nxp/2+1,y,i,1)=bfgrad2(1,y,i,1)
            bfgrad1(nxp/2+1,y,i,2)=bfgrad1(1,y,i,2)
            bfgrad2(nxp/2+1,y,i,2)=bfgrad2(1,y,i,2)
            write(82,*)y,bfgrad1(1,y,i,1),bfgrad1(1,y,i,2)
         end do
      end do

      end subroutine bfgrad
