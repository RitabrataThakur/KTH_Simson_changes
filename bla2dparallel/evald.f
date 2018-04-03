c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/evald.f $
c $LastChangedDate: 2007-10-29 14:23:13 +0100 (Mon, 29 Oct 2007) $
c $LastChangedBy: qiang@MECH.KTH.SE $
c $LastChangedRevision: 786 $
c
c ***********************************************************************
      subroutine evald(dj1ur,dj1ui,dj1uh1,dj1uh2,dj1uh3,
     &     dj1omr,dj1omi,dj1omh,
     &     d2vr,d2vi,domyr,domyi,
     &     d2omyr,d2omyi,d2vh,d2vh2,domyh,d2omyh,
     &     bbeta,k2i,alfa,ibc,xb)
c
c     Evaluates the derivatives of the 
c     particular and homogeneous solutions at y=1
c     with the nonlinear terms calculated in step 2 as a driving force
c     
c     Note that the use of ibc is only to avoid unnecessary
c     parts of the code. The code should work if all
c     lines are executed regardless of the value of ibc
c
      implicit none

      include 'par.f'
c
c     Evaluated derivatives at the boundary y=1 
c
      real dj1ur(memnx*mbz,3,5),dj1ui(memnx*mbz,3,5)
      real dj1uh1(memnx*mbz,3,5),dj1uh2(memnx*mbz,3,5)
      real dj1uh3(memnx*mbz,3,4)
      real dj1omr(memnx*mbz,5),dj1omi(memnx*mbz,5)
      real dj1omh(memnx*mbz,5)
c
c     Field variables
c
      real d2vr(memnx*mbz,nyp),d2vi(memnx*mbz,nyp)
      real domyr(memnx*mbz,nyp),domyi(memnx*mbz,nyp)
      real d2omyr(memnx*mbz,nyp),d2omyi(memnx*mbz,nyp)
      real d2vh(memnx*mbz,nyp)
      real d2vh2(memnx*mbz,nyp)
      real domyh(memnx*mbz,nyp),d2omyh(memnx*mbz,nyp)
c
c     Wavenumbers
c
      real bbeta(nx/2*mbz),k2i(nx/2*mbz),alfa(nx/2*mbz)
      
      integer ibc
c
c     Local variables
c
      integer j,xz,nxz,y,xb,x
      real fy,ffy
      
      nxz=memnx*mbz
c
c     We make use of the following properties of the variables at y=1 :
c
c     vr=vi=dvr=dvi=0
c     vh=0,dvh=2
c     vh2=2,dvh2=0
c     omyr=omyi=0, omyh=2
c
c     As a bc :
c     dj1q(:,i,j+1) is the jth derivative of the q
c     where q is one of vr,vi,vh1,vh2,omr,omi or omh
c     calculate the derivatives of each cheb poly
c     at x=1 :
c
      do j=0,4
         do xz=1,nxz
            dj1ur(xz,2,j+1)=0.0
            dj1ui(xz,2,j+1)=0.0
            dj1uh1(xz,2,j+1)=0.0
            dj1uh2(xz,2,j+1)=0.0
            dj1omr(xz,j+1)=0.0
            dj1omi(xz,j+1)=0.0
            dj1omh(xz,j+1)=0.0
         end do
      end do
      do y=1,ny
         do xz=1,nxz
            dj1ur(xz,2,3)=dj1ur(xz,2,3)+d2vr(xz,y)
            dj1ui(xz,2,3)=dj1ui(xz,2,3)+d2vi(xz,y)
            dj1uh1(xz,2,3)=dj1uh1(xz,2,3)+d2vh(xz,y)
            dj1uh2(xz,2,3)=dj1uh2(xz,2,3)+d2vh2(xz,y)
         end do
      end do
      do xz=1,nxz
        dj1omh(xz,1)=2.0
        dj1uh2(xz,2,1)=2.0
        dj1uh1(xz,2,2)=2.0
      end do
      if (ibc.ne.20.and.ibc.ne.120) then
         do y=1,ny
            fy=real(y-1)
            fy=fy*fy
            ffy=fy*(fy-1.)/3.
            do xz=1,nxz
c
c     d(Tk(y))dy (y=1) = k*k
c
               dj1ur(xz,2,4)=dj1ur(xz,2,4)+d2vr(xz,y)*fy
               dj1ui(xz,2,4)=dj1ui(xz,2,4)+d2vi(xz,y)*fy
               dj1uh1(xz,2,4)=dj1uh1(xz,2,4)+d2vh(xz,y)*fy
               dj1uh2(xz,2,4)=dj1uh2(xz,2,4)+d2vh2(xz,y)*fy
c
c     d2(Tk(y))dy2 (y=1) = 1/3*(k*k)*(k*k-1)
c
               dj1ur(xz,2,5)=dj1ur(xz,2,5)+d2vr(xz,y)*ffy
               dj1ui(xz,2,5)=dj1ui(xz,2,5)+d2vi(xz,y)*ffy
               dj1uh1(xz,2,5)=dj1uh1(xz,2,5)+d2vh(xz,y)*ffy
               dj1uh2(xz,2,5)=dj1uh2(xz,2,5)+d2vh2(xz,y)*ffy
               dj1omr(xz,2)=dj1omr(xz,2)+domyr(xz,y)
               dj1omi(xz,2)=dj1omi(xz,2)+domyi(xz,y)
               dj1omh(xz,2)=dj1omh(xz,2)+domyh(xz,y)
               dj1omr(xz,3)=dj1omr(xz,3)+d2omyr(xz,y)
               dj1omi(xz,3)=dj1omi(xz,3)+d2omyi(xz,y)
               dj1omh(xz,3)=dj1omh(xz,3)+d2omyh(xz,y)
               dj1omr(xz,4)=dj1omr(xz,4)+d2omyr(xz,y)*fy
               dj1omi(xz,4)=dj1omi(xz,4)+d2omyi(xz,y)*fy
               dj1omh(xz,4)=dj1omh(xz,4)+d2omyh(xz,y)*fy
            end do
         end do
      end if
c      
c     Now generate the derivatives on u and w
c     u=1/k2(i*alfa*dv-i*beta*omy)
c     w=1/k2(i*beta*dv+i*alfa*omy)
c
      do j=1,4
         do xz=1,nxz
            x=xz+xb
            dj1ur(xz,1,j)=
     &           (-alfa(x)*dj1ui(xz,2,j+1)+
     &           bbeta(x)*dj1omi(xz,j))*k2i(x)
            dj1ui(xz,1,j)=
     &           (alfa(x)*dj1ur(xz,2,j+1)-
     &           bbeta(x)*dj1omr(xz,j))*k2i(x)
            dj1ur(xz,3,j)=
     &           (-bbeta(x)*dj1ui(xz,2,j+1)-
     &           alfa(x)*dj1omi(xz,j))*k2i(x)
            dj1ui(xz,3,j)=
     &           (bbeta(x)*dj1ur(xz,2,j+1)+
     &           alfa(x)*dj1omr(xz,j))*k2i(x)
            dj1uh1(xz,1,j)=alfa(x)*dj1uh1(xz,2,j+1)*k2i(x)
            dj1uh1(xz,3,j)=bbeta(x)*dj1uh1(xz,2,j+1)*k2i(x)
            dj1uh2(xz,1,j)=alfa(x)*dj1uh2(xz,2,j+1)*k2i(x)
            dj1uh2(xz,3,j)=bbeta(x)*dj1uh2(xz,2,j+1)*k2i(x)
            dj1uh3(xz,1,j)=-bbeta(x)*dj1omh(xz,j)*k2i(x)
            dj1uh3(xz,2,j)=0.0
            dj1uh3(xz,3,j)=alfa(x)*dj1omh(xz,j)*k2i(x)
         end do
      end do

      end subroutine evald
