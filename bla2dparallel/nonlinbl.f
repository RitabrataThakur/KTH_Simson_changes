c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/nonlinbl.f $
c $LastChangedDate: 2011-05-19 11:13:29 +0200 (Thu, 19 May 2011) $
c $LastChangedBy: pschlatt@MECH.KTH.SE $
c $LastChangedRevision: 1666 $
c
c ***********************************************************************
      subroutine nonlinbl(amp,campw,kx,kz,nwave,cflp,vext,cext,iext,
     &     ur,ui,tc,xsc,zsc,xs,yb,rot,spat,fltype,
     &     fstart,fend,nst,
     &     bu1,bu2,
     &     loctyp,fp1,fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7,
     &     fpds8,fpdds1,fpdds2,fpdds3,fpdds4,fpdds5,
     &     tripf,txsc,tysc,ttzc,tx0,
     &     osmod,osnumb,
     &     osur,osui,osvr,osvi,oswr,oswi,
     &     xc1,xc2,fring1,fring2,cphas1,cphas2,sphas1,sphas2,
     &     evr,evi,eur,eui,ewr,ewi,afw,alfaw,betaw,ampob,amp2d,
     &     ev2r,ev2i,eu2r,eu2i,afw2d,alf2d,
     &     fbla,nbla,dybla,rlam,cdev,re,xblc,
     &     alfa,beta,eta,deta,xl,zl,prex,prez,pres,prea,
     &     it,icfl,iamp,u2r,u2i,om2r,om2i,th2r,th2i,wr,wi,
     &     iles,gur,gui,chi_scaled,prt,
     &     my_node_world,my_node_z,my_node_x,my_comm_z,my_comm_x,
     &     realg1,realg2,realg3,realg4,realg5,realg6,
     &     suction,vsuc,ampst,streak,betast,omegast,
     &     ndxst,uust_r,uust_i,vvst_r,vvst_i,wwst_r,wwst_i,
     &     tsmoo,tsmst,tsmend,iampst,phist,
     &     isfd,sfd_chi,sfd_delta,sfdur,sfdui,an,bn,wbci,
     &     pert,lin,bom1,bom2,x0,spanv,rv1r,rv1i,do_press_nl,
     &     tcom,tser,gsr,gsi,fmhdr,fmhdi,imhd,du2r,du2i,
     &     b0,mhd_n,
     &     boxr_z,boxi_z,boxr,boxi,wr_z,wi_z,wr_x,wi_x,
     &     bf3,u2bf3_r,u2bf3_i,ybp,
     &     adjoint,bfgrad1,bfgrad2,u2xr,u2xi,u2zr,u2zi,
     &     boxrz_z,boxiz_z,boxrx_z,boxix_z,boxrs_z,boxis_z,
     &     tempr,tempi,aaar,aaai)
c      
c     Calculates the nonlinear terms, accumulates amp and cfl for a xz-box
c     i.e. nby times xz-planes
c
      implicit none

      include 'par.f'
#ifdef MPI
      include 'mpif.h'
#endif
      logical spat
      integer nwave,yb,it,icfl,iamp,kx(nwave),kz(nwave),iext,ybp
      complex campw(nyp,4,nwave)
      real cflp,cflp_bf,amp(nyp,20)
      real vext(nyp,2,6),cext(nyp,2,6,2)
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real rot,alfa(nx/2*mby),beta(nz)
      integer fltype,loctyp,nst
      real fp1,fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7
      real fpds8,fpdds1,fpdds2,fpdds3,fpdds4,fpdds5
      logical tripf
      real txsc,tysc,ttzc(nzp+2,4),tx0
      real fstart,fend
      real evr(nyp),evi(nyp),eur(nyp),eui(nyp),ewr(nyp),ewi(nyp)
      real afw,alfaw,betaw,ampob,amp2d
      real ev2r(nyp),ev2i(nyp),eu2r(nyp),eu2i(nyp),afw2d,alf2d
      real eta(nyp),deta(nyp),xl,zl,tc,xsc,zsc,xs,cdev,re,xblc
      real rlam
      real fbla(mbla,7+3*scalar)
      real bu1(nxp/2+1,nyp,3+scalar),bu2(nxp/2+1,nyp,3+scalar)
      real prex(nxp+15)
      real prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)
      real u2r((nxp/2+1)*mby,nzd_new/nprocz,3)
      real u2i((nxp/2+1)*mby,nzd_new/nprocz,3)
      real om2r((nxp/2+1)*mby,nzd_new/nprocz,3)
      real om2i((nxp/2+1)*mby,nzd_new/nprocz,3)
      real du2r((nxp/2+1)*mby,nzd,3,3)
      real du2i((nxp/2+1)*mby,nzd,3,3)
      real wr((nxp/2+1)*mby,nzd),wi((nxp/2+1)*mby,nzd)
      real x0,spanv
      real d1f,d2f,d3f,d2g,dybla
      integer nbla
c
c     Temporal forcing
c
      real c1,xx,etabl,ybl,h3u
      real tempr((nxp/2+1)*mby,nzd_new/nprocz,2)
      real tempi((nxp/2+1)*mby,nzd_new/nprocz,2)
c
c     Linearized-perturbation
c
      logical pert,lin
      real bom1(nxp/2+1,nyp,3),bom2(nxp/2+1,nyp,3)
c
c     Adjoint
c
      logical adjoint
      real bfgrad1(nxp/2+1,nyp,3,2),bfgrad2(nxp/2+1,nyp,3,2)
      real u2xr((nxp/2+1)*mby,nzd_new/nprocz,3)
      real u2xi((nxp/2+1)*mby,nzd_new/nprocz,3)
      real u2zr((nxp/2+1)*mby,nzd_new/nprocz,3)
      real u2zi((nxp/2+1)*mby,nzd_new/nprocz,3)
c
c     Scalar
c
      real th2r((nxp/2+1)*mby,nzd_new/nprocz,4*scalar)
      real th2i((nxp/2+1)*mby,nzd_new/nprocz,4*scalar)
      integer ith
c
c     Pressure
c
      real rv1r(nx/2,nzc,2*pressure)
      real rv1i(nx/2,nzc,2*pressure)
      logical do_press_nl
c
c
c     OS modes
c
      logical osmod
      integer osnumb
      real osur(osnf,nyp), osvr(osnf,nyp), oswr(osnf,nyp)
      real osui(osnf,nyp), osvi(osnf,nyp), oswi(osnf,nyp)
      real cphas1(osnf,nxp/2,nzd_new/nprocz)
      real cphas2(osnf,nxp/2,nzd_new/nprocz)
      real sphas1(osnf,nxp/2,nzd_new/nprocz)
      real sphas2(osnf,nxp/2,nzd_new/nprocz)
 
      real fring1(nxp/2),fring2(nxp/2)
      real xc1(nxp/2),xc2(nxp/2)

      logical sym
      integer i,j,xy,z,zp,nxy,y,y1,x,nxp2,npl,npp,xb,xp
      real h1u,h2u,h1e,h2e
c
c     Suction boundary layer and jet
c
      logical suction
      real vsuc
      integer wbci
c
c     Streak generation
c
      logical streak
      real iampst,ampst(2),tsmoo(4),tsmst(2),tsmend(2)
      real betast(2),omegast(2),phist
      real uust_r(nyp,180,2),vvst_r(nyp,180,2),wwst_r(nyp,180,2)
      real uust_i(nyp,180,2),vvst_i(nyp,180,2),wwst_i(nyp,180,2)
      integer ndxst
c
c     LES
c
      real gur(memnx,memny,memnz,5)
      real gui(memnx,memny,memnz,5)
      real gsr(memnx,memny,memnz,scalar)
      real gsi(memnx,memny,memnz,scalar)
      integer iles,ll
      real chi_scaled,prt
c
c     MHD
c
      real fmhdr(memnx,memny,memnz,2)
      real fmhdi(memnx,memny,memnz,2)
      integer imhd
      real cr1,cr2,cr3
      real ci1,ci2,ci3
      real b0(3),mhd_n,fact
      real ri(scalar)
c
c     SFD
c
      integer isfd
      real sfd_chi,sfd_delta
      real sfdur(memnx,memny,memnz,6)
      real sfdui(memnx,memny,memnz,6)
      real an,bn

c
c     3D baseflow
c
      logical bf3
      real u2bf3_r((nxp/2+1)*mby,nzd_new/nprocz,nyp/nprocz+1,6)
      real u2bf3_i((nxp/2+1)*mby,nzd_new/nprocz,nyp/nprocz+1,6)

c
c     MPI
c
      integer my_node_world
      integer realg1,realg2,realg3,realg4,realg5,realg6
      integer my_node_z,my_node_x,my_comm_z,my_comm_x
#ifdef MPI
      integer ierror
      real boxr_z(memnx,nzd_new,6)
      real boxi_z(memnx,nzd_new,6)
      real boxr(memnx,nzd_new,4*scalar)
      real boxi(memnx,nzd_new,4*scalar)
      real wr_z(memnx,nzd_new)
      real wi_z(memnx,nzd_new)
      real wr_x(nxp/2+1,nzd_new/nprocz)
      real wi_x(nxp/2+1,nzd_new/nprocz)
      real aaar(nxp/2+1,nzd_new,2)
      real aaai(nxp/2+1,nzd_new,2)
c
c     Work arrays needed for adjoint
c
      real boxrz_z(memnx,nzd_new,6)
      real boxiz_z(memnx,nzd_new,6)
      real boxrx_z(memnx,nzd_new,6)
      real boxix_z(memnx,nzd_new,6)
      real boxrs_z(memnx,nzd_new,6)
      real boxis_z(memnx,nzd_new,6)
#endif
c
c     Functions
c     
      real, external :: cubip


      real t1,t2
      real tcom,tser
c
c     If running channels without pressure boundary conditions
c
c      do_press_nl = .false.


#ifdef MPI
c      call mpi_barrier(mpi_comm_world,ierror)
#endif
      call wall_time(t1)


      if (mby.gt.1.and.pert) then
         write(*,*)'mby>1 and pert----> not implemented'
         call stopnow(1009)
      end if

      npp=min(nyp,yb+mby-1)
      npl=min(mby,nyp-yb+1)
      nxp2=nxp/2+1
      nxy=nxp2*npl
      xb=mod(my_node_world,nprocx)*memnx
#ifdef MPI
      if (nzd_new.gt.nzp) then
         do i=1,3
            do z=nzp+1,nzd_new
               do x=1,memnx
                  boxr_z(x,z,i)=0.
                  boxi_z(x,z,i)=0.
                  if (adjoint) then
                     boxrx_z(x,z,i)=0.
                     boxix_z(x,z,i)=0.
                     boxrz_z(x,z,i)=0.
                     boxiz_z(x,z,i)=0.
                     boxrs_z(x,z,i)=0.
                     boxis_z(x,z,i)=0.
                  end if
               end do
            end do
         end do
      end if

c      write(*,*),shape(boxrx_z),shape(boxix_z),shape(boxrz_z),
c     &     shape(boxiz_z),shape(boxrs_z),shape(boxis_z)
#endif
c
c     u2r,u2i,om2r,om2i are on the dealiasing grid
c     middle is padded, i.e. nonzero values:
c     x = 1..nx/2
c          nx/2
c     z = 1..nz/2, nzp-nz/2+2..nzp
c           nz/2      nz/2-1
c     (if ipad=0 z=1..nz/2, nz-nz/2+2)
c                  nz/2       nz/2-1

c
c     Get the velocities (padded)
c
      if (nproc.eq.1) then
c
c     All code for nproc = 1 removed, it can be found in older versions!!!
c
      else
#ifdef MPI
c
c     Get velocities into boxr_z(:,:,1-3)
c
         do i=1,3
            call getpxz_z(boxr_z(1,1,i),boxi_z(1,1,i),yb,i,1,ur,ui,
     &           realg1,realg2,my_node_z,my_comm_z,my_node_world)
         end do
c
c     Save velocities if computing the adjoint
c
         if (adjoint) then
            do x=1,memnx
               do z=1,nzd_new
                  do i=1,6
                     boxrs_z(x,z,i) = boxr_z(x,z,i)
                     boxis_z(x,z,i) = boxi_z(x,z,i)
                  end do
               end do
            end do

         end if
c
c     Get scalars into boxr(:,:,1+...)
c     and derivatives into boxr(:,:,3+...)
c
         do ith=1,scalar
            call getpxz_z(boxr(1,1,1+4*(ith-1)),boxi(1,1,1+4*(ith-1)),yb
     &           ,8+pressure+3*(ith-1),1,ur,ui,
     &           realg1,realg2,my_node_z,my_comm_z,my_node_world)
            call getpxz_z(boxr(1,1,3+4*(ith-1)),boxi(1,1,3+4*(ith-1)),yb
     &           ,9+pressure+3*(ith-1),1,ur,ui,
     &           realg1,realg2,my_node_z,my_comm_z,my_node_world)
         end do

c         call mpi_barrier(mpi_comm_world,ierror)
         call wall_time(t2)
         tcom = tcom + (t2-t1)
         t1=t2             
c
c     Compute the wall-normal vorticity
c
ccccc symmetry is not considered here
         do z=1,nz
            zp=nzp-nz+z
            do y=yb,npp
               y1=(y-yb)*nxp2
               do xp=1,memnx
                  x=xp+xb
                  xy=xp+y1
                  if (z.le.nz/2) then
                     boxr_z(xy,z,5)=-beta(z)*boxi_z(xy,z,1)
     &                    +alfa(x)*boxi_z(xy,z,3)
                     boxi_z(xy,z,5)= beta(z)*boxr_z(xy,z,1)
     &                    -alfa(x)*boxr_z(xy,z,3)
                  else
                     boxr_z(xy,zp,5)=-beta(z)*boxi_z(xy,zp,1)
     &                    +alfa(x)*boxi_z(xy,zp,3)
                     boxi_z(xy,zp,5)= beta(z)*boxr_z(xy,zp,1)
     &                    -alfa(x)*boxr_z(xy,zp,3) 
                  end if
               end do
            end do
         end do
      
c
c     Compute the derivatives of theta th2(:,:,1) in 2,3,4
c
ccccc symmetry is not considered here
         do ith=1,scalar
            do z=1,nz
               zp=nzp-nz+z
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do xp=1,memnx
                     x=xp+xb
                     xy=xp+y1
                     if (z.le.nz/2) then
                        boxr(xy,z,4+4*(ith-1))=-beta(z)
     &                       *boxi(xy,z,1+4*(ith-1))
                        boxi(xy,z,4+4*(ith-1))= beta(z)
     &                       *boxr(xy,z,1+4*(ith-1))

                        boxr(xy,z,2+4*(ith-1))=-alfa(x)
     &                       *boxi(xy,z,1+4*(ith-1))
                        boxi(xy,z,2+4*(ith-1))= alfa(x)
     &                       *boxr(xy,z,1+4*(ith-1))
                     else
                        boxr(xy,zp,4+4*(ith-1))=-beta(z)
     &                       *boxi(xy,zp,1+4*(ith-1))
                        boxi(xy,zp,4+4*(ith-1))= beta(z)
     &                       *boxr(xy,zp,1+4*(ith-1))
                        
                        boxr(xy,zp,2+4*(ith-1))=-alfa(x)
     &                       *boxi(xy,zp,1+4*(ith-1))
                        boxi(xy,zp,2+4*(ith-1))= alfa(x)
     &                       *boxr(xy,zp,1+4*(ith-1))
                     end if
                  end do
               end do
            end do
         end do
c     
c     And pad/oddball
c     
         if (nfxd.eq.1.and.nfyd.eq.0.and.nfzd.eq.1) then
            do z=(nz+1)/2+1,nzp+1-nz/2
               do x=1,memnx
                  boxr_z(x,z,5)=0.0
                  boxi_z(x,z,5)=0.0
               end do
            end do
         else
            if (mod(nz,2).eq.0) then
               do x=1,memnx
                  boxr_z(x,nz/2+1,5)=0.0
                  boxi_z(x,nz/2+1,5)=0.0
               end do
            end if
         end if
         do ith=1,scalar
            if (nfxd.eq.1.and.nfyd.eq.0.and.nfzd.eq.1) then
               do z=(nz+1)/2+1,nzp+1-(nz+1)/2
                  do xy=1,memnx*npl
                     boxr(xy,z,2+4*(ith-1))=0.0
                     boxi(xy,z,2+4*(ith-1))=0.0
                     
                     boxr(xy,z,4+4*(ith-1))=0.0
                     boxi(xy,z,4+4*(ith-1))=0.0
                  end do
               end do
            else
               do z=1,nzd_new
                  do y=yb,npp
                     y1=(y-yb)*nxp2
                     do x=nx/2+1,nxp2
                        xy=x+y1
                        boxr(xy,z,2+4*(ith-1))=0.0
                        boxi(xy,z,2+4*(ith-1))=0.0
                     
                        boxr(xy,z,4+4*(ith-1))=0.0
                        boxi(xy,z,4+4*(ith-1))=0.0
                     end do
                  end do
               end do
            end if
         end do
         
c         call mpi_barrier(mpi_comm_world,ierror)
         call wall_time(t2)
         tser = tser + (t2-t1)
         t1=t2
c
c     Get vorticities into boxr_z(:,:,4) and boxr_z(:,:,6)

         call getpxz_z(boxr_z(1,1,4),boxi_z(1,1,4),yb,4,1,ur,ui,
     &        realg1,realg2,my_node_z,my_comm_z,my_node_world)
         call getpxz_z(boxr_z(1,1,6),boxi_z(1,1,6),yb,5,1,ur,ui,
     &        realg1,realg2,my_node_z,my_comm_z,my_node_world)

c         call mpi_barrier(mpi_comm_world,ierror)
         call wall_time(t2)
         tcom = tcom + (t2-t1)
         t1=t2

         if (pressure.eq.1.and.do_press_nl) then         
            if (yb.eq.1) then
               call getpxz_x(tempr(1,1,1),tempi(1,1,1),yb,2,1,
     &              boxr_z(1,1,2),boxi_z(1,1,2),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
            end if
            if (yb.eq.nyp) then
               call getpxz_x(tempr(1,1,2),tempi(1,1,2),yb,2,1,
     &              boxr_z(1,1,2),boxi_z(1,1,2),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
            end if
         end if
c     
c     Backward Fourier transform in z of velocities, vorticities and scalars/derivatives
c     
         do i=1,6    
            call vcfftb(boxr_z(1,1,i),boxi_z(1,1,i),wr_z,wi_z,
     &           nzp,memnx,memnx,1,prez)
         end do
         do ith=1,scalar
            do i=1,4
               call vcfftb(boxr(1,1,i+4*(ith-1)),boxi(1,1,i+4*(ith-1)),
     &              wr_z,wi_z,nzp,memnx,memnx,1,prez)
            end do
         end do

c         call mpi_barrier(mpi_comm_world,ierror)
         call wall_time(t2)
         tser = tser + (t2-t1)
         t1=t2      
c
c     Get velocities/vorticities into u2r(:,:,1-3) and om2r(:,:,1-3)

         do i=1,3
            call getpxz_x(u2r(1,1,i),u2i(1,1,i),yb,i,1,
     &           boxr_z(1,1,i),boxi_z(1,1,i),
     &           realg3,realg4,my_node_x,my_comm_x,my_node_world)
            call getpxz_x(om2r(1,1,i),om2i(1,1,i),yb,i+3,1,
     &           boxr_z(1,1,i+3),boxi_z(1,1,i+3),
     &           realg3,realg4,my_node_x,my_comm_x,my_node_world)
         end do
c
c     Get scalars/derivatives into th2r(:,:,1-4)
c         
         do ith=1,scalar
            do i=1,4
               call getpxz_x(th2r(1,1,i+4*(ith-1)),th2i(1,1,i+4*(ith-1))
     &              ,yb,8+pressure+3*(ith-1),1,
     &              boxr(1,1,i+4*(ith-1)),boxi(1,1,i+4*(ith-1)),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
            end do
         end do

c         call mpi_barrier(mpi_comm_world,ierror)
         call wall_time(t2)
         tcom = tcom + (t2-t1)
         t1=t2
c     
c     Backward Fourier transform in x
c     
         do i=1,3
            call vrfftb(u2r(1,1,i),u2i(1,1,i),wr_x,wi_x,
     &           nxp,nzd_new/nprocz,1,nxp/2+1,prex)
            call vrfftb(om2r(1,1,i),om2i(1,1,i),wr_x,wi_x,
     &           nxp,nzd_new/nprocz,1,nxp/2+1,prex)
         end do
          do ith=1,scalar
            do i=1,4
               call vrfftb(th2r(1,1,i+4*(ith-1)),th2i(1,1,i+4*(ith-1)),
     &              wr_x,wi_x,nxp,nzd_new/nprocz,1,nxp/2+1,prex)
            end do
         end do

c         call mpi_barrier(mpi_comm_world,ierror)
         call wall_time(t2)
         tser = tser + (t2-t1)
         t1=t2

#endif
      end if

c
c     Save the velocities for MHD
c     for later use when computing the Lorentz force
c
      if (imhd.eq.1) then
               write(*,*)'Not Implemented'
               call stopnow(33445566)
         do ll=1,3
            do z=1,nz/2
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     du2r(xy,z,ll,2)=u2r(xy,z,ll)
                     du2i(xy,z,ll,2)=u2i(xy,z,ll)
                  end do
               end do
            end do
            do z=nz/2+1,nz
               zp=nzp-nz+z
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     du2r(xy,zp,ll,2)=u2r(xy,zp,ll)
                     du2i(xy,zp,ll,2)=u2i(xy,zp,ll)
                  end do
               end do
            end do
         end do
      end if
c
c     Computations necessary for adjoint code
c
      if (adjoint) then
c
c     Derivatives of velocity in x and z directions, using the saved
c     velocities in boxis_r and boxrs_r
c
         do z=1,nz
            zp=nzp-nz+z
            do y=yb,npp
               y1=(y-yb)*nxp2
               do xp=1,memnx
                  x=xp+xb
                  xy=xp+y1
                  do i=1,3
                     if (z.le.nz/2) then
                        boxrz_z(xy,z,i)=-beta(z)*boxis_z(xy,z,i)
                        boxiz_z(xy,z,i)=beta(z)*boxrs_z(xy,z,i)
                        boxrx_z(xy,z,i)=-alfa(x)*boxis_z(xy,z,i)
                        boxix_z(xy,z,i)=alfa(x)*boxrs_z(xy,z,i)
                     else
                        boxrz_z(xy,zp,i)=-beta(z)*boxis_z(xy,zp,i)
                        boxiz_z(xy,zp,i)=beta(z)*boxrs_z(xy,zp,i)
                        boxrx_z(xy,zp,i)=-alfa(x)*boxis_z(xy,zp,i)
                        boxix_z(xy,zp,i)=alfa(x)*boxrs_z(xy,zp,i)
                     end if
                  end do
               end do
            end do
         end do
c
c     Add pad/oddball
c
         if (nfxd.eq.1.and.nfyd.eq.0.and.nfzd.eq.1) then
            do z=(nz+1)/2+1,nzp+1-nz/2
               do x=1,memnx
                  do i=1,6
                     boxrx_z(x,z,i)=0.0
                     boxix_z(x,z,i)=0.0
                     boxrz_z(x,z,i)=0.0
                     boxiz_z(x,z,i)=0.0
                  end do
               end do
            end do
         else
            if (mod(nz,2).eq.0) then
               do x=1,memnx
                  do i=1,6
                     boxrx_z(x,nz/2+1,i)=0.0
                     boxix_z(x,nz/2+1,i)=0.0
                     boxrz_z(x,nz/2+1,i)=0.0
                     boxiz_z(x,nz/2+1,i)=0.0
                  end do
               end do
            end if
         end if
c
c     Now do the Fourier synthesis (backward transform) in z first
c
         do i=1,3    
            call vcfftb(boxrx_z(1,1,i),boxix_z(1,1,i),wr_z,wi_z,
     &           nzp,memnx,memnx,1,prez)
            call vcfftb(boxrz_z(1,1,i),boxiz_z(1,1,i),wr_z,wi_z,
     &           nzp,memnx,memnx,1,prez)
         end do
c
c     Get gradients into u2xr(:,:,1-3) and u2xi(:,:,1-3)
c
         do i=1,3
            call getpxz_x(u2xr(1,1,i),u2xi(1,1,i),yb,i,1,
     &           boxrx_z(1,1,i),boxix_z(1,1,i),
     &           realg3,realg4,my_node_x,my_comm_x,my_node_world)
            call getpxz_x(u2zr(1,1,i),u2zi(1,1,i),yb,i,1,
     &           boxrz_z(1,1,i),boxiz_z(1,1,i),
     &           realg3,realg4,my_node_x,my_comm_x,my_node_world)
         end do
c
c     Backward Fourier transform in x
c
         do i=1,3
            call vrfftb(u2xr(1,1,i),u2xi(1,1,i),wr_x,wi_x,
     &           nxp,nzd_new/nprocz,1,nxp/2+1,prex)
            call vrfftb(u2zr(1,1,i),u2zi(1,1,i),wr_x,wi_x,
     &           nxp,nzd_new/nprocz,1,nxp/2+1,prex)
         end do
c
c     Now we have the gradients needed for the adjoint code on the fine
c     grid padded/without oddball
c
      end if

c
c     Now we have all velocities and vorticities on the
c     fine grid padded/without oddball
c
      
c
c     Compute/accumulate the amplitudes (PHYSICAL SPACE!!!)
c
      if (mod(it-1-nst,iamp).eq.0.and.iamp.gt.0) then
         call boxamp(amp,campw,kx,kz,nwave,
     &        u2r,u2i,om2r,om2i,th2r,th2i,
     &        yb,alfa,beta,my_node_x,my_node_z,ybp,
     &        my_node_world)
      end if
c     
c     Pressure: store v at the wall and in the free stream for 
c               Neumann boundary condition at the last full time step
c
c     ATTENTION: MPI 
c      if (pressure.eq.1.and.mod(it-1,nst).eq.0) then
c      do_press_nl=.true.

      if (pressure.eq.1.and.do_press_nl) then
         if (yb.eq.1.or.yb.eq.(nby-1)*mby+1) then
c
c     i is 1 for yb=1 and i is 2 for yb=(nby-1)*mby+1 (usually nyp)
c     
            i = min(2,yb)
            y1 = (nyp-yb)*nxp2*(1-1/yb)
            do z = 1,1/nz+nz/2*(2-nfzsym)
               zp = z+z/(nz/2+1)*(nzp-nz)
               do x=1,nx/2
                  rv1r(x,z,i) = u2r(x+y1,zp,2)
                  rv1i(x,z,i) = u2i(x+y1,zp,2)
c
c     For Dirichlet boundary conditions one would use
c     rv1r(x,z,i)=-alfa(x)*u2i(x+y1,zp,1)-beta(z)*u2i(x+y1,zp,3)
c     rv1i(x,z,i)= alfa(x)*u2r(x+y1,zp,1)+beta(z)*u2r(x+y1,zp,3)
c
               end do
            end do
         end if

#ifdef MPI
         if (nproc.gt.1) then
            if (yb.eq.1.or.yb.eq.(nby-1)*mby+1) then
               i = min(2,yb)
               zp=mod(my_node_world,nprocx)*nzd_new/nprocz
               do z=1,nzd_new/nprocz
                  do x=1,nxp/2+1
                     aaar(x,z+zp,i)=tempr(x,z,i)
                     aaai(x,z+zp,i)=tempi(x,z,i)
                  end do
               end do
            end if
c
c     Communicate planes if necessary
c
c     yb=1
c     
            if (yb.le.nprocz) then
               do i=0,nprocx-1
                  call mpi_bcast(aaar(1,1+i*nzd_new/nprocz,1),
     &                 (nxp/2+1)*nzd_new/nprocz,
     &                 mpi_double_precision,i,mpi_comm_world,ierror)
                  call mpi_bcast(aaai(1,1+i*nzd_new/nprocz,1),
     &                 (nxp/2+1)*nzd_new/nprocz,
     &                 mpi_double_precision,i,mpi_comm_world,ierror)
               end do
            end if
c
c     yb=nyp
c
            if (yb.gt.(nyp/nprocz)*nprocz) then
               do i=mod(nyp-1,nprocz)*nprocx,mod(nyp-1,nprocz)*nprocx
     &              +nprocx-1
                  j=mod(i,nprocx)
                  call mpi_bcast(aaar(1,1+j*nzd_new/nprocz,2),
     &                 (nxp/2+1)*nzd_new/nprocz,
     &                 mpi_double_precision,i,mpi_comm_world,ierror)
                  call mpi_bcast(aaai(1,1+j*nzd_new/nprocz,2),
     &                 (nxp/2+1)*nzd_new/nprocz,
     &                 mpi_double_precision,i,mpi_comm_world,ierror)
               end do
            end if

            if (yb.eq.1.or.yb.eq.(nby-1)*mby+1) then
               i = min(2,yb)
               y1 = (nyp-yb)*nxp2*(1-1/yb)
               do z = 1,1/nz+nz/2*(2-nfzsym)
                  zp = z+z/(nz/2+1)*(nzp-nz)
                  do x=1,nx/2
                     rv1r(x,z,i) = aaar(x+y1,zp,i)
                     rv1i(x,z,i) = aaai(x+y1,zp,i)
                  end do
               end do
            end if
c
c     Communicate planes if necessary
c
c     yb=1
c     
            if (yb.le.nprocz) then
                i=0
               call mpi_bcast(rv1r(1,1,1),nx/2*nzc,mpi_double_precision,
     &              i,mpi_comm_world,ierror)
               call mpi_bcast(rv1i(1,1,1),nx/2*nzc,mpi_double_precision,
     &              i,mpi_comm_world,ierror)
               
            end if
c
c     yb=nyp
c
            if (yb.gt.(nyp/nprocz)*nprocz) then
                i=mod(nyp-1,nprocz)*nprocx
               call mpi_bcast(rv1r(1,1,2),nx/2*nzc,mpi_double_precision,
     &              i,mpi_comm_world,ierror)
               call mpi_bcast(rv1i(1,1,2),nx/2*nzc,mpi_double_precision,
     &              i,mpi_comm_world,ierror)
            end if
         end if
#endif
      end if  


      if (nproc.eq.1) then
c
c     All code for nproc = 1 removed, it can be found in older versions!!!
c
      end if    
c
c     We are now in the physical space (on dealiasing grid)
c     x=1..nxp
c     z=1..nzp
c     y = yb..yb+mby-1 (y=1 is the upper boundary)
c

c
c     Compute CFL number
c
      if (mod(it-1,icfl).eq.0) then
         if (pert) then
            call boxcflBF(cflp_bf,bu1,bu2,yb,deta,xl,zl)
         else
            cflp_bf=0.
         end if
         if (.not.lin) then 
            call boxcfl(cflp,u2r,u2i,yb,deta,xl,zl,wbci)
         else
            cflp=0.
         end if
         cflp=cflp+cflp_bf
      end if
c
c     Compute extrema
c
      if (mod(it-1,iext).eq.0.and.iext.gt.0) call
     &     boxext(vext,cext,u2r,u2i,om2r,om2i,xs,yb,xl,zl)
c      
c     Calculate the full advection term (rotational form)
c     and save it in om2r,om2i
c     in physical space with or without dealiasing
c     Rotational form: H_i = eps_ijk * u_j * omega_k
c     Convective form: N_i = u_j du_i/dx_j
c     Relation:        N_i = .5*d(u_ju_j)/dx_i - H_i
c
      if (.true.) then
         if (.not.pert) then
c
c     Nonlinear term for the velocities
c
            do z=1,nzd_new/nprocz
               do xy=1,nxy
                  h1u = u2r(xy,z,2)*(om2r(xy,z,3)+2.*rot)-
     &                 u2r(xy,z,3)*om2r(xy,z,2)
                  h2u = u2r(xy,z,3)*om2r(xy,z,1)-
     &                 u2r(xy,z,1)*(om2r(xy,z,3)+2.*rot)
                  om2r(xy,z,3) = u2r(xy,z,1)*om2r(xy,z,2)-
     &                 u2r(xy,z,2)*om2r(xy,z,1)
                  om2r(xy,z,1) = h1u
                  om2r(xy,z,2) = h2u
                  h1e = u2i(xy,z,2)*(om2i(xy,z,3)+2.*rot)
     &                 -u2i(xy,z,3)*om2i(xy,z,2)
                  h2e = u2i(xy,z,3)*om2i(xy,z,1)
     &                 -u2i(xy,z,1)*(om2i(xy,z,3)+2.*rot)
                  om2i(xy,z,3) = u2i(xy,z,1)*om2i(xy,z,2)
     &                 -u2i(xy,z,2)*om2i(xy,z,1)
                  om2i(xy,z,1) = h1e
                  om2i(xy,z,2) = h2e
               end do
            end do           
c     
c     Nonlinear term for the scalar
c
            do ith=1,scalar
c               ri(ith) = 0.15
c               ri(ith) = 0.71/8.*1800.
c               ri(ith) = 0.
               do z=1,nzd_new/nprocz
                  do xy=1,nxy
c
c     This is for an active scalar
c
c                     if (ri(ith).ne.0.) then
c                        om2r(xy,z,2) = om2r(xy,z,2) + 
c     &                       th2r(xy,z,1+4*(ith-1))*ri(ith)
c                        om2i(xy,z,2) = om2i(xy,z,2) +
c     &                       th2i(xy,z,1+4*(ith-1))*ri(ith)
c                     end if

                     h1u= u2r(xy,z,1) * th2r(xy,z,2+4*(ith-1)) +
     &                    u2r(xy,z,2) * th2r(xy,z,3+4*(ith-1)) +
     &                    u2r(xy,z,3) * th2r(xy,z,4+4*(ith-1))
                     
                     h2u= u2i(xy,z,1) * th2i(xy,z,2+4*(ith-1)) +
     &                    u2i(xy,z,2) * th2i(xy,z,3+4*(ith-1)) +
     &                    u2i(xy,z,3) * th2i(xy,z,4+4*(ith-1))
                     
                     th2r(xy,z,4+4*(ith-1)) = -h1u
                     th2i(xy,z,4+4*(ith-1)) = -h2u
                  end do
               end do
            end do
         else
            if (scalar.gt.0) then
c
c     Perturbation formulation for scalar not implemented
c
               call stopnow(45353)
            end if

            if (.not.bf3) then
               if (lin) then
                  if (.not.adjoint) then
c
c     Perturbation formulation, linear, 2D base flow
c
                     do z=1,nzd_new/nprocz
                        do xy=1,nxy
                           h1u=     u2r(xy,z,2)*(bom1(xy,yb,3)+2.*rot)-
     &                          u2r(xy,z,3)*bom1(xy,yb,2)+
     &                          bu1(xy,yb,2)*(om2r(xy,z,3)+2.*rot)-
     &                          bu1(xy,yb,3)*om2r(xy,z,2)
                           
                           h2u=     u2r(xy,z,3)*bom1(xy,yb,1)-
     &                          u2r(xy,z,1)*(bom1(xy,yb,3)+2.*rot)+
     &                          bu1(xy,yb,3)*om2r(xy,z,1)-
     &                          bu1(xy,yb,1)*(om2r(xy,z,3)+2.*rot)
                           
                           om2r(xy,z,3)=u2r(xy,z,1)*bom1(xy,yb,2)-
     &                          u2r(xy,z,2)*bom1(xy,yb,1)+
     &                          bu1(xy,yb,1)*om2r(xy,z,2)-
     &                          bu1(xy,yb,2)*om2r(xy,z,1)
                           
                           om2r(xy,z,1)=h1u
                           om2r(xy,z,2)=h2u
                           
                           h1e=     u2i(xy,z,2)*(bom2(xy,yb,3)+2.*rot)-
     &                          u2i(xy,z,3)*bom2(xy,yb,2)+
     &                          bu2(xy,yb,2)*(om2i(xy,z,3)+2.*rot)-
     &                          bu2(xy,yb,3)*om2i(xy,z,2)
                           
                           h2e=     u2i(xy,z,3)*bom2(xy,yb,1)-
     &                          u2i(xy,z,1)*(bom2(xy,yb,3)+2.*rot)+
     &                          bu2(xy,yb,3)*om2i(xy,z,1)-
     &                          bu2(xy,yb,1)*(om2i(xy,z,3)+2.*rot)
                           om2i(xy,z,3)=u2i(xy,z,1)*bom2(xy,yb,2)-
     &                          u2i(xy,z,2)*bom2(xy,yb,1)+
     &                          bu2(xy,yb,1)*om2i(xy,z,2)-
     &                          bu2(xy,yb,2)*om2i(xy,z,1)
                           
                           om2i(xy,z,1) = h1e
                           om2i(xy,z,2) = h2e
                        end do
                     end do
                  else
c
c     Adjoint equations for 2-D base flow would go here
c
                  end if
               else

c     
c     Perturbation formulation, nonlinear, 2D base flow
c     
                  do z=1,nzd_new/nprocz
                     do xy=1,nxy
                        h1u=     u2r(xy,z,2)*(bom1(xy,yb,3)+2.*rot)-
     &                       u2r(xy,z,3)*bom1(xy,yb,2)+
     &                       bu1(xy,yb,2)*(om2r(xy,z,3)+2.*rot)-
     &                       bu1(xy,yb,3)*om2r(xy,z,2)+
     &                       u2r(xy,z,2)*(om2r(xy,z,3)+2.*rot)-
     &                       u2r(xy,z,3)*om2r(xy,z,2)
                           
                        h2u=     u2r(xy,z,3)*bom1(xy,yb,1)-
     &                       u2r(xy,z,1)*(bom1(xy,yb,3)+2.*rot)+
     &                       bu1(xy,yb,3)*om2r(xy,z,1)-
     &                       bu1(xy,yb,1)*(om2r(xy,z,3)+2.*rot)+
     &                       u2r(xy,z,3)*om2r(xy,z,1)-
     &                       u2r(xy,z,1)*(om2r(xy,z,3)+2.*rot)
                           
                        om2r(xy,z,3)=u2r(xy,z,1)*bom1(xy,yb,2)-
     &                       u2r(xy,z,2)*bom1(xy,yb,1)+
     &                       bu1(xy,yb,1)*om2r(xy,z,2)-
     &                       bu1(xy,yb,2)*om2r(xy,z,1)+
     &                       u2r(xy,z,1)*om2r(xy,z,2)-
     &                       u2r(xy,z,2)*om2r(xy,z,1)
                           
                        om2r(xy,z,1)=h1u
                        om2r(xy,z,2)=h2u
                           
                        h1e=     u2i(xy,z,2)*(bom2(xy,yb,3)+2.*rot)-
     &                       u2i(xy,z,3)*bom2(xy,yb,2)+
     &                       bu2(xy,yb,2)*(om2i(xy,z,3)+2.*rot)-
     &                       bu2(xy,yb,3)*om2i(xy,z,2)+
     &                       u2i(xy,z,2)*(om2i(xy,z,3)+2.*rot)-
     &                       u2i(xy,z,3)*om2i(xy,z,2)
                           
                        h2e=     u2i(xy,z,3)*bom2(xy,yb,1)-
     &                       u2i(xy,z,1)*(bom2(xy,yb,3)+2.*rot)+
     &                       bu2(xy,yb,3)*om2i(xy,z,1)-
     &                       bu2(xy,yb,1)*(om2i(xy,z,3)+2.*rot)+
     &                       u2i(xy,z,3)*om2i(xy,z,1)-
     &                       u2i(xy,z,1)*(om2i(xy,z,3)+2.*rot)
                        om2i(xy,z,3)=u2i(xy,z,1)*bom2(xy,yb,2)-
     &                       u2i(xy,z,2)*bom2(xy,yb,1)+
     &                       bu2(xy,yb,1)*om2i(xy,z,2)-
     &                       bu2(xy,yb,2)*om2i(xy,z,1)+
     &                       u2i(xy,z,1)*om2i(xy,z,2)-
     &                       u2i(xy,z,2)*om2i(xy,z,1)
                           
                        om2i(xy,z,1) = h1e
                        om2i(xy,z,2) = h2e
                     end do
                  end do                    
               end if
            else
               if (lin) then
                  if (.not. adjoint) then
c     
c     Perturbation, linearized, 3D base flow
c     
                     
                     do z=1,nzd_new/nprocz
                        do xy=1,nxy
                           h1u=   u2r(xy,z,2)*(u2bf3_r(xy,z,ybp,6)
     &                        +2.*rot)-
     &                        u2r(xy,z,3)*u2bf3_r(xy,z,ybp,5)+
     &                        u2bf3_r(xy,z,ybp,2)*(om2r(xy,z,3)+2.*rot)-
     &                        u2bf3_r(xy,z,ybp,3)*om2r(xy,z,2)
                           
                           h2u=   u2r(xy,z,3)*u2bf3_r(xy,z,ybp,4)-
     &                        u2r(xy,z,1)*
     &                        (u2bf3_r(xy,z,ybp,6)+2.*rot)+
     &                        u2bf3_r(xy,z,ybp,3)*om2r(xy,z,1)-
     &                        u2bf3_r(xy,z,ybp,1)*(om2r(xy,z,3)+2.*rot)
                           
                           om2r(xy,z,3)=u2r(xy,z,1)*u2bf3_r(xy,z,ybp,5)-
     &                          u2r(xy,z,2)*u2bf3_r(xy,z,ybp,4)+
     &                          u2bf3_r(xy,z,ybp,1)*om2r(xy,z,2)-
     &                          u2bf3_r(xy,z,ybp,2)*om2r(xy,z,1)
                           
                           
                           om2r(xy,z,1)=h1u
                           om2r(xy,z,2)=h2u
                           
                           h1e=     u2i(xy,z,2)*(u2bf3_i(xy,z,ybp,6)
     &                        +2.*rot)-
     &                        u2i(xy,z,3)*u2bf3_i(xy,z,ybp,5)+
     &                        u2bf3_i(xy,z,ybp,2)*(om2i(xy,z,3)+2.*rot)-
     &                        u2bf3_i(xy,z,ybp,3)*om2i(xy,z,2)
                           
                           h2e=     u2i(xy,z,3)*u2bf3_i(xy,z,ybp,4)-
     &                        u2i(xy,z,1)*
     &                        (u2bf3_i(xy,z,ybp,6)+2.*rot)+
     &                        u2bf3_i(xy,z,ybp,3)*om2i(xy,z,1)-
     &                        u2bf3_i(xy,z,ybp,1)*(om2i(xy,z,3)+2.*rot)
                           
                           om2i(xy,z,3)=u2i(xy,z,1)*u2bf3_i(xy,z,ybp,5)-
     &                          u2i(xy,z,2)*u2bf3_i(xy,z,ybp,4)+
     &                          u2bf3_i(xy,z,ybp,1)*om2i(xy,z,2)-
     &                          u2bf3_i(xy,z,ybp,2)*om2i(xy,z,1)
                           
                           om2i(xy,z,1) = h1e
                           om2i(xy,z,2) = h2e
                        end do
                     end do
                     
                  else
c
c     Adjoint equations for 3D base flow
c
                     do z=1,nzd_new/nprocz
                        do xy=1,nxy
                           
                           h1u= 2.*(u2bf3_r(xy,z,ybp,1)*u2xr(xy,z,1)+
     &                          u2bf3_r(xy,z,ybp,2)*u2xr(xy,z,2)+
     &                          u2bf3_r(xy,z,ybp,3)*u2xr(xy,z,3))-
     &                          (u2bf3_r(xy,z,ybp,2)*om2r(xy,z,3)-
     &                          u2bf3_r(xy,z,ybp,3)*om2r(xy,z,2))
                           
                           h2u= 2.*(u2bf3_r(xy,z,ybp,1)
     &                          *(u2xr(xy,z,2)-om2r(xy,z,3))+
     &                          u2bf3_r(xy,z,ybp,2)*(-u2xr(xy,z,1)
     &                          -u2zr(xy,z,3))+u2bf3_r(xy,z,ybp,3)
     &                          *(u2zr(xy,z,2)+om2r(xy,z,1)))-
     &                          (u2bf3_r(xy,z,ybp,3)*om2r(xy,z,1)-
     &                          u2bf3_r(xy,z,ybp,1)*om2r(xy,z,3))
               
                           om2r(xy,z,3)= 2.*(u2bf3_r(xy,z,ybp,1)*
     &                          u2zr(xy,z,1)+
     &                          u2bf3_r(xy,z,ybp,2)*u2zr(xy,z,2)+
     &                          u2bf3_r(xy,z,ybp,3)*u2zr(xy,z,3))-
     &                          (u2bf3_r(xy,z,ybp,1)*om2r(xy,z,2)-
     &                          u2bf3_r(xy,z,ybp,2)*om2r(xy,z,1))
                           
                           om2r(xy,z,1)=h1u
                           om2r(xy,z,2)=h2u
                           
                           h1e= 2.*(u2bf3_i(xy,z,ybp,1)*u2xi(xy,z,1)+
     &                          u2bf3_i(xy,z,ybp,2)*u2xi(xy,z,2)+
     &                          u2bf3_i(xy,z,ybp,3)*u2xi(xy,z,3))-
     &                          (u2bf3_i(xy,z,ybp,2)*om2i(xy,z,3)-
     &                          u2bf3_i(xy,z,ybp,3)*om2i(xy,z,2))
                           
                           h2e= 2.*(u2bf3_i(xy,z,ybp,1)
     &                          *(u2xi(xy,z,2)-om2i(xy,z,3))+
     &                          u2bf3_i(xy,z,ybp,2)
     &                          *(-u2xi(xy,z,1)-u2zi(xy,z,3))+
     &                          u2bf3_i(xy,z,ybp,3)
     &                          *(u2zi(xy,z,2)+om2i(xy,z,1)))-
     &                          (u2bf3_i(xy,z,ybp,3)*om2i(xy,z,1)-
     &                          u2bf3_i(xy,z,ybp,1)*om2i(xy,z,3))
                           
                           om2i(xy,z,3)= 2.*(u2bf3_i(xy,z,ybp,1)*
     &                          u2zi(xy,z,1)+
     &                          u2bf3_i(xy,z,ybp,2)*u2zi(xy,z,2)+
     &                          u2bf3_i(xy,z,ybp,3)*u2zi(xy,z,3))-
     &                          (u2bf3_i(xy,z,ybp,1)*om2i(xy,z,2)-
     &                          u2bf3_i(xy,z,ybp,2)*om2i(xy,z,1))
                           
                           om2i(xy,z,1) = h1e
                           om2i(xy,z,2) = h2e
                        end do
                     end do
                  end if
c
c     Adjoint equations for 3D base flow end here
c
               else
                     
c     Perturbation, nonlinear, 3D base flow
c     
                  do z=1,nzd_new/nprocz
                     do xy=1,nxy
                        h1u=     u2r(xy,z,2)*(u2bf3_r(xy,z,ybp,6)
     &                       +2.*rot)-
     &                       u2r(xy,z,3)*u2bf3_r(xy,z,ybp,5)+
     &                       u2bf3_r(xy,z,ybp,2)*(om2r(xy,z,3)+2.*rot)-
     &                       u2bf3_r(xy,z,ybp,3)*om2r(xy,z,2)+
     &                       u2r(xy,z,2)*(om2r(xy,z,3)+2.*rot)-
     &                       u2r(xy,z,3)*om2r(xy,z,2)
                           
                        h2u=     u2r(xy,z,3)*u2bf3_r(xy,z,ybp,4)-
     &                       u2r(xy,z,1)*
     &                       (u2bf3_r(xy,z,ybp,6)+2.*rot)+
     &                       u2bf3_r(xy,z,ybp,3)*om2r(xy,z,1)-
     &                       u2bf3_r(xy,z,ybp,1)*(om2r(xy,z,3)+2.*rot)+
     &                       u2r(xy,z,3)*om2r(xy,z,1)-
     &                       u2r(xy,z,1)*(om2r(xy,z,3)+2.*rot)
                           
                        om2r(xy,z,3)=u2r(xy,z,1)*u2bf3_r(xy,z,ybp,5)-
     &                       u2r(xy,z,2)*u2bf3_r(xy,z,ybp,4)+
     &                       u2bf3_r(xy,z,ybp,1)*om2r(xy,z,2)-
     &                       u2bf3_r(xy,z,ybp,2)*om2r(xy,z,1)+
     &                       u2r(xy,z,1)*om2r(xy,z,2)-
     &                       u2r(xy,z,2)*om2r(xy,z,1)
                           
                        om2r(xy,z,1)=h1u
                        om2r(xy,z,2)=h2u
                           
                        h1e=     u2i(xy,z,2)*(u2bf3_i(xy,z,ybp,6)
     &                       +2.*rot)-
     &                       u2i(xy,z,3)*u2bf3_i(xy,z,ybp,5)+
     &                       u2bf3_i(xy,z,ybp,2)*
     &                       (om2i(xy,z,3)+2.*rot)-
     &                       u2bf3_i(xy,z,ybp,3)*om2i(xy,z,2)+
     &                       u2i(xy,z,2)*(om2i(xy,z,3)+2.*rot)-
     &                       u2i(xy,z,3)*om2i(xy,z,2)
                           
                        h2e=     u2i(xy,z,3)*u2bf3_i(xy,z,ybp,4)-
     &                       u2i(xy,z,1)*
     &                       (u2bf3_i(xy,z,ybp,6)+2.*rot)+
     &                       u2bf3_i(xy,z,ybp,3)*om2i(xy,z,1)-
     &                       u2bf3_i(xy,z,ybp,1)*
     &                       (om2i(xy,z,3)+2.*rot)+
     &                       u2i(xy,z,3)*om2i(xy,z,1)-
     &                       u2i(xy,z,1)*(om2i(xy,z,3)+2.*rot)
                        om2i(xy,z,3)=u2i(xy,z,1)*u2bf3_i(xy,z,ybp,5)-
     &                       u2i(xy,z,2)*u2bf3_i(xy,z,ybp,4)+
     &                       u2bf3_i(xy,z,ybp,1)*om2i(xy,z,2)-
     &                       u2bf3_i(xy,z,ybp,2)*om2i(xy,z,1)+
     &                       u2i(xy,z,1)*om2i(xy,z,2)-
     &                       u2i(xy,z,2)*om2i(xy,z,1)
                           
                        om2i(xy,z,1) = h1e
                        om2i(xy,z,2) = h2e
                     end do
                  end do
               end if
            end if
         end if
      else
c     
c     Suction boundary layer
c     
         if (pert) then
c     
c     Perturbation formulation for suction not implemented
c
            call stopnow(534523)
         end if
               write(*,*)'Not Implemented'
               call stopnow(4455667)
               do z=1,nzpc
            do xy=1,nxy
               h1u = (u2r(xy,z,2)-vsuc)*(om2r(xy,z,3)+2.*rot)-
     &              u2r(xy,z,3)*om2r(xy,z,2)
               h2u = u2r(xy,z,3)*om2r(xy,z,1)-
     &              u2r(xy,z,1)*(om2r(xy,z,3)+2.*rot)
               om2r(xy,z,3) = u2r(xy,z,1)*om2r(xy,z,2)-
     &              (u2r(xy,z,2)-vsuc)*om2r(xy,z,1)
               om2r(xy,z,1) = h1u
               om2r(xy,z,2) = h2u
               h1e = (u2i(xy,z,2)-vsuc)*(om2i(xy,z,3)+2.*rot)
     &              -u2i(xy,z,3)*om2i(xy,z,2)
               h2e = u2i(xy,z,3)*om2i(xy,z,1)
     &              -u2i(xy,z,1)*(om2i(xy,z,3)+2.*rot)
               om2i(xy,z,3) = u2i(xy,z,1)*om2i(xy,z,2)
     &              -(u2i(xy,z,2)-vsuc)*om2i(xy,z,1)
               om2i(xy,z,1) = h1e
               om2i(xy,z,2) = h2e
            end do
         end do
         do ith=1,scalar
            do z=1,nzpc
               do xy=1,nxy
                  h1u= u2r(xy,z,1)*th2r(xy,z,2+4*(ith-1))+
     &                 (u2r(xy,z,2)-vsuc)*th2r(xy,z,3+4*(ith-1))+
     &                 u2r(xy,z,3)*th2r(xy,z,4+4*(ith-1))
                  
                  h2u= u2i(xy,z,1)*th2i(xy,z,2+4*(ith-1))+
     &                 (u2i(xy,z,2)-vsuc)*th2i(xy,z,3+4*(ith-1))+
     &                 u2i(xy,z,3)*th2i(xy,z,4+4*(ith-1))

                  th2r(xy,z,4+4*(ith-1)) = -h1u
                  th2i(xy,z,4+4*(ith-1)) = -h2u
               end do
            end do
         end do
      end if
c
c     Fringe region for spatial simulations (or temporal simulations
c     with fringe, like fltype 9)
c
      if (spat) then
         call fring(om2r,om2i,u2r,u2i,tc,xsc,zsc,xl,zl,yb,
     &        fstart,fend,bu1,bu2,
     &        osmod,osnumb,
     &        osur,osui,osvr,osvi,oswr,oswi,
     &        evr,evi,eur,eui,ewr,ewi,afw,alfaw,betaw,ampob,amp2d,
     &        ev2r,ev2i,eu2r,eu2i,afw2d,alf2d,
     &        xc1,xc2,fring1,fring2,cphas1,cphas2,sphas1,sphas2,
     &        th2r,th2i,     
     &        ampst,streak,betast,omegast,
     &        ndxst,uust_r,uust_i,vvst_r,vvst_i,wwst_r,wwst_i,
     &        tsmoo,tsmst,tsmend,iampst,phist,pert,vsuc,my_node_world)
      end if
c
c     Localized volume force
c
      if (loctyp.ge.1) then
         call locf(om2r,om2i,yb,xl,zl,xsc,zsc,eta,tc,loctyp,
     &        fp1,fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7,fpds8,
     &        fpdds1,fpdds2,fpdds3,fpdds4,fpdds5,wr,wi,th2r,th2i,
     &        my_node_world)
      end if
c
c     Trip forcing (on the velocities)
c
      if (tripf) then 
         call trip(om2r,om2i,yb,xl,xsc,eta,txsc,tysc,ttzc,tx0,
     &        my_node_world)
      end if

      if (nproc.eq.1) then
c
c     All code for nproc = 1 removed, it can be found in older versions!!!
c
      end if

c
c     Adding a volume force to wavenumber zero (mean component)
c     to get the correct boundary layer growth for parallel flows.
c     Note that we are in Fourier-Chebyshev space but have not yet
c     normalized for the Fourier transform.
c     This means that the force has to be multiplied by nxp*nzp.
c     The time-independent forcing term (f1,0,f3) on the right hand side
c     of the Navier-Stokes equations needs to preserve the initial
c     parallel profile if no disturbances are present.
c     f1=-c*dU/dx-1/Re*d2U/d2y=-f'''/(2*x0)
c     f3=-spanv/Re*d2W/d2y
c     f is similarity solution of
c     f'''+f*f''+2*m/(m+1)*(1-f'*F')=0
c     g''+fg'=0
c     where
c     f'=df/ds
c     s=y*sqrt((m+1)*U_0/(2*nu*x0)), x0 inflow location from the leading edge,
c     scaled by displacement thickness, c a reference
c     speed (usually a group velocity). These terms are easily computed with
c     the nonlinear terms over here. There is other choice, which is time-
c     dependent forcing.(See Spalart & Yang, 1987, JFM) 
c
      if (.true.) then
c
c     Forcing only if Blasius base flow is considered
c
         if (fltype.eq. 3.or.fltype.eq.9.or.
     &       fltype.eq.-1.or.fltype.eq.-2) then
               write(*,*)'Not Implemented'
               call stopnow(5566778)
            c1=real(nxp)*real(nzp)
 1          xx=(xblc/x0)**rlam
            do y=1,npl
               ybl=1.+eta(y+yb-1)
               etabl=ybl*sqrt((rlam+1.)*re*xx/(2.*xblc))
               h1u = cdev*rlam*xx/xblc
               h2u = 0.5*(rlam-1)*sqrt(0.5*(rlam+1.)*re)*
     &               ybl*cdev*(xx/xblc)**1.5
               h3u = 0.5*(rlam+1.)/xblc
               d1f=cubip(etabl,fbla(1,2),dybla,nbla)
               d2f=cubip(etabl,fbla(1,3),dybla,nbla)
               d3f=cubip(etabl,fbla(1,4),dybla,nbla)
               om2r(1+(y-1)*(nxp/2+1),1,1)=om2r(1+(y-1)*(nxp/2+1),1,1)-
     &              (d1f*h1u + d2f*h2u + d3f*xx**2*h3u)*c1
               if (fltype.eq.-2.or.fltype.eq.9) then
                  d2g=cubip(etabl,fbla(1,7),dybla,nbla)
                  om2r(1+(y-1)*(nxp/2+1),1,3)=om2r(1+(y-1)*
     &                 (nxp/2+1),1,3)-spanv*h3u*d2g*c1
               end if
            end do
         end if
      end if
c
c     Additional forcing terms
c     Note that these could also be added in the linear step since
c     they do not involve nonlinear products. Thereby, a communication
c     step could be avoided making everything a bit more efficient.
c


#ifdef MPI
      do i=1,3
c
c     Forward Fourier transform of nonlinear term H_i in x
c     and put it onto boxr_z
c
         call vrfftf(om2r(1,1,i),om2i(1,1,i),wr_x,wi_x,
     &        nxp,nzd_new/nprocz,1,nxp/2+1,prex)

c         call mpi_barrier(mpi_comm_world,ierror)
         call wall_time(t2)
         tser = tser + (t2-t1)
         t1=t2

         call putpxz_x(om2r(1,1,i),om2i(1,1,i),yb,i,1,
     &        boxr_z(1,1,i),boxi_z(1,1,i),
     &        realg3,realg4,my_node_x,my_comm_x,my_node_world)
         
c         call mpi_barrier(mpi_comm_world,ierror)
         call wall_time(t2)
         tcom = tcom + (t2-t1)
         t1=t2

c
c     Forward Fourier transform of nonlinear term H_i in z
c
         call vcfftf(boxr_z(1,1,i),boxi_z(1,1,i),wr_z,wi_z,
     &        nzp,memnx,memnx,1,prez)
         
c         call mpi_barrier(mpi_comm_world,ierror)
         call wall_time(t2)
         tser = tser + (t2-t1)
         t1=t2

c            call putpxz_z(boxr_z(1,1,i),boxi_z(1,1,i),yb,i,1,ur,ui,
c     &           realg1,realg2,my_node_z,my_comm_z,my_node_world)

c         call mpi_barrier(mpi_comm_world,ierror)
         call wall_time(t2)
         tcom = tcom + (t2-t1)
         t1=t2

      end do
      do ith=1,scalar
c     
c     Forward Fourier transform of nonlinear term H_i in x
c
         call vrfftf(th2r(1,1,4+4*(ith-1)),th2i(1,1,4+4*(ith-1)),
     &        wr_x,wi_x,nxp,nzd_new/nprocz,1,nxp/2+1,prex)
         
c         call mpi_barrier(mpi_comm_world,ierror)
         call wall_time(t2)
         tser = tser + (t2-t1)
         t1=t2

         call putpxz_x(th2r(1,1,4+4*(ith-1)),th2i(1,1,4+4*(ith-1)),yb
     &        ,8+pressure+3*(ith-1),1,
     &        boxr(1,1,4+4*(ith-1)),boxi(1,1,4+4*(ith-1)),
     &        realg3,realg4,my_node_x,my_comm_x,my_node_world)
         
c         call mpi_barrier(mpi_comm_world,ierror)
         call wall_time(t2)
         tcom = tcom + (t2-t1)
         t1=t2

c
c     Forward Fourier transform of nonlinear term H_i in z
c
         call vcfftf(boxr(1,1,4+4*(ith-1)),boxi(1,1,4+4*(ith-1)),
     &        wr_z,wi_z,nzp,memnx,memnx,1,prez)
         
c         call mpi_barrier(mpi_comm_world,ierror)
         call wall_time(t2)
         tser = tser + (t2-t1)
         t1=t2

c            call putpxz_z(boxr(1,1,4+4*(ith-1)),boxi(1,1,4+4*(ith-1)),yb
c     &           ,8+pressure+3*(ith-1),1,ur,ui,
c     &           realg1,realg2,my_node_z,my_comm_z,my_node_world)

c         call mpi_barrier(mpi_comm_world,ierror)
         call wall_time(t2)
         tcom = tcom + (t2-t1)
         t1=t2
         
      end do
#endif
      


c
c-----------------------------------------------------------------
c
c     MHD Lorentz force
c
      if (imhd.eq.1) then
               write(*,*)'Not Implemented'
               call stopnow(55667788)
c     
c     Reuse du2r array for the forcing term
c
         do ll=1,2
            if (nproc.eq.1) then
               call getxz(du2r(1,1,ll,1),du2i(1,1,ll,1),
     &              yb,ll,1,fmhdr,fmhdi)
            else
#ifdef MPI
               call getpxz(du2r(1,1,ll,1),du2i(1,1,ll,1),yb,ll,1,
     &              fmhdr,fmhdi,realg1,realg2,my_node_world)
#endif
            end if
         end do
c     
c     Add the MHD force to the nonlinear term
c     (note: Factors nxp*nzp are added due to non-normalised FFT)
c     
         fact = mhd_n*nxp*nzp
         do z=1,nz/2
            do y=yb,npp
               y1=(y-yb)*nxp2
               do x=1,nx/2
                  xy=x+y1
                  
                  cr1 =  alfa(x)*du2i(xy,z,1,1)
     &                 +du2r(xy,z,2,2)*b0(3)-du2r(xy,z,3,2)*b0(2)
                  ci1 = -alfa(x)*du2r(xy,z,1,1)
     &                 +du2i(xy,z,2,2)*b0(3)-du2i(xy,z,3,2)*b0(2)
                  
                  cr2 = -du2r(xy,z,2,1)
     &                 +du2r(xy,z,3,2)*b0(1)-du2r(xy,z,1,2)*b0(3)
                  ci2 = -du2i(xy,z,2,1)
     &                 +du2i(xy,z,3,2)*b0(1)-du2i(xy,z,1,2)*b0(3)
                  
                  cr3 =  beta(z)*du2i(xy,z,1,1)
     &                 +du2r(xy,z,1,2)*b0(2)-du2r(xy,z,2,2)*b0(1)
                  ci3 = -beta(z)*du2r(xy,z,1,1)
     &                 +du2i(xy,z,1,2)*b0(2)-du2i(xy,z,2,2)*b0(1)
c     
c     Lorentz force
c     
                  om2r(xy,z,1) = om2r(xy,z,1) +
     &                 fact*(cr2*b0(3) - cr3*b0(2))
                  om2i(xy,z,1) = om2i(xy,z,1) +
     &                 fact*(ci2*b0(3) - ci3*b0(2))
                  
                  om2r(xy,z,2) = om2r(xy,z,2) +
     &                 fact*(cr3*b0(1) - cr1*b0(3))
                  om2i(xy,z,2) = om2i(xy,z,2) +
     &                 fact*(ci3*b0(1) - ci1*b0(3))
                  
                  om2r(xy,z,3) = om2r(xy,z,3) +
     &                 fact*(cr1*b0(2) - cr2*b0(1))
                  om2i(xy,z,3) = om2i(xy,z,3) +
     &                 fact*(ci1*b0(2) - ci2*b0(1))
                  

                  end do
               end do
            end do
c     
c     NOTE: +1 could be replaced by +2 (since these are the oddball modes)
c     
            do z=nz/2+1,nz
               zp=nzp-nz+z
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1


                  cr1 =  alfa(x)*du2i(xy,zp,1,1)
     &                 +du2r(xy,zp,2,2)*b0(3)-du2r(xy,zp,3,2)*b0(2)
                  ci1 = -alfa(x)*du2r(xy,zp,1,1)
     &                 +du2i(xy,zp,2,2)*b0(3)-du2i(xy,zp,3,2)*b0(2)
                  
                  cr2 = -du2r(xy,zp,2,1)
     &                 +du2r(xy,zp,3,2)*b0(1)-du2r(xy,zp,1,2)*b0(3)
                  ci2 = -du2i(xy,zp,2,1)
     &                 +du2i(xy,zp,3,2)*b0(1)-du2i(xy,zp,1,2)*b0(3)
                  
                  cr3 =  beta(z)*du2i(xy,zp,1,1)
     &                 +du2r(xy,zp,1,2)*b0(2)-du2r(xy,zp,2,2)*b0(1)
                  ci3 = -beta(z)*du2r(xy,zp,1,1)
     &                 +du2i(xy,zp,1,2)*b0(2)-du2i(xy,zp,2,2)*b0(1)
c     
c     Lorentz force
c     
                  om2r(xy,zp,1) = om2r(xy,zp,1) +
     &                 fact*(cr2*b0(3) - cr3*b0(2))
                  om2i(xy,zp,1) = om2i(xy,zp,1) +
     &                 fact*(ci2*b0(3) - ci3*b0(2))
                  
                  om2r(xy,zp,2) = om2r(xy,zp,2) +
     &                 fact*(cr3*b0(1) - cr1*b0(3))
                  om2i(xy,zp,2) = om2i(xy,zp,2) +
     &                 fact*(ci3*b0(1) - ci1*b0(3))
                  
                  om2r(xy,zp,3) = om2r(xy,zp,3) +
     &                 fact*(cr1*b0(2) - cr2*b0(1))
                  om2i(xy,zp,3) = om2i(xy,zp,3) +
     &                 fact*(ci1*b0(2) - ci2*b0(1))

                  end do
               end do
            end do

      end if

c     
c
c--------------------------------------------------------------
c
c     Large-eddy simulation ( LES )
c
      if (iles.eq.1) then
c
c     ADM-RT model:
c     Get filtered velocities
c
         do ll=1,3
            if (nproc.eq.1) then
               call getxz(u2r(1,1,1),u2i(1,1,1),yb,ll,1,gur,gui)
            else
#ifdef MPI
               call getpxz_z(boxr_z(1,1,3+ll),boxi_z(1,1,3+ll),yb,ll,1,
     &              gur,gui,realg1,realg2,my_node_z,my_comm_z,
     &              my_node_world)

               call getpxz_x(om2r(1,1,ll),om2i(1,1,ll),yb,ll+3,1,
     &              boxr_z(1,1,ll),boxi_z(1,1,ll),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
               
               call getpxz_x(u2r(1,1,ll),u2i(1,1,ll),yb,ll,1,
     &              boxr_z(1,1,3+ll),boxi_z(1,1,3+ll),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
#endif
            end if
c
c     Add the relaxation term to the nonlinear term
c     (note: Factors nxp*nzp are added due to non-normalised FFT)
c
            if (nproc.eq.1) then
c
c     Code removed for nproc = 1, see older versions if you need it!!!
c
            else
#ifdef MPI
               do z=1,nzd_new/nprocz
                  do y=yb,npp
                     y1=(y-yb)*nxp2
                     do x=1,nx/2
                        xy=x+y1
                        om2r(xy,z,ll) = om2r(xy,z,ll) - 
     &                       chi_scaled*u2r(xy,z,ll)*nxp*nzp
                        om2i(xy,z,ll) = om2i(xy,z,ll) - 
     &                       chi_scaled*u2i(xy,z,ll)*nxp*nzp
                     end do
                  end do
               end do

#endif
            end if
            
         end do
         
         do ith=1,scalar
            if (nproc.eq.1) then
               call getxz(u2r(1,1,1),u2i(1,1,1),yb,ith,1,gsr,gsi)
            else
#ifdef MPI
               call getpxz_z(boxr(1,1,1+4*(ith-1)),
     &              boxi(1,1,1+4*(ith-1)),yb,ith,1,gsr,gsi,
     &              realg1,realg2,my_node_z,my_comm_z,my_node_world)

               call getpxz_x(u2r(1,1,1),u2i(1,1,1),yb,ith,1,
     &              boxr(1,1,1+4*(ith-1)),boxi(1,1,1+4*(ith-1)),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)

               call getpxz_x(th2r(1,1,4+4*(ith-1)),th2i(1,1,4+4*(ith-1))
     &              ,yb,8+pressure+3*(ith-1),1,
     &              boxr(1,1,4+4*(ith-1)),boxi(1,1,4+4*(ith-1)),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
#endif
            end if
            if (nproc.eq.1) then
c
c     Code for nproc = 1 removed, see older versions if you need it!!!
c
            else
#ifdef MPI
               do z=1,nzd_new/nprocz
                  do y=yb,npp
                     y1=(y-yb)*nxp2
                     do x=1,nx/2
                        xy=x+y1    
                        th2r(xy,z,4*(ith-1)+4) = th2r(xy,z,4*(ith-1)+4) 
     &                       - chi_scaled*u2r(xy,z,1)*nxp*nzp/prt
                        th2i(xy,z,4*(ith-1)+4) = th2i(xy,z,4*(ith-1)+4) 
     &                       - chi_scaled*u2i(xy,z,1)*nxp*nzp/prt
                     end do
                  end do
               end do
#endif
            end if
         end do
      else if (iles.eq.2.or.iles.eq.3) then
c
c     (HPF) Eddy-viscosity model:
c     get the SGS force
c
         do ll=1,3
            if (nproc.eq.1) then
               call getxz(u2r(1,1,1),u2i(1,1,1),yb,ll,1,gur,gui)
            else
#ifdef MPI
               call getpxz_z(boxr_z(1,1,3+ll),boxi_z(1,1,3+ll),yb,ll,1,
     &              gur,gui,realg1,realg2,my_node_z,my_comm_z,
     &              my_node_world)

               call getpxz_x(u2r(1,1,ll),u2i(1,1,ll),yb,ll,1,
     &              boxr_z(1,1,3+ll),boxi_z(1,1,3+ll),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)

               call getpxz_x(om2r(1,1,ll),om2i(1,1,ll),yb,ll,1,
     &              boxr_z(1,1,ll),boxi_z(1,1,ll),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
#endif
            end if

            if (nproc.eq.1) then
c
c     Code removed for nproc = 1, see older versions if you need it!!!
c
            else
#ifdef MPI
               do z=1,nzd_new/nprocz
                  do y=yb,npp
                     y1=(y-yb)*nxp2
                     do x=1,nx/2
                        xy=x+y1
                        om2r(xy,z,ll) = om2r(xy,z,ll) - 
     &                       u2r(xy,z,ll)*nxp*nzp
                        om2i(xy,z,ll) = om2i(xy,z,ll) - 
     &                       u2i(xy,z,ll)*nxp*nzp
                     end do
                  end do
               end do
#endif
            end if
         end do
c
c     (HPF) Eddy-diffusivity model for the scalars
c
         do ith=1,scalar

            if (nproc.eq.1) then
               call getxz(u2r(1,1,1),u2i(1,1,1),yb,ith,1,gsr,gsi)
            else
#ifdef MPI
               call getpxz_z(boxr(1,1,1+4*(ith-1)),
     &              boxi(1,1,1+4*(ith-1)),yb,ith,1,gsr,gsi,
     &              realg1,realg2,my_node_z,my_comm_z,my_node_world)

               call getpxz_x(u2r(1,1,1),u2i(1,1,1),yb,ith,1,
     &              boxr(1,1,1+4*(ith-1)),boxi(1,1,1+4*(ith-1)),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)

               call getpxz_x(th2r(1,1,4+4*(ith-1)),th2i(1,1,4+4*(ith-1))
     &              ,yb,8+pressure+3*(ith-1),1,
     &              boxr(1,1,4+4*(ith-1)),boxi(1,1,4+4*(ith-1)),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
#endif
            end if
            if (nproc.eq.1) then
c
c     Code removed for nproc = 1, see older versions if you need it!
c
            else
#ifdef MPI
               do z=1,nzd_new/nprocz
                  do y=yb,npp
                     y1=(y-yb)*nxp2
                     do x=1,nx/2
                        xy=x+y1    
                        th2r(xy,z,4*(ith-1)+4) = th2r(xy,z,4*(ith-1)+4) 
     &                       - u2r(xy,z,1)*nxp*nzp
                        th2i(xy,z,4*(ith-1)+4) = th2i(xy,z,4*(ith-1)+4) 
     &                       - u2i(xy,z,1)*nxp*nzp
                     end do
                  end do
               end do
#endif
            end if
         end do
         
      end if
c
c     SFD (Selective Frequency Damping)
c
c     Memory layout:
c     sfdur(:,:,:,1-3): temporally filtered velocity field
c     sfdur(:,:,:,4-6): unfiltered-filtered velocity field (RHS^n-1)
c
c     u2r(:,:,1): temporally filtered velocity
c     u2r(:,:,2): unfiltered velocity
c     u2r(:,:,3): unfiltered-filtered velocity (RHS^n)
c     u2r(:,:,1): new filtered velocity
c     u2r(:,:,2): RHS^n-1
c
c     Get filtered and unfiltered velocities into u2r(:,:,1) and u2r(:,:,2)
c
      if (isfd.eq.1) then
         do ll=1,3
            if (nproc.eq.1) then
               call stopnow(54454)
               call getxz(u2r(1,1,1),u2i(1,1,1),yb,ll,1,sfdur,sfdui)
               call getxz(u2r(1,1,2),u2i(1,1,2),yb,ll,1,   ur,   ui)

            else
#ifdef MPI
c
c     Get filtered solution in u2r(:,:,1)
c
               call getpxz_z(boxr_z(1,1,4),boxi_z(1,1,4),yb,ll,1,
     &              sfdur,sfdui,realg1,realg2,my_node_z,my_comm_z,
     &              my_node_world)
               call getpxz_x(u2r(1,1,1),u2i(1,1,1),yb,4,1,
     &              boxr_z(1,1,4),boxi_z(1,1,4),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
c
c     Get unfiltered solution in u2r(:,:,2)
c
               call getpxz_z(boxr_z(1,1,4),boxi_z(1,1,4),yb,ll,1,
     &              ur,ui,realg1,realg2,my_node_z,my_comm_z,
     &              my_node_world)
               call getpxz_x(u2r(1,1,2),u2i(1,1,2),yb,4,1,
     &              boxr_z(1,1,4),boxi_z(1,1,4),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
c
c     Get nonlinear term in om2r
c
               call getpxz_x(om2r(1,1,ll),om2i(1,1,ll),yb,ll+3,1,
     &              boxr_z(1,1,ll),boxi_z(1,1,ll),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
#endif

            end if
c
c     Add the SFD relaxation term to the nonlinear term
c     (note: Factors nxp*nzp are added due to non-normalised FFT)
c
            do z=1,nzd_new/nprocz
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     u2r(xy,z,3) = u2r(xy,z,2)-u2r(xy,z,1)
                     u2i(xy,z,3) = u2i(xy,z,2)-u2i(xy,z,1)

                     om2r(xy,z,ll) = om2r(xy,z,ll) - 
     &                    sfd_chi*u2r(xy,z,3)*nxp*nzp
                     om2i(xy,z,ll) = om2i(xy,z,ll) - 
     &                    sfd_chi*u2i(xy,z,3)*nxp*nzp
                  end do
               end do
            end do
c
c     Do the RK integration (explicit) of the filtered term
c     Add G^n in u2r(:,:,1)
c
            do z=1,nzd_new/nprocz
               do y=yb,npp
                  y1=(y-yb)*nxp2
                  do x=1,nx/2
                     xy=x+y1
                     u2r(xy,z,1) = u2r(xy,z,1) +
     &                    an/sfd_delta*u2r(xy,z,3)
                     u2i(xy,z,1) = u2i(xy,z,1) +
     &                    an/sfd_delta*u2i(xy,z,3)
                  end do
               end do
            end do
            if (bn.ne.0) then
               ! get G^n-1 onto u2r(:,:,2)
               if (nproc.eq.1) then
                  call getxz(u2r(1,1,2),u2i(1,1,2),yb,ll+3,1,
     &                 sfdur,sfdui)
               else
#ifdef MPI
               call getpxz_z(boxr_z(1,1,4),boxi_z(1,1,4),yb,ll+3,1,
     &              sfdur,sfdui,realg1,realg2,my_node_z,my_comm_z,
     &              my_node_world)
               call getpxz_x(u2r(1,1,2),u2i(1,1,2),yb,4,1,
     &              boxr_z(1,1,4),boxi_z(1,1,4),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
#endif
               end if
               ! add G^n-1
               do z=1,nzd_new/nprocz
                  do y=yb,npp
                     y1=(y-yb)*nxp2
                     do x=1,nx/2
                        xy=x+y1
                        u2r(xy,z,1) = u2r(xy,z,1) +
     &                       bn/sfd_delta*u2r(xy,z,2)
                        u2i(xy,z,1) = u2i(xy,z,1) +
     &                       bn/sfd_delta*u2i(xy,z,2)
                     end do
                  end do
               end do
            end if
c
c     Put G^n and filtered u
c
            if (nproc.eq.1) then
               call putxz(u2r(1,1,1),u2i(1,1,1),yb,ll  ,sfdur,sfdui)
               call putxz(u2r(1,1,3),u2i(1,1,3),yb,ll+3,sfdur,sfdui)
            else
#ifdef MPI
               call putpxz_x(u2r(1,1,1),u2i(1,1,1),yb,1,1,
     &              boxr_z(1,1,4),boxi_z(1,1,4),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
               call putpxz_z(boxr_z(1,1,4),boxi_z(1,1,4),yb,ll,1,
     &              sfdur,sfdui,
     &              realg1,realg2,my_node_z,my_comm_z,my_node_world)
               call putpxz_x(u2r(1,1,3),u2i(1,1,3),yb,1,1,
     &              boxr_z(1,1,4),boxi_z(1,1,4),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
               call putpxz_z(boxr_z(1,1,4),boxi_z(1,1,4),yb,ll+3,1,
     &              sfdur,sfdui,
     &              realg1,realg2,my_node_z,my_comm_z,my_node_world)




#endif
            end if
         end do
      end if
c
c     Put the planes of the nonlinear term back onto the velocities
c     (here, the truncation to the normal grid occurs)
c     NOTE: THE ODDBALL MODES ARE STILL IN THERE!!!!
c



#ifdef MPI
c      call mpi_barrier(mpi_comm_world,ierror)
#endif
      call wall_time(t2)
      tser = tser + (t2-t1)
      t1=t2

      
      if (nproc.eq.1) then
c
c     Code removed for nproc = 1, see older versions if you need it!
c
      else
#ifdef MPI
         do i=1,3
c
c     Forward Fourier transform of nonlinear term H_i in x
c
c            call mpi_barrier(mpi_comm_world,ierror)
            call wall_time(t2)
            tser = tser + (t2-t1)
            t1=t2
            if (iles.ne.0.or.isfd.ne.0) then
               call putpxz_x(om2r(1,1,i),om2i(1,1,i),yb,i,1,
     &              boxr_z(1,1,i),boxi_z(1,1,i),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
            end if
c            call mpi_barrier(mpi_comm_world,ierror)
            call wall_time(t2)
            tcom = tcom + (t2-t1)
            t1=t2

c
c     Forward Fourier transform of nonlinear term H_i in z
c
c            call mpi_barrier(mpi_comm_world,ierror)
            call wall_time(t2)
            tser = tser + (t2-t1)
            t1=t2

            call putpxz_z(boxr_z(1,1,i),boxi_z(1,1,i),yb,i,1,ur,ui,
     &           realg1,realg2,my_node_z,my_comm_z,my_node_world)

c            call mpi_barrier(mpi_comm_world,ierror)
            call wall_time(t2)
            tcom = tcom + (t2-t1)
            t1=t2

         end do
         do ith=1,scalar
c
c     Forward Fourier transform of nonlinear term H_i in x
c
c            call mpi_barrier(mpi_comm_world,ierror)
            call wall_time(t2)
            tser = tser + (t2-t1)
            t1=t2
            if (iles.ne.0) then
               call putpxz_x(th2r(1,1,4+4*(ith-1)),th2i(1,1,4+4*(ith-1))
     &              ,yb,8+pressure+3*(ith-1),1,
     &              boxr(1,1,4+4*(ith-1)),boxi(1,1,4+4*(ith-1)),
     &              realg3,realg4,my_node_x,my_comm_x,my_node_world)
            end if
c            call mpi_barrier(mpi_comm_world,ierror)
            call wall_time(t2)
            tcom = tcom + (t2-t1)
            t1=t2

c
c     Forward Fourier transform of nonlinear term H_i in z
c
c            call mpi_barrier(mpi_comm_world,ierror)
            call wall_time(t2)
            tser = tser + (t2-t1)
            t1=t2

            call putpxz_z(boxr(1,1,4+4*(ith-1)),boxi(1,1,4+4*(ith-1)),yb
     &           ,8+pressure+3*(ith-1),1,ur,ui,
     &           realg1,realg2,my_node_z,my_comm_z,my_node_world)

c            call mpi_barrier(mpi_comm_world,ierror)
            call wall_time(t2)
            tcom = tcom + (t2-t1)
            t1=t2

         end do
#endif
      end if
c
c     Now the nonlinear term is in Fourier/Fourier/Real space.
c     Note that everything is multiplyed by nxp*nzp.
c

      end subroutine nonlinbl
