c ***********************************************************************
c
c $HeadURL: https://www2.mech.kth.se/svn/simson/trunk/bla/pre_nlin_3d.f $
c $LastChangedDate: 2008-08-26 09:14:54 +0200 (Tue, 26 Aug 2008) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 1283 $
c
c ***********************************************************************
      subroutine pre_nlin_3d(ur,ui,my_node_x,my_comm_x,my_node_world,
     & my_node_z,my_comm_z,
     & prex,pres,prez,prea,
     & realg1,realg2,realg3,realg4,
     & u2bf3_r,u2bf3_i,yb,ybp,
     & boxi_z,boxr_z,wr_x,wr_z)
c
c     Transform to physical space (vel and vort) and return on the
c     dealiased grid.
c
      implicit none

      include 'par.f'

      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real prex(nxp+15)
      real prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)
      real u2bf3_r((nxp/2+1)*mby,nzd_new/nprocz,nyp/nprocz+1,6)
      real u2bf3_i((nxp/2+1)*mby,nzd_new/nprocz,nyp/nprocz+1,6)
      integer realg1,realg2,realg3,realg4
      integer my_node_x,my_comm_x,my_node_world
      integer my_node_z,my_comm_z
      integer yb,ybp,i
      integer x,z
 
#ifdef MPI
      integer ierror
      real boxr_z(memnx,nzd_new,6)
      real boxi_z(memnx,nzd_new,6)
c      real boxr(memnx,nzd_new,4*scalar)
c      real boxi(memnx,nzd_new,4*scalar)
      real wr_z(memnx,nzd_new)
      real wi_z(memnx,nzd_new)
      real wr_x(nxp/2+1,nzd_new/nprocz)
      real wi_x(nxp/2+1,nzd_new/nprocz)
      real aaar(nxp/2+1,nzd_new,2)
      real aaai(nxp/2+1,nzd_new,2)
#endif


#ifdef MPI


         do i=1,6
            call getpxz_z(boxr_z(1,1,i),boxi_z(1,1,i),yb,i,1,ur,ui,
     &           realg1,realg2,my_node_z,my_comm_z,my_node_world)
         end do

c     
c     Backward Fourier transform in z 
c
         do i=1,6    
            call vcfftb(boxr_z(1,1,i),boxi_z(1,1,i),wr_z,wi_z,
     &          nzp,memnx,memnx,1,prez)
         end do

c
c     Get velocities/vorticities of base flow into u2r(:,:,1-6)
c

         do i=1,6
            call getpxz_x(u2bf3_r(1,1,ybp,i),u2bf3_i(1,1,ybp,i),yb,i,1,
     &          boxr_z(1,1,i),boxi_z(1,1,i),
     &          realg3,realg4,my_node_x,my_comm_x,my_node_world)
         end do

c
c     Backward Fourier transform in x
c

c     BEGIN DEBUG, CLEAN UP LATER

c         if (my_node_world.eq.0.and.ybp.eq.1) then
c        
c            write(*,*) 'before fft'
c            write(*,*) 'u real',u2bf3_r(1,1,1,1)
c            write(*,*) 'u imag',u2bf3_i(1,1,1,1)
c            write(*,*) 'u real2',u2bf3_r(2,1,1,1)
c            write(*,*) 'u imag2',u2bf3_i(2,1,1,1)
         
c         end if   

c     END DEBUG


         do i=1,6
            call vrfftb(u2bf3_r(1,1,ybp,i),u2bf3_i(1,1,ybp,i),wr_x,wi_x,
     &          nxp,nzd_new/nprocz,1,nxp/2+1,prex)
         end do

c     BEGIN DEBUG, CLEAN UP LATER


c         if (my_node_world.eq.0.and.ybp.eq.1) then
c            write(*,*) 'after fft'
c            write(*,*) 'u real',u2bf3_r(1,1,1,1)
c            write(*,*) 'u imag',u2bf3_i(1,1,1,1)
c            write(*,*) 'u real2', u2bf3_r(2,1,1,1)
c            write(*,*) 'u imag2', u2bf3_i(2,1,1,1)
         
c         end if
         
c         if (my_node_world.eq.0) then

c            write(*,*), 'head node,ybp,ureal',ybp,u2bf3_r(1,1,1,1)
c         end if

c     END DEBUG

#endif

      end subroutine pre_nlin_3d
