c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/filter_cutoff.f $
c $LastChangedDate: 2008-05-05 11:09:08 +0200 (Mon, 05 May 2008) $
c $LastChangedBy: qiang@MECH.KTH.SE $
c $LastChangedRevision: 1175 $
c
c ***********************************************************************
      subroutine filter_cutoff(gu2r,gu2i,ll,npl,nxy,prexn,prezn,wr,wi,
     &     yb,realg5,realg6,my_node_x,my_comm_x,my_node_world,
     &     gboxr_z,gboxi_z,wr_z,wi_z,wr_x,wi_x)
c
c     Spectral cutoff (lowpass) filter at w_c=pi/2
c     filters ll variables of a box in x/z direction
c     input and output is in physical space
c
      implicit none

      include 'par.f'

      integer i,x,y,z,ll,npl,xy,nxy,yb

      real gu2r((nxp/2+1)*mby,nzd_new/nprocz,ll)
      real gu2i((nxp/2+1)*mby,nzd_new/nprocz,ll)
      real gboxr_z(memnx,nzd_new,6)
      real gboxi_z(memnx,nzd_new,6)
      real prexn(nx+15),prezn(nz*2+15)
      real wr(nxp/2+1,mby,nzd_new/nprocz),wi(nxp/2+1,mby,nzd_new/nprocz)
      real wr_z(memnx,nzd_new)
      real wi_z(memnx,nzd_new)
      real wr_x(nxp/2+1,nzd_new/nprocz)
      real wi_x(nxp/2+1,nzd_new/nprocz)

      integer realg5,realg6,my_node_x,my_comm_x,my_node_world
      integer xs,xe
      integer zs,ze

      do i=1,ll
         if (nproc.eq.1) then
c
c     Forward transform
c
            call vrfftf(gu2r(1,1,i),gu2i(1,1,i),wr,wi,
     &           nx,nzpc*mby,1,nxp/2+1,prexn)
            call vcfftf(gu2r(1,1,i),gu2i(1,1,i),wr,wi,nz,
     &           nxy,(nxp/2+1)*mby,1,prezn)
c
c     Spectral filter
c
            xs = nx/4+1
            xe = nx/2+1
            do z=1,nz
               do y=1,npl
                  do x=xs,xe
                     xy=x+(y-1)*(nxp/2+1)
                     gu2r(xy,z,i) = 0.
                     gu2i(xy,z,i) = 0.
                  end do
               end do
            end do
            zs = nz/4+2
            ze = nz-nz/4 
            do z=zs,ze
               do y=1,npl
                  do x=1,nx/2+1
                     xy=x+(y-1)*(nxp/2+1)
                     gu2r(xy,z,i) = 0.
                     gu2i(xy,z,i) = 0.
                  end do
               end do
            end do
c
c     Backtransform and scaling
c
            call vcfftb(gu2r(1,1,i),gu2i(1,1,i),wr,wi,nz,
     &           nxy,(nxp/2+1)*mby,1,prezn)
            call vrfftb(gu2r(1,1,i),gu2i(1,1,i),wr,wi,
     &           nx,nzpc*mby,1,nxp/2+1,prexn)
            do z=1,nzc
               do y=1,npl
                  do x=1,nx/2
                     xy=x+(y-1)*(nxp/2+1)
                     gu2r(xy,z,i) = gu2r(xy,z,i) / real(nx*nz)
                     gu2i(xy,z,i) = gu2i(xy,z,i) / real(nx*nz)
                  end do
               end do
            end do
         else
#ifdef MPI
            call vrfftf(gu2r(1,1,i),gu2i(1,1,i),wr_x,wi_x,
     &           nx,memnz,1,nxp/2+1,prexn) 
            call putpxz_x(gu2r(1,1,i),gu2i(1,1,i),yb,i,0,
     &           gboxr_z(1,1,i),gboxi_z(1,1,i),
     &           realg5,realg6,my_node_x,my_comm_x,my_node_world)
            call vcfftf(gboxr_z(1,1,i),gboxi_z(1,1,i),wr_z,wi_z,
     &           nz,memnx,memnx,1,prezn)
c
c     Spectral filter
c
            zs = nz/4+2
            ze = nz-nz/4 
            do z=zs,ze
               do y=1,npl
                  do x=1,memnx
                     xy=x+(y-1)*(nxp/2+1)
                     gboxr_z(xy,z,i) = 0.
                     gboxi_z(xy,z,i) = 0.
                  end do
               end do
            end do
            call getpxz_x(gu2r(1,1,i),gu2i(1,1,i),yb,i,0,
     &           gboxr_z(1,1,i),gboxi_z(1,1,i),
     &           realg5,realg6,my_node_x,my_comm_x,my_node_world)
            xs = nx/4+1
            xe = nx/2+1
            do z=1,memnz
               do y=1,npl
                  do x=xs,xe
                     xy=x+(y-1)*(nxp/2+1)
                     gu2r(xy,z,i) = 0.
                     gu2i(xy,z,i) = 0.
                  end do
               end do
            end do
            call putpxz_x(gu2r(1,1,i),gu2i(1,1,i),yb,i,0,
     &           gboxr_z(1,1,i),gboxi_z(1,1,i),
     &           realg5,realg6,my_node_x,my_comm_x,my_node_world)
c
c     Backtransform and scaling
c
            call vcfftb(gboxr_z(1,1,i),gboxi_z(1,1,i),wr_z,wi_z,
     &           nz,memnx,memnx,1,prezn)
            call getpxz_x(gu2r(1,1,i),gu2i(1,1,i),yb,i,0,
     &           gboxr_z(1,1,i),gboxi_z(1,1,i),
     &           realg5,realg6,my_node_x,my_comm_x,my_node_world)
            call vrfftb(gu2r(1,1,i),gu2i(1,1,i),wr_x,wi_x,
     &           nx,memnz,1,nxp/2+1,prexn)

            do z=1,memnz
               do y=1,npl
                  do x=1,nx/2
                     xy=x+(y-1)*(nxp/2+1)
                     gu2r(xy,z,i) = gu2r(xy,z,i) / real(nx*nz)
                     gu2i(xy,z,i) = gu2i(xy,z,i) / real(nx*nz)
                  end do
               end do
            end do
#endif
         end if
      end do
      
      end subroutine filter_cutoff








