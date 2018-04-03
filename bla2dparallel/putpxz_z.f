      subroutine putpxz_z(boxr_z,boxi_z,yb,iii,ipad,ur,ui,
     &     realg1,realg2,my_node,my_comm,my_node_world)
c      
c     Get an xz box from ur,ui with MPI communication
c     Pad for dealiazing if ipad = 1
c
      implicit none

      include 'par.f'
      include 'mpif.h'

      integer yb,iii,ipad

      real boxr_z(memnx,nzd_new)
      real boxi_z(memnx,nzd_new)
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      integer realg1,realg2,my_node,my_comm,my_node_world

      integer x,z,nzpad,yyp,ierror
c
c     This routine shall only be called if more than one 
c     MPI processors are used, otherwise use getxz.f
c
      if (nproc.eq.1) call stopnow(3898)
c
c     Base index (yyp=yb for node 0)
c
      yyp = yb-my_node
c
c     Copy to the coarse grid
c
      do z=nz/2+1,nz
         do x=1,memnx
            boxr_z(x,z)=boxr_z(x,z+nzp-nz)
            boxi_z(x,z)=boxi_z(x,z+nzp-nz)
         end do
      end do
c
c     Do the global communication
c
      call mpi_alltoall(boxr_z,1,
     &     realg2,ur(1,yyp,1,iii),1,realg1,
     &     my_comm,ierror)
      call mpi_alltoall(boxi_z,1,
     &     realg2,ui(1,yyp,1,iii),1,realg1,
     &     my_comm,ierror)
      

      end subroutine putpxz_z
