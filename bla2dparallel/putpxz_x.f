      subroutine putpxz_x(boxr_x,boxi_x,yb,iii,ipad,boxr_z,boxi_z,
     &     realg3,realg4,my_node,my_comm,my_node_world)
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
      real boxr_x(nxp/2+1,nzd_new/nprocz)
      real boxi_x(nxp/2+1,nzd_new/nprocz)

      integer realg3,realg4,my_node,my_comm,my_node_world

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
c     Do the global communication
c
      call mpi_alltoall(boxr_x,1,
     &     realg4,boxr_z,1,realg3,
     &     my_comm,ierror)
      call mpi_alltoall(boxi_x,1,
     &     realg4,boxi_z,1,realg3,
     &     my_comm,ierror)


      end subroutine putpxz_x

