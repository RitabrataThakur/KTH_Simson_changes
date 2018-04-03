      subroutine getpxz_x(boxr_x,boxi_x,yb,iii,ipad,boxr_z,boxi_z,
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

      integer x,z,nzpad,ierror
c
c     This routine shall only be called if more than one 
c     MPI processors are used, otherwise use getxz.f
c
      if (nproc.eq.1) call stopnow(3898)
c
c     Do the global communication
c
      call mpi_alltoall(boxr_z,1,
     &     realg3,boxr_x,1,realg4,
     &     my_comm,ierror)
      call mpi_alltoall(boxi_z,1,
     &     realg3,boxi_x,1,realg4,
     &     my_comm,ierror)

c
c     Pad with zeros
c
      if (ipad.eq.1) then
         do z=1,nzd_new/nprocz
            do x=nx/2+1,nxp/2+1
               boxr_x(x,z)= 0.0
               boxi_x(x,z)= 0.0
            end do
         end do
      end if
      
 
      end subroutine getpxz_x
