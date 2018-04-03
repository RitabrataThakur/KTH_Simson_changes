      subroutine getpxz_z(boxr_z,boxi_z,yb,iii,ipad,ur,ui,
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
c     Do the global communication
c
      call mpi_alltoall(ur(1,yyp,1,iii),1,
     &     realg1,boxr_z,1,realg2,
     &     my_comm,ierror)
      call mpi_alltoall(ui(1,yyp,1,iii),1,
     &     realg1,boxi_z,1,realg2,
     &     my_comm,ierror)

c
c     Copy to the upper part of boxr_z,boxi_z
c
      if (nfzsym.eq.0) then
         nzpad=(nzp-nz)*ipad
         if (nzpad.ne.0) then
            do z=nz/2+1,nz
               do x=1,memnx
                  boxr_z(x,z+nzpad)=boxr_z(x,z)
                  boxi_z(x,z+nzpad)=boxi_z(x,z)
               end do
            end do
         end if
      end if
c     
c     Pad with zeros
c
      if (ipad.eq.1) then
         do z=(nz+1)/2+1,nzp+1-nz/2
            do x=1,memnx
               boxr_z(x,z)=0.0
               boxi_z(x,z)=0.0
            end do
         end do
c     
c     Oddball zeroing only
c     
      else
         if (mod(nz,2).eq.0) then
            do x=1,memnx
               boxr_z(x,nz/2+1)=0.0
               boxi_z(x,nz/2+1)=0.0
            end do
         end if
      end if
      

      end subroutine getpxz_z
