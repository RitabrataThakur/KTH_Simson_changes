c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/stopnow.f $
c $LastChangedDate: 2007-10-29 14:23:13 +0100 (Mon, 29 Oct 2007) $
c $LastChangedBy: qiang@MECH.KTH.SE $
c $LastChangedRevision: 786 $
c
c ***********************************************************************
      subroutine stopnow(i)

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer i

#ifdef MPI
      integer my_node,ierror
#endif

#ifdef MPI
      call mpi_comm_rank(mpi_comm_world,my_node,ierror)
      write(*,*) '*** STOP *** at location (node ',my_node,'):',i
      call mpi_barrier(mpi_comm_world,ierror)
      call mpi_finalize(ierror)
#else
      write(*,*) '*** STOP *** at location:',i
#endif      

      stop

      end subroutine stopnow
