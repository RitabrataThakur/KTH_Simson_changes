c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/bla/stopnow.f $
c $LastChangedDate: 2007-11-12 13:19:01 +0100 (Mon, 12 Nov 2007) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 850 $
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
