      subroutine w_mean_u(tav_u,my_node,eta,tav_cnt)

      implicit none

      include 'par.f'
#ifdef MPI
      include 'mpif.h'
#endif

      real tav_u(nyp)
      real temp
      integer my_node,k,nypp,tav_cnt
      real eta(nyp)
#ifdef MPI
      integer ybp, yb, z, itag
      integer req1,ierror,npget
      integer status1(mpi_status_size), status2(mpi_status_size)
#endif

#ifdef MPI

      if(nproc.gt.1) then
         nypp = nyp/nproc+1
         if(my_node.gt.0) then
            
            do ybp=1,nypp
               yb =(ybp-1)*nproc+my_node+1
c               write(*,*) 'Inside this loop'
c               write(*,*) my_node, yb
               if (yb.gt.nyp) then 
                  continue
               else
                  itag = yb
                  call mpi_isend(tav_u(yb),1,mpi_double_precision,0,
     &                 itag,mpi_comm_world,req1,ierror)
                  call mpi_wait(req1,status2,ierror)
               end if
            end do
         else
            
            do k = 1,nproc-1
               npget=my_node+k
               do ybp=1,nypp
                  yb =(ybp-1)*nproc+npget+1
c                  write(*,*) 'IN this loop yb: ', yb
                  if (yb.gt.nyp) then 
                     continue
                  else
c                     write(*,*) my_node,yb
                     itag = yb
                     call mpi_irecv(tav_u(yb),1,mpi_double_precision,
     &                    npget,itag,mpi_comm_world,req1,ierror)
                     call mpi_wait(req1,status2,ierror)
c                     write(*,*) my_node,yb,tav_u(yb)
                  end if
               end do
            end do
         end if
      end if
#endif         
      
      if(my_node.eq.0) then
         open(unit=89,file='U-mean_tav.dat')
c         write(*,*) '###############'
c         write(*,*) '###############'
c         write(*,*) '###############'
c         write(*,*) '###############'         
c         write(*,*) 'tav_cnt: ', tav_cnt
         do k = 1,nyp
            temp = tav_u(k)/real(tav_cnt)
            write(89,*) eta(k), temp
         end do
         close(unit=89)
c         tav_u = 0.*tav_u
      end if

c      write(*,*) 'In w_mean_u.f. ###############'
c      write(*,*) '###############'
c      write(*,*) '###############'
c      write(*,*) 'tav_cnt : ', tav_cnt
c      write(*,*) 'my_node : ', my_node
c      write(*,*) '###############'
c      write(*,*) '###############'

      return
      
      end subroutine
