C
C     Jose: For a measure of "entropy"
C     S = sum(sum(f(y,z,x0))) for a given streamwise location 
      
c
     
      subroutine entropy_j(om_entj,om_prevr,om_previ,ent_cal,ent,
     &     en_ty_fl,eta,dz,my_node,wint)

      implicit none

      include 'par.f'
#ifdef MPI
      include 'mpif.h'
#endif
      
      real om_entj(nyp,nzd,(nxp/2+1))
      real om_prevr(nyp,nzd,(nxp/2+1)),om_previ(nyp,nzd,(nxp/2+1))
      real ent(nyp), dz, dy, diff2, sum
      integer ent_cal, en_ty_fl,my_node, i, j, k
      real eta(nyp)
      real wint(nyp)
#ifdef MPI
      integer nypp, ybp, yb, z, itag
      integer req1,ierror,npget
      integer status1(mpi_status_size), status2(mpi_status_size)
#endif
C
c     First we have to communicate all the values of om_entj to one processor.
c     Come to this a little later on.
C     The communication for the case when en_ty_fl = 1 is not done properly yet. 

      if(en_ty_fl.eq.1.and.ent_cal.eq.0) then
         
c     First gather the previous vorticity vector into 1 place. 
c     This will be done only once as there is only one reference state.
#ifdef MPI
         if (nproc.gt.1) then
            nypp = nyp/nproc+1
            if(my_node.gt.0) then
               do ybp=1,nypp
                  yb =(ybp-1)*nproc+my_node+1
                  if (yb.gt.nyp) then 
                     continue
                  else
                     do z=1,nzd
                        do i = 1,nx/2
                           itag = yb + (z-1)*nyp + (i-1)*nyp*nz
                           call mpi_isend(om_prevr(yb,z,i),1,
     &                          mpi_double_precision,0,itag,
     &                          mpi_comm_world,req1,ierror)
                           call mpi_wait(req1,status1,ierror)
                           itag = yb + (z-1)*nyp + (i-1)*nyp*nz + 1
                           call mpi_isend(om_previ(yb,z,i),1,
     &                          mpi_double_precision,0,itag,
     &                          mpi_comm_world,req1,ierror)
                           call mpi_wait(req1,status2,ierror)
                        end do
                     end do
                  end if
               end do
               
            else
               do k = 1,nproc-1
                  npget=my_node+k
                  do ybp=1,nypp
                     yb =(ybp-1)*nproc+npget+1
                     if (yb.gt.nyp) then 
                        continue
                     else
                        do z=1,nzd
                           do i = 1,nx/2
                              itag = yb + (z-1)*nyp + (i-1)*nyp*nz + 2
                              call mpi_recv(om_prevr(yb,z,i),1,
     &                             mpi_double_precision,npget,itag,
     &                             mpi_comm_world,status1,ierror)
                              itag = yb + (z-1)*nyp + (i-1)*nyp*nz + 3
                              call mpi_recv(om_previ(yb,z,i),1,
     &                             mpi_double_precision,npget,itag,
     &                             mpi_comm_world,status2,ierror)
                           end do
                        end do
                     end if
                  end do
               end do
            end if
         end if
#endif
      end if

      if(en_ty_fl.eq.1) then
#ifdef MPI
         if (nproc.gt.1) then
            nypp = nyp/nproc+1
            if(my_node.gt.0) then
               do ybp=1,nypp
                  yb =(ybp-1)*nproc+my_node+1
                  if (yb.gt.nyp) then 
                     continue
                  else
                     do z=1,nzd
                        do i = 1,6
                           itag = yb + (z-1)*nyp + (i-1)*nyp*nz + 4
                           call mpi_isend(om_entj(yb,z,i),1,
     &                          mpi_double_precision,0,itag,
     &                          mpi_comm_world,req1,ierror)
                           call mpi_wait(req1,status2,ierror)
                        end do
                     end do
                  end if
               end do
               
            else
               do k = 1,nproc-1
                  npget=my_node+k
                  do ybp=1,nypp
                     yb =(ybp-1)*nproc+npget+1
                     if (yb.gt.nyp) then 
                        continue
                     else
                        do z=1,nzd
                           do i = 1,6
                              itag = yb + (z-1)*nyp + (i-1)*nyp*nz + 4
                              call mpi_recv(om_entj(yb,z,i),1,
     &                             mpi_double_precision,npget,itag,
     &                             mpi_comm_world,status1,ierror)
                           end do
                        end do
                     end if
                  end do
               end do
            end if
         end if
#endif
      end if
C     
c     After all the values are brought to one single processor.
c     Check what happens when the loop is done from j = 1,nzd.
C     Don't conclude anything too early. 

      if(en_ty_fl.eq.1) then
         ent = 0.
         if (my_node.eq.0) then
c            if(ent_cal.le.1) then
c               om_prev = om_entj
c            else
            do i = 1,6
               sum = 0.
               do j = 1,nzd
                  do k = 1,nyp
                     diff2 = (om_entj(k,j,i) - om_prevr(k,j,i))**2
                     if((k.eq.1).or.(k.eq.nyp)) then
                        dy = abs(1.0 - eta(2))
                     else
                        dy = 0.5*abs(eta(k-1) - eta(k+1))
                     end if
                     sum = sum + diff2*dy*dz
                  end do
               end do
               ent(i) = sum
            end do
c     om_prev = om_entj
         end if
C         end if
      elseif(en_ty_fl.eq.2) then
#ifdef MPI
         if (nproc.gt.1) then
            nypp = nyp/nproc+1
            if(my_node.gt.0) then
               do ybp=1,nypp
                  yb =(ybp-1)*nproc+my_node+1
                  if (yb.gt.nyp) then 
                     continue
                  else
                     itag = yb + nypp
                     call mpi_isend(ent(yb),1,
     &                    mpi_double_precision,0,itag,
     &                    mpi_comm_world,req1,ierror)
                     call mpi_wait(req1,status2,ierror)                  
                  end if
               end do               
            else
               do k = 1,nproc-1
                  npget=my_node+k
                  do ybp=1,nypp
                     yb =(ybp-1)*nproc+npget+1
                     if (yb.gt.nyp) then 
                        continue
                     else
                        itag = yb + nypp
                        call mpi_recv(ent(yb),1,
     &                       mpi_double_precision,npget,itag,
     &                       mpi_comm_world,status1,ierror)
                     end if
                  end do
               end do
            end if
         end if
#endif
         sum = 0.d0 
         if(my_node.eq.0) then
c            sum = 0.d0 
            do k = 1,nyp
               sum = sum + ent(k)*wint(k)
            end do
         end if
         ent = 0.
         ent(1) = sum/2.d0
      end if

c#ifdef MPI
c      if (nproc.gt.1) then
c         do i = 1,6
c            call mpi_bcast(entyz(i),1,mpi_double_precision,0,
c     &           mpi_comm_world,ierror)
c         end do
c      end if
c#endif
      
      return
      
      end subroutine
