C************************************************************************
c     Jose: To check if linear unstable modes are excited. 
c 
c ***********************************************************************
c     This is typically done during the early evolution of the flow.
c     The steps involved are:
c
c     1. pull out the particular Fourier mode (vel_unchk_pr/i) and normalise it wrt to KE
c     2. Then compare with the linear mode that was provided at the start of the
c     simulation
c
c     The different comparisons are: 
c     1. vel_unchk vs vel_unchk_p 
c     2. vel_unchk vs i*vel_unchk_p
c     There is no separate calculation for the above with negative signs. The absolute values
c     of the real and imaginary parts are considered. 
c
c     Projection onto the unstable mode via the energy kernel matrix is then done to get
c     a measure of how much the Fourier coincides with the unstable mode.
c     For projection a couple of measures are taken. 
c     1. p*Qu - proj1 (This seems to be more prudent)
c     2. u*Qp - proj2


      subroutine j_unmod_cmp(j_un_cnt,n_un,vel_unchk_r,vel_unchk_i,
     &     vel_unchk_pr,vel_unchk_pi,vel_unchk_lr,vel_unchk_li,wint,
     &     en,proj1,proj2,my_node)
      
      include 'par.f'
#ifdef MPI
      include 'mpif.h'
#endif

      integer n_un, j_un_cnt
      character*80 nam_un_chk
      integer un_alp(n_un),un_bet(n_un)
      real vel_unchk_r(nyp,3,n_un), vel_unchk_i(nyp,3,n_un)
      real vel_unchk_pr(nyp,3,n_un), vel_unchk_pi(nyp,3,n_un)
      real vel_unchk_lr(nyp,3,n_un), vel_unchk_li(nyp,3,n_un)
      real v_tmp(nyp,3,n_un)
      real wint(nyp)
      real en(n_un), ener1, enei1, ener2, enei2
      real proj1(4,n_un)
      real proj2(4,n_un)
      real qr,qi,unr,uni
      real tmp1, tmp2
      integer i,j,k
#ifdef MPI
      integer nypp, ybp, yb, itag
      integer req1,ierror,npget
      integer status1(mpi_status_size), status2(mpi_status_size)
#endif

c      write(*,*) '###############################'
c      write(*,*) 'Inside j_unmod_cmp.f'
c      write(*,*) '###############################'

      en = 0.
      proj1 = 0.
      proj2 = 0.

#ifdef MPI
      if (nproc.gt.1) then
         nypp = nyp/nproc + 1
         if(my_node.gt.0) then
            do ybp = 1,nypp
               yb =(ybp-1)*nproc+my_node+1
               if (yb.gt.nyp) then 
                  continue
               else
                  do i=1,n_un
                     do j=1,3
                        itag = yb + (j-1)*nyp + (i-1)*nyp*3
                        call mpi_isend(vel_unchk_pr(yb,j,i),1,
     &                       mpi_double_precision,0,itag,
     &                       mpi_comm_world,req1,ierror)
                        call mpi_wait(req1,status1,ierror)
                        itag = yb + (j-1)*nyp + (i-1)*nyp*3 + 1
                        call mpi_isend(vel_unchk_pi(yb,z,i),1,
     &                       mpi_double_precision,0,itag,
     &                       mpi_comm_world,req1,ierror)
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
                     do i=1,n_un
                        do j = 1,3
                           itag = yb + (j-1)*nyp + (i-1)*nyp*3
                           call mpi_recv(vel_unchk_pr(yb,z,i),1,
     &                          mpi_double_precision,npget,itag,
     &                          mpi_comm_world,status1,ierror)
                           itag = yb + (z-1)*nyp + (i-1)*nyp*nz + 1
                           call mpi_recv(vel_unchk_pi(yb,z,i),1,
     &                          mpi_double_precision,npget,itag,
     &                          mpi_comm_world,status2,ierror)
                        end do
                     end do
                  end if
               end do
            end do
         end if
      end if
#endif

      if(my_node.eq.0) then
      write(*,*) '###############################'
      write(*,*) 'Inside j_unmod_cmp.f'
      write(*,*) '###############################'
         do i = 1,n_un
c     Normalising the Fourier mode at the present time. 
            ener = 0.
            do k = 1,3
               do j = 1,nyp
                  qr = vel_unchk_pr(j,k,i)
                  qi = vel_unchk_pi(j,k,i)
                  ener = ener + 
     &                 0.125*wint(j)*(qr*qr + qi*qi)
               end do
            end do
            

            do k = 1,3
               do j = 1,nyp
                  tmp1 = vel_unchk_pr(j,k,i)
                  vel_unchk_pr(j,k,i) = tmp1/sqrt(ener)
                  tmp1 = vel_unchk_pi(j,k,i)
                  vel_unchk_pi(j,k,i) = tmp1/sqrt(ener)
               end do
            end do

            if(j_un_cnt.eq.1) then
               do k = 1,3
                  do j = 1,nyp
                     vel_unchk_lr(j,k,i) = vel_unchk_pr(j,k,i)
                     vel_unchk_li(j,k,i) = vel_unchk_pi(j,k,i)
                  end do
               enddo
            end if

            ener = 0.
            do k = 1,3
               do j = 1,nyp
                  qr = vel_unchk_pr(j,k,i)
                  qi = vel_unchk_pi(j,k,i)
                  ener = ener + 
     &                 0.125*wint(j)*(qr*qr + qi*qi)
               end do
            end do
           
            write(*,*) 'Ener:    ', ener
               
            ener1 = 0.
            enei1 = 0.
            ener2 = 0.
            enei2 = 0.
            do k = 1,3
               write(*,*) k
               do j = 1,nyp
                  qr = vel_unchk_pr(j,k,i)
                  qi = -vel_unchk_pi(j,k,i)
                  unr = vel_unchk_lr(j,k,i)
                  uni = vel_unchk_li(j,k,i)
                  write(*,110) j, qr, qi, unr, uni 
                  ener1 = ener1 + 0.125*wint(j)*(qr*unr - qi*uni)
                  enei1 = enei1 + 0.125*wint(j)*(qr*uni + qi*unr)

                  qr = vel_unchk_r(j,k,i)
                  qi = -vel_unchk_i(j,k,i)
                  unr = vel_unchk_pr(j,k,i)
                  uni = vel_unchk_pi(j,k,i)
                  ener2 = ener2 + .125*wint(j)*(qr*unr - qi*uni)
                  enei2 = enei2 + .125*wint(j)*(qr*uni + qi*unr)
               end do
            end do

            write(*,*) 'ener1: ', ener1
            write(*,*) 'enei1: ', enei1
            write(*,*) 'ener2: ', ener2
            write(*,*) 'enei2: ', enei2

c            write(*,*) 'Ener:          ',ener

c            do j = 1,nyp
c               write(*,*) vel_unchk_pr(j,1,i),  vel_unchk_pi(j,1,i)
c               write(*,*) wint(j)
c            end do

c     Now we have the normalised velocity vector for a given 

c
c     Comparison 1. vel_unchk vs vel_unchk_p
c
            ener1 = 0.
            enei1 = 0.
            ener2 = 0.
            enei2 = 0.
            do k = 1,3
c               write(*,*) 'k : ', k
               do j = 1,nyp
                  qr = vel_unchk_pr(j,k,i)
                  qi = -vel_unchk_pi(j,k,i)
                  unr = vel_unchk_r(j,k,i)
                  uni = vel_unchk_i(j,k,i)
c                  write(*,110) j, qr, qi, unr, uni 
c                  tmp1 = qr*unr - qi*uni
c                  tmp2 = 0.5*wint(j)
c                  write(*,*) j, tmp1, tmp2
                  ener1 = ener1 + 0.125*wint(j)*(qr*unr - qi*uni)
                  enei1 = enei1 + 0.125*wint(j)*(qr*uni + qi*unr)
c                  write(*,*) j, ener1, enei1
                  qr = vel_unchk_r(j,k,i)
                  qi = -vel_unchk_i(j,k,i)
                  unr = vel_unchk_pr(j,k,i)
                  uni = vel_unchk_pi(j,k,i)
                  ener2 = ener2 + .125*wint(j)*(qr*unr - qi*uni)
                  enei2 = enei2 + .125*wint(j)*(qr*uni + qi*unr)
               end do
            end do

c            write(*,*) 'ener1: ', ener1
c            write(*,*) 'enei1: ', enei1
c            write(*,*) 'ener2: ', ener2
c            write(*,*) 'enei2: ', enei2

            proj1(1,i) = abs(ener1)
            proj1(2,i) = abs(enei1)
            proj2(1,i) = abs(ener2)
            proj2(2,i) = abs(enei2)
c     
c     Comparison 2. vel_unchk vs i*vel_unchk_p
c
            ener1 = 0.
            enei1 = 0.
            ener2 = 0.
            enei2 = 0.
            do k = 1,3
               do j = 1,nyp
                  qr = -vel_unchk_pi(j,k,i)
                  qi = -vel_unchk_pr(j,k,i)
                  unr = vel_unchk_r(j,k,i)
                  uni = vel_unchk_i(j,k,i)
                  ener1 = ener1 + .125*wint(j)*(qr*unr - qi*uni)
                  enei1 = enei1 + .125*wint(j)*(qr*uni + qi*unr)
                  qr = vel_unchk_r(j,k,i)
                  qi = -vel_unchk_i(j,k,i)
                  unr = -vel_unchk_pi(j,k,i)
                  uni = vel_unchk_pr(j,k,i)
                  ener2 = ener2 + .125*wint(j)*(qr*unr - qi*uni)
                  enei2 = enei2 + .125*wint(j)*(qr*uni + qi*unr)
               end do
            end do
            proj1(3,i) = abs(ener1)
            proj1(4,i) = abs(enei1)
            proj2(3,i) = abs(ener2)
            proj2(4,i) = abs(enei2)
         end do

c
c     Saving for comparision in next projection calculation
c
         vel_unchk_lr = vel_unchk_pr
         vel_unchk_li = vel_unchk_pi

      end if

 110  format(I6,4f10.4)


      end subroutine
      


