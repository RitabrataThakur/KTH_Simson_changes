C************************************************************************
c     Jose: For reading some parameters for evaluating entropy measure
c 
c ***********************************************************************

      subroutine r_jent(j_en_flag,en_ty_fl,nam_entj,i_jen,t_je)
      
      include 'par.f'

      logical j_en_flag
      real t_je, zer
      integer i_jen,en_ty_fl
      character*80 nam_entj
      
      zer = 0.d0

      open(unit=27,status='old',file='entropy.i')

    
C     Check if entropy is to be measured.
      read(27,*) j_en_flag
      
      if(j_en_flag) then
         read(27,*) en_ty_fl
c     Read the type of calculation to be made
c     en_ty_fl : 1 = Entropy calculation for 6 X planes.
c     en_ty_fl : 2 = Entropy calculation averaged over the box
         read(27,*) nam_entj
C     Read the interval for evaluating the entropy
         read(27,*) i_jen

C     Read the start time of entropy calculation
         read(27,*) t_je

C     Read the number of points at which velocities are to be saved
C     These points all correspond to the (0,0) in the x-z plane
C         read(27,*) n_vel

C     Evaluate the y indices where the velocities are to be saved.
c         if((n_vel.gt.0).and.(n_vel.le.10).and.(my_node.eq.0)) then
c            open(unit=28,status='new',file='traj_coord.dat')

c            indstep = (nyp-1)/(n_vel+1)
c            do i=1,n_vel,1
c               en_yind(i) = 1 + i*indstep
cc     Write the coordinates of the chosen point for tracking velocities
c               write(28,*) en_yind(i), zer, eta(en_yind(i)), zer
c            end do
c            close(28)
c         end if
      end if

      close(27)
      
      return
      
      end subroutine

      
