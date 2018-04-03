C************************************************************************
c     Jose: To check if linear unstable modes are excited. 
c     Parameters of the unstable mode, and the unstable mode are read here
c ***********************************************************************

      subroutine r_un_chk(j_unch,j_un_chk,n_un,nam_un_chk,un_alp,un_bet,
     &     vel_unchk_r,vel_unchk_i,wint)

      include 'par.f'

      logical j_un_chk, stp
      integer n_un,j_unch
      character*80 nam_un_chk, nam_tmp
      integer un_alp(j_unch),un_bet(j_unch)
      real vel_unchk_r(nyp,3,j_unch), vel_unchk_i(nyp,3,j_unch)
      real wint(nyp)
      real qr,qi, ener
      integer i,j,k

      open(unit=10,status='old',file='unstab_chk.i')

c     Flag for unstable mode check
      read(10,*) j_un_chk
      if(j_un_chk) then
c     Number of unstable modes to be checked.
         read(10,*) n_un
         if(n_un.gt.j_unch) then
            write(ios,*) 'Reduce n_un or change j_unch in bla.f'
            write(ios,*) 'Program stop'
            stp = .true.
         end if
         
         if(n_un.eq.0) stp =.true.

         if(stp) call stopnow(8913)

c     Read the wavenumbers and the unstable mode
        
         ener = 0.0
         do i=1,n_un
            read(10,*) nam_tmp
            read(10,*) un_alp(i)
            read(10,*) un_bet(i)
            
            open(unit=13,status='old',file=nam_tmp)
            do j = 1,nyp
               read (13,*) vel_unchk_r(j,1,i),vel_unchk_i(j,1,i),
     &              vel_unchk_r(j,2,i),vel_unchk_i(j,2,i),
     &              vel_unchk_r(j,3,i),vel_unchk_i(j,3,i)
               qr = vel_unchk_r(j,1,i)
               qi = vel_unchk_i(j,1,i)
               ener = ener + 0.125*wint(j)*(qr*qr + qi*qi)
               qr = vel_unchk_r(j,2,i)
               qi = vel_unchk_i(j,2,i)
               ener = ener + 0.125*wint(j)*(qr*qr + qi*qi)
               qr = vel_unchk_r(j,3,i)
               qi = vel_unchk_i(j,3,i)
               ener = ener + 0.125*wint(j)*(qr*qr + qi*qi)
            end do

c
c     Jose: I am leaving this here. This can be a reference later on
c     for checking what is the norm that I have used. 
            write(*,*) 'Mode : ', i, 'Energy = ', ener
           
            close(unit=13)
         end do

c     Output file where check results are to be saved
         read(10,*) nam_un_chk         
         
      end if

      close(unit=10)

      end subroutine
