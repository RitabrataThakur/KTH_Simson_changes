c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/bla/bfvort3.f $
c $LastChangedDate: 2007-12-13 22:41:14 +0100 (Thu, 13 Dec 2007) $
c $LastChangedBy: pschlatt@MECH.KTH.SE $
c $LastChangedRevision: 1069 $
c
c ***********************************************************************
      subroutine bfvort3(ur,ui,alfa,beta,my_node_x,my_node_z,prey)
c
c     Compute vorticities for 3d base flow
c
      implicit none

      include 'par.f'
#ifdef MPI
      include 'mpif.h'
      integer ierror
#endif
      real ur(memnx,memny,memnz,6),ui(memnx,memny,memnz,6)
      real alfa(nx/2*mbz),beta(nz)
      real app1(nyp),app2(nyp)
      real prey(nyp*2+15),w(nxp/2+1,nyp)
      integer my_node_x,my_node_z,i,x,z,nn,zbp,y

      integer zp,xb,xp,nxp2,xy,y1,npp,yb

c
c     y-derivatives
c
      nn=0
      do i=1,3,2
         nn=nn+1
         do x=1,memnx
            do z=1,memnz
               do y=1,nyp
                  app1(y)=ur(x,y,z,i)
                  app2(y)=ui(x,y,z,i)
               end do
               call vchbf(app1,w,nyp,1,1,1,prey)
               call vchbf(app2,w,nyp,1,1,1,prey)
               call rdcheb(app1,nyp,1,1)
               call rdcheb(app2,nyp,1,1)
               call vchbb(app1,w,nyp,1,1,1,prey)
               call vchbb(app2,w,nyp,1,1,1,prey)
               do y=1,nyp
                  ur(x,y,z,7-i)=(-1)**nn*app1(y)*(2./real(nyp-1))
                  ui(x,y,z,7-i)=(-1)**nn*app2(y)*(2./real(nyp-1))


c     BEGIN DEBUG, CLEAN UP LATER

c                  if (i.eq.3.and.x.eq.1.and.y.eq.1) then
c                     if (my_node_x.eq.0) then

c                        write(*,*) 'znode,z,y=1,ur,ui',my_node_z,z,
c     &                       ur(x,y,z,7-i),ui(x,y,z,7-i)
c                     end if
c                  end if

c     END DEBUG

               end do
            end do
         end do
      end do
c
c     Vorticities
c


      do zbp=1,memnz
         z=my_node_z*memnz+zbp

c     BEGIN DEBUG, CLEAN UP LATER

c            write(*,*) 'nodex,nodez', my_node_x,my_node_z
c            write(*,*) 'z,zbp,beta(z)',z,zbp,beta(z)
  
c     END DEBUG
       
          do xp=1,memnx
             x=my_node_x*memnx+xp
            do y=1,nyp
c     omega_y
               ur(xp,y,zbp,5)=-beta(z)*ui(xp,y,zbp,1)+
     &              alfa(x)*ui(xp,y,zbp,3)
               ui(xp,y,zbp,5)= beta(z)*ur(xp,y,zbp,1)-
     &              alfa(x)*ur(xp,y,zbp,3)
c     omega_x
               ur(xp,y,zbp,4)=ur(xp,y,zbp,4)+beta(z)*ui(xp,y,zbp,2)
               ui(xp,y,zbp,4)=ui(xp,y,zbp,4)-beta(z)*ur(xp,y,zbp,2)
c     omega_z
               ur(xp,y,zbp,6)=ur(xp,y,zbp,6)-alfa(x)*ui(xp,y,zbp,2)
               ui(xp,y,zbp,6)=ui(xp,y,zbp,6)+alfa(x)*ur(xp,y,zbp,2)

c               if (y.eq.1) then

c                  write(*,*) 'y=1,ur,ui',ur(x,y,zbp,6),ui(x,y,zbp,6)
c               end if

            end do
         end do

         call mpi_barrier(mpi_comm_world,ierror)

      end do
      
      end subroutine bfvort3
