c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/wdiscp.f $
c $LastChangedDate: 2007-11-02 09:31:18 +0100 (Fri, 02 Nov 2007) $
c $LastChangedBy: qiang@MECH.KTH.SE $
c $LastChangedRevision: 789 $
c
c ***********************************************************************
      subroutine wdiscp(ur,ui,re,pr,m1,xl,zl,t,xs,dstar,fltype,
     &     bstart,bslope,rlam,spanv,namnut,gall,
     &     boxr,boxi,urx,alfa,zs,beta,my_node_world)
c
c     Writes m variables from ur,ui to file namnut
c
      implicit none

      include 'par.f'
#ifdef MPI
      include 'mpif.h'
#endif

      character*80 namnut
      integer fltype
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real re,pr(scalar),xl,zl,t,xs,dstar
      real bstart,bslope,rlam,spanv,m1(scalar)
      logical gall
      real urx(nx)
      real boxr(memnx,mbz,nyp),boxi(memnx,mbz,nyp)
      real alfa(nx/2*mbz),zs,beta(nz)
      real urtemp(memnx),uitemp(memnx)
      integer x,y,z,i,zb,ii,iii,zb_t
      real uw0low
c
c     MPI
c
      integer my_node_world,zbp,ip
#ifdef MPI
      integer ierror
      integer status1(mpi_status_size),status2(mpi_status_size)

      if (nproc.gt.1) call mpi_barrier(mpi_comm_world,ierror)
#endif      

      if (pressure.eq.0) then
         if (my_node_world.eq.0) then
            write(*,*) 'Pressure is not included in the simulation'
         end if
         call stopnow(453634)
      end if


      if (my_node_world.eq.0) then
c
c     Write file header
c     Data is shifted so that xs,zs=0 for all fields
c     Any existing file will be overwritten
c
         open(unit=11,file=trim(namnut)//'.p',form='unformatted')
         rewind(11)
         if (scalar.ge.1) then
            write(11) re,.false.,xl,zl,t,0.,(pr(i),m1(i),i=1,scalar)
         else
            write(11) re,.false.,xl,zl,t,0.
         end if
         write(11) nx,nyp,nzc,nfzsym
         write(11) fltype,dstar
         if (fltype.lt.0) write(11) rlam
         if (fltype.ge.4) write(11) bstart,bslope,rlam,spanv
      end if

      do i=1,1
c
c     Choose which field to write
c
         ii = 8

        
         do ip=0,nprocz-1
            do zb_t=1,memnz
               do y=1,nyp
                  do iii=0,nprocx-1
#ifdef MPI
                     if (nproc.gt.1.and.my_node_world.eq.
     &                    ip*nprocx+iii.and.(ip+iii).ne.0) then
c     
c     Send complete field to processor 0
c
                        call mpi_send(ur(1,y,zb_t,ii),memnx,
     &                       mpi_double_precision,0,ip+iii+100,
     &                       mpi_comm_world,ierror)
                        call mpi_send(ui(1,y,zb_t,ii),memnx,
     &                       mpi_double_precision,0,ip+iii+200,
     &                       mpi_comm_world,ierror)
                     else
                        do x=1,memnx
                           urtemp(x)=ur(x,y,zb_t,ii)
                           uitemp(x)=ui(x,y,zb_t,ii)
                        end do
                     end if
c     
c     Receive individual fields from processors >0
c
                     if (my_node_world.eq.0.and.nproc.gt.1
     &                    .and.(ip+iii).ne.0) then
                        call mpi_recv(urtemp,memnx,
     &                       mpi_double_precision,
     &                       ip*nprocx+iii,ip+iii+100,
     &                       mpi_comm_world,status1,ierror)
                        call mpi_recv(uitemp,memnx,
     &                       mpi_double_precision,
     &                       ip*nprocx+iii,ip+iii+200,
     &                       mpi_comm_world,status2,ierror)
                     end if
c
c     Wait until communication is finished
c
                     if (nproc.gt.1) then 
                        call mpi_barrier(mpi_comm_world,ierror)
                     end if
#endif
                     if (nproc.eq.1) then
                        do x=1,memnx
                           urtemp(x)=ur(x,y,zb_t,ii)
                           uitemp(x)=ui(x,y,zb_t,ii)
                        end do
                     end if
c
c     Shift data so that xs,zs=0 for all fields
c
                     zb=zb_t+ip*memnz
                     if (my_node_world.eq.0) then
                        call xysh(urtemp,uitemp,xs,zs,
     &                       alfa,beta,zb,iii)
                        do x=1,memnx
                           urx(2*iii*memnx+2*x-1)= urtemp(x)
                           urx(2*iii*memnx+2*x)  = uitemp(x)
                        end do
                     end if
                  end do
                 
                  if (my_node_world.eq.0) then    
                     write(11) urx
                  end if
               end do
            end do
         end do
      end do


     
      
      if (my_node_world.eq.0) then
c
c     Close the file and flush...
c
         call cflush(11)
         close(unit=11)
      end if

      end subroutine wdiscp
