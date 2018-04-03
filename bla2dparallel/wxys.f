c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/wxys.f $
c $LastChangedDate: 2011-06-05 23:09:12 +0200 (Sun, 05 Jun 2011) $
c $LastChangedBy: pschlatt@MECH.KTH.SE $
c $LastChangedRevision: 1678 $
c
c ***********************************************************************
      subroutine wxys(xys,xysth,my_comm_z,my_comm_x,
     &     totcorr,corr,ncorr,corrf,
     &     totcorr_ms,corr_ms,
     &     totcorr_x,corr_x,ncorr_x,corrf_x,
     &     totcorr_ms_x,corr_ms_x,
     &     my_node_world,
     &     series,serf,nser,totseries,
     &     sumw,re,xl,zl,t,dstar,fltype,
     &     bstart,bslope,rlam,spanv,namxys,wxy,
     &     corrnam,corrxx,corryy,corrwsave,
     &     corrx,corry,corrxi,corryi,pr,m1,corrt,
     &     corrnam_x,corrzz,corryy_x,corrwsave_x,
     &     corrz,corry_x,corrzi,corryi_x,corrt_x,
     &     nsamp,namser,sercount,mhd_n,b0)
c
c     Sum the statistics xys and the correlations corr
c     from all the processors and writes it to a file
c
c     The global reduction could alternatively be implemented via
c     MPI_REDUCE operations.
c
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      include 'par.f'
      integer y,i,j,x,yb,z,xb
      real xys     (nx,nyp/nprocz+1,nxys)
      real xysth   (nx,nyp/nprocz+1,nxysth,scalar)
      real corr(nzc+2,mcorr)
      real totcorr(nzc+2)
      real totcorr_ms(4),corr_ms(4,mcorr)
      logical corrf
      real corr_x(nx+2,mcorr)
      real totcorr_x(nx+2)
      real totcorr_ms_x(4),corr_ms_x(4,mcorr)
      logical corrf_x
      integer ncorr,ncorr_x,ith

      real pr(scalar),m1(scalar)
      integer fltype
      real sumw,sumww
      real wxy(nx,nyp)
      real mhd_n,b0(3)
c
c     Two-point correlation
c
      real corrxx(mcorr),corryy(mcorr)
      real corrx(mcorr),corry(mcorr)
      integer corrxi(mcorr),corryi(mcorr)
      real cd(nzc+2)
      integer corrt(mcorr)
      real corrwsave(nzc+15),cw(nzc+2)
      real corrzz(mcorr),corryy_x(mcorr)
      real corrz(mcorr),corry_x(mcorr)
      integer corrzi(mcorr),corryi_x(mcorr)
      real cd_x(nx+2)
      integer corrt_x(mcorr)
      real corrwsave_x(nx+15),cw_x(nx+2)
c
c     Time series
c
      real totseries(msamp,0:mser)
      integer nser,nsamp,sercount
      logical serf
      character(len=80) namser

      real re,xl,zl,t,dstar,bstart,bslope,rlam,spanv

      character*80 namxys
      character*80 corrnam,corrnam_x
c
c     Time series
c
      real series(msamp,0:mser)
c
c     MPI
c
      integer my_node_world,my_comm_z,my_comm_x
#ifdef MPI
      real txys(nx),temp(nx)
      real tcorr(nzc+2),tcorr_ms(4)
      real tser(msamp,mser)
      integer ierror,ip,iii
      integer status1(mpi_status_size)
#endif
      xb=mod(my_node_world,nprocx)*memnx
c
c     Velocity, pressure and scalar xy-statistics
c     *******************************************
c
      call mpi_barrier(mpi_comm_world,ierror)

      if (my_node_world.eq.0) then
c
c     Open statistics file
c
         write(*,*) '** doing statistics'
         open(unit=19,file=namxys,form='unformatted')
         rewind(19)
         
c
c     Write header
c
         write(19) re,.false.,xl,zl,t,0.,(pr(i),m1(i),i=1,scalar)
c     Item 'A' added to identify if line is included or not when reading data
         write(19) 'A',mhd_n,b0
         write(19) nx,nyp,nzc,nfzsym
         write(19) fltype,dstar
         if (fltype.lt.0) write(19) rlam
         if (fltype.ge.6) write(19) bstart,bslope,rlam,spanv
c     
c     Write number of statistic quantities
c     
         write(19) sumw,nxys,nxysth,scalar

      end if
c
c     Loop over the velocity statistics
c
      do i=1,nxys
#ifdef MPI
         if (nproc.gt.1) then
            do y=1,nyp/nprocz+1

              
               call mpi_reduce(xys(1,y,i),txys,nx,
     &              mpi_double_precision,mpi_sum,0,
     &              my_comm_x,ierror)
               

               do ip=0,nprocz-1
                  yb=ip+1+(y-1)*nprocz
                  
                  if (yb.le.nyp) then
                     
                     if (ip.ge.1) then

                        if (my_node_world.eq.ip*nprocz) then
                           call mpi_ssend(txys,nx,
     &                          mpi_double_precision,0,
     &                          ip,
     &                          mpi_comm_world,ierror)
                           
                        end if
                        
                        
                        if (my_node_world.eq.0) then
                           call mpi_recv(txys,nx,
     &                          mpi_double_precision,
     &                          ip*nprocz,ip,
     &                          mpi_comm_world,status1,ierror)
                           
                           
                        end if
                        
                     else
c
c     keep what you have in txys, will later be overwritten
c
                     end if
                  end if

                  call mpi_barrier(mpi_comm_world,ierror)


c     
c     Divide accumulated statistics by sum of weights
c     and write velocity statistics
c     
                  if (my_node_world.eq.0) then
                     if (yb.le.nyp) then
                        do x=1,nx
                           txys(x)=txys(x)*(1./sumw)
                           wxy(x,yb) = txys(x)
                        end do
                     end if
                  end if
               end do
            end do
         end if
#endif

         if (nproc.eq.1) then
c
c     Put own statistics on wxy
c
            do yb=1,nyp
               do x=1,nx
                  wxy(x,yb)=xys(x,yb,i)
               end do   
            end do
c     
c     Divide accumulated statistics by sum of weights
c     and write velocity statistics
c     
            do y=1,nyp
               do x=1,nx
                  wxy(x,y)=wxy(x,y)*(1./sumw)
               end do
            end do  
         end if

         if (my_node_world.eq.0) then
            write(19) wxy
         end if
      end do
c
c     Loop over scalars and scalar statistics
c
      do ith=1,scalar
         do i=1,nxysth
#ifdef MPI
            if (nproc.gt.1) then
               do y=1,nyp/nprocz+1
                  do ip=0,nprocz-1
                     yb=ip+1+(y-1)*nprocz
                     txys=0.
                     do iii=0,nprocx-1
                        if (yb.le.nyp) then
c     
c     Send the individual statistics to processor 0
c     
                           if (my_node_world.eq.ip*nprocx+iii
     &                          .and.(ip+iii).ne.0) then
                              call mpi_ssend(xysth(1,y,i,ith),nx,
     &                             mpi_double_precision,0,
     &                             ip*nprocx+iii+100,
     &                             mpi_comm_world,ierror)
                              
                           end if
c     
c     Receive individual statistics in wxys and put on wxy
c     
                           if (my_node_world.eq.0) then
                              if ((ip+iii).ne.0) then
                                 call mpi_recv(temp,nx,
     &                                mpi_double_precision,
     &                                ip*nprocx+iii,ip*nprocx+iii+100,
     &                                mpi_comm_world,status1,ierror)
                              else
                                 do x=1,nx
                                    temp(x)=xysth(x,y,i,ith)
                                 end do
                              end if
                              do x=1,nx
                                 txys(x)=temp(x)+txys(x)
                              end do
                           end if
                        end if

c                        call mpi_barrier(mpi_comm_world,ierror)

                     end do
c     
c     Divide accumulated statistics by sum of weights
c     and write velocity statistics
c     
                     if (my_node_world.eq.0) then
                        if (yb.le.nyp) then
                           do x=1,nx
                              txys(x)=txys(x)*(1./sumw)
                              wxy(x,yb) = txys(x)
                           end do
                        end if
                     end if
                  end do
               end do
            end if
#endif
            if (nproc.eq.1) then
c
c     Put own scalar statistics on wxy
c
               do yb=1,nyp
                  do x=1,nx
                     wxy(x,yb)=xysth(x,yb,i,ith)
                  end do   
               end do
c     
c     Divide accumulated statistics by sum of weights
c     and write scalar statistics
c     
               do y=1,nyp
                  do x=1,nx
                     wxy(x,y)=wxy(x,y)*(1./sumw)
                  end do
               end do  
            end if
            if (my_node_world.eq.0) then
               write(19) wxy
            end if
         end do
      end do
c     
c     Close statistics file
c     
      if (my_node_world.eq.0) then
         close(unit=19)
      end if
c      
c     Two-point correlations in z
c     ***************************
c
      if (corrf) then
         if (my_node_world.eq.0) then
            write(*,*) '** doing two-point correlations in z'
            open(unit=19,file=corrnam,form='unformatted')
            rewind(19)
c     
c     Write header
c     
            write(19) -re*dstar,.false.,xl/dstar,zl/dstar,t/dstar,0.,
     &           (pr(i),m1(i),i=1,scalar)
            write(19) nx,nyp,nzc,nfzsym
            write(19) fltype,dstar
            if (fltype.ge.6) write(19) bstart/dstar,bslope/dstar
            write(19) sumw/dstar,scalar
            write(19) ncorr
            
         end if
c     
c     Send the individual statistics to processor 0
c     
         do i=1,ncorr
c
c     if having saved a whole x plane...
c
            if (corrxi(i).eq.0) then
               sumww = sumw*nx
            else
               sumww = sumw
            end if

            if (nproc.eq.1) then
               do z=1,nzc+2
                  totcorr(z) = corr(z,i)
               end do
               do z=1,4
                  totcorr_ms(z) = corr_ms(z,i)
               end do
            else 
               if (xb.eq.0) then
#ifdef MPI
                  call mpi_reduce(corr(1,i),totcorr,nzc+2,
     &                 mpi_double_precision,mpi_sum,0,my_comm_z,
     &                 ierror)
                  call mpi_reduce(corr_ms(1,i),totcorr_ms,4,
     &                 mpi_double_precision,mpi_sum,0,my_comm_z,
     &                 ierror)
#endif
               end if
            end if
           
            if (my_node_world.eq.0) then
c     
c     Postprocess and write the correlations
c     
c     
c     Do the backward transform
c     
               do j=1,nzc+2
                  cd(j)=totcorr(j)
               end do
               call vrfftb(cd(1),cd(2),cw(1),cw(2),
     &              nzc,1,2,1,corrwsave)
c
c     Write coordinate info
c
               write(19) corrx(i),corry(i),
     &              corrxi(i),corryi(i),
     &              corrxx(i),corryy(i),corrt(i)
c
c     Write mean and square
c
               write(19) (totcorr_ms(z)*(1./sumww),z=1,4)
c
c     And write data
c
               write(19) (cd(z)*(1./sumww),z=1,nzc)
            
            end if
         end do
         
         if (my_node_world.eq.0) then
            close(19)
         end if
      end if


c      
c     Two-point correlations in x
c     ***************************
c
      if (corrf_x) then
         if (my_node_world.eq.0) then
            write(*,*) '** doing two-point correlations in x'
            open(unit=19,file=corrnam_x,form='unformatted')
            rewind(19)
c     
c     Write header
c     
            write(19) -re*dstar,.false.,xl/dstar,zl/dstar,t/dstar,0.,
     &           (pr(i),m1(i),i=1,scalar)
            write(19) nx,nyp,nzc,nfzsym
            write(19) fltype,dstar
            if (fltype.ge.6) write(19) bstart/dstar,bslope/dstar
            write(19) sumw/dstar,scalar
            write(19) ncorr_x
            
         end if
c     
c     Send the individual statistics to processor 0
c     
         do i=1,ncorr_x
            if (corrzi(i).eq.0) then
               sumww = sumw*nz
            else
               sumww = sumw
            end if

            if (nproc.eq.1) then
               do z=1,nx+2
                  totcorr_x(z) = corr_x(z,i)
               end do
               do z=1,4
                  totcorr_ms_x(z) = corr_ms_x(z,i)
               end do
            else 
c               if (xb.eq.0) then
#ifdef MPI
                  call mpi_reduce(corr_x(1,i),totcorr_x,nx+2,
     &                 mpi_double_precision,mpi_sum,0,mpi_comm_world,
     &                 ierror)
                  call mpi_reduce(corr_ms_x(1,i),totcorr_ms_x,4,
     &                 mpi_double_precision,mpi_sum,0,mpi_comm_world,
     &                 ierror)
#endif
c               end if
            end if
           
            if (my_node_world.eq.0) then
c     
c     Postprocess and write the correlations
c     
c     
c     Do the backward transform
c     
               do j=1,nx+2
                  cd_x(j)=totcorr_x(j)
               end do
               call vrfftb(cd_x(1),cd_x(2),cw_x(1),cw_x(2),
     &              nx,1,2,1,corrwsave_x)
c
c     Write coordinate info
c
               write(19) corrz(i),corry_x(i),
     &              corrzi(i),corryi_x(i),
     &              corrzz(i),corryy_x(i),corrt_x(i)
c
c     Write mean and square
c
               write(19) (totcorr_ms_x(z)*(1./sumww),z=1,4)
c
c     And write data
c
               write(19) (cd_x(z)*(1./sumww),z=1,nx)
            
            end if
         end do
         
         if (my_node_world.eq.0) then
            close(19)
         end if
      end if





c
c     Time Series
c     ***********
c
      if (serf) then

         if (my_node_world.eq.0) then
            write(*,*) '** doing time series'
         end if

         if (sercount.gt.nsamp) then
            if (my_node_world.eq.0) then
               write(*,*) 'sercount larger than nsamp'
               write(*,*) sercount,nsamp
            end if
            call stopnow(165645)
         end if
c
c     Collect all series in totseries
c
         if (nproc.eq.1) then
            do i=0,nser
               do j=1,msamp
                  totseries(j,i) = series(j,i)
               end do
            end do
         else
#ifdef MPI
            call mpi_reduce(series(1,1),totseries(1,1),msamp*nser,
     &           mpi_double_precision,mpi_sum,0,mpi_comm_world,
     &           ierror)
            do i=1,sercount
               totseries(i,0) = series(i,0)
            end do
#endif
         end if
c
c     Write series to disk
c
         if (my_node_world.eq.0) then
            open(unit=19,file=namser,form='formatted',position='append')
            do i=1,sercount
               totseries(i,0) = totseries(i,0)/dstar
               write(19,'(4000e25.16e3)') (totseries(i,j),j=0,nser)
            end do
            close(19)
         end if
      end if
      
      end subroutine wxys
