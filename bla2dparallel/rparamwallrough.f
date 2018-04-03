c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/rparamwallrough.f $
c $LastChangedDate: 2007-10-29 14:23:13 +0100 (Mon, 29 Oct 2007) $
c $LastChangedBy: qiang@MECH.KTH.SE $
c $LastChangedRevision: 786 $
c
c ***********************************************************************
      subroutine rparamwallrough(wallroughinput,updat,every,
     &     monitor,monfile,v_wall,taylor4,h_rough,hstart,hend,hrise,
     &     hfall,zfunc,qq,psi,rghfil,roughfile,pert,my_node)
c
c     Reads in the parameters needed for the modification
c     (projection based on Taylor expansion) of the boundary
c     conditions at the "wall", if there is surface roughness.
c
c     The file which to read from is specified in rparambl.f.
c
c     Called by bla.f
c
      implicit none

      include 'par.f'

      integer every,zfunc,my_node
      real h_rough,hstart,hend,hrise,hfall,qq,psi
      logical updat,v_wall,taylor4,rghfil,monitor,pert
      character(80) wallroughinput,monfile,roughfile
c
c     Initialization
c
      updat=.false.
      every=0
      monitor=.false.
      v_wall=.false.
      taylor4=.false.
      h_rough=0.
      hstart=0.
      hend=0.
      hrise=0.
      hfall=0.
      zfunc=0
      qq=0
      rghfil=.false.


      if (my_node.eq.0) then
         write(*,*)
         write(*,*) '--> Wall roughness included. Opening  ',
     &        trim(wallroughinput)
         write(*,*)
      end if
c
c     Read the wall-roughness parameters
c
      open(unit=15,file=wallroughinput,status='old')

      read(15,*) updat

      if (updat) then
         if (pert) then
            updat=.false.
            if (my_node .eq. 0) then
               write(*,*)
               write(*,*) '--> Note: Update of roughness conditions' 
               write(*,*) '    not possible in perturbation mode'
               write(*,*)
            end if
         end if

         read(15,*) every
         if (.not.pert) then
            if (my_node.eq.0) then
               write(*,*)
               write(*,*) 'Wall-roughn. BCs updated every ',
     &              every,'th dt'
               write(*,*)       
            end if
         end if

         read(15,*) monitor
         if (monitor) read(15,*) monfile
      end if

      read(15,*) v_wall
      read(15,*) taylor4
      read(15,*) h_rough
      read(15,*) hstart
      read(15,*) hend
      read(15,*) hrise
      read(15,*) hfall
      read(15,*) zfunc

      if (zfunc .eq. 1) then
         read(15,*) qq
         read(15,*) psi
      end if

      read(15,*) rghfil
      if (rghfil) read(15,*) roughfile

      close(15)

      
      end subroutine  rparamwallrough
