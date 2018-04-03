c ***********************************************************************
c
c $HeadURhttps://www2.mech.kth.se/svn/simson/trunk/bls/bls.f $
c $LastChangedDate: 2007-11-23 22:05:23 +0100 (Fri, 23 Nov 2007) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 975 $
c
c ***********************************************************************
      program bls
c
c     Start program for boundary layer simulations.
c     A temperature field is written along with the velocity
c     components. The file format of the input table has also changed and
c     should be produced by the program fsc.f in this directory.
c
c     Note!
c     The localized perturbation earlier generated in sta now corresponds
c     to ditype=4 and that ipoly=4 in sta corresponds to ipoly=5
c
c     Flow types:
c     -3 Asymptotic suction boundary layer???
c     -2 Temporal Falkner-Skan-Cooke
c     -1 Temporal Falkner-Skan
c      1 Temporal Poiseuille flow
c      2 Temporal Couette flow
c      3 Temporal Blasius
c      4 Spatial Poiseuille flow
c      5 Spatial Couette flow
c      6 Spatial Blasius
c      7 Spatial Falkner-Skan
c      8 Spatial Falkner-Skan-Cooke
c      9 Spatial parallel Blasius/Falkner-Skan/Falkner-Skan-Cooke
c     10 Asymptotic suction boundary layer (with fringe)???
c     11 Wall jet flow???
c
      implicit none

      include 'par.f'
c
c     Main storage
c     
c     Q. Why is it 3+scalar?         
      complex ur(memnx,memny,memnz,3+scalar)
c
c     Step 2 storage
c
      real u2(nx+2,nz,3+scalar),w(nx+2,nz)
      complex pxy(nx/2,nyp)
      complex urx(nx/2)
c
c     xy storage
c
      complex uxy(nx/2,nyp,3,2)
c
c     File parameters
c
      character*80 namnut,namnin
      logical filein
c
c     Geometrics
c
      real alfa(nx/2),beta(nz),eta(nyp),deta(nyp),xl,zl,xs
      real pi
      parameter (pi = 3.1415926535897932385)
c
c     FFT preprocessing data
c
      real prex(nx+15),prey(nyp*2+15)
      real prez(nz*2+15)
c
c     Disturbance parameters
c
      integer ipoly,ith
      real amp,theta,xscale,zscale,yscale,xloc0
      logical locdi
c
c     Type of disturbance
c
      integer ditype
c
c     Noise parameters
c
      integer nxn,nyn,nzn,seed
      real ed
      logical noise
c
c     Time
c
      real t
c
c     Flow
c
      integer fltype
      logical spat
      real re,thgrad(scalar)
c
c     Loop indices etc
c
      integer x,y,z,i,j,k
      real sym
c
c     Coordinates
c
      real x1,y1,z1,xp,zp,rsquare,r1
c
c     Localized disturbance temporaries
c
      real w1,poly,dpoly,scale1,scalef
c
c     Oblique waves +- fi deg
c
      real ystart,yrise,yend,yfall,walfa,wbeta,k2
      real yfunc,dyfunc,step,dstep,energy
      logical waves
c
c     Gaussian initial condition
c
      real y0,yscale1
      real gaussy,dgaussy
      real gauss,dgauss
      logical gaussian
c
c     New parameters for BL
c
      real etab(nyp),xlb,zlb,reb,h2,dstar,dstar2,ushift
      real cubip
      integer nbla
      real fbla(mbla),dfbla(mbla),d2fbla(mbla),d3fbla(mbla),gbl(mbla)
      real dgbl(mbla),th(mbla,scalar),dth(mbla,scalar)
      real d2gbl(mbla),d2th(mbla,scalar)
      real dybla
c
c     Spatially growing boundary layer
c
      real x0,bstart,blength,u0low,u0lowin
      real bu(nx+2,nyp,3+scalar),wb(nx+2,nyp)
      real xbl1(nx),xbl2(nx)
      real bst(nx),bc(nx)
c
c     Asymptotic suction boundary layer
c
      logical asbl
      real ysuc,usuc(nyp)
c
c     Falkner-Skan-Cooke
c
      real rlam,spanv
      character*80 namflo
c
c     Input file data
c
      real rein,xlin,zlin,dstarin,bslin,bstin,rlamin,spanin
      integer flin
c
c     OS velocity eigenfunctions for temporal simulation.
c     eps2, eps3= disturbance amplitudes.
c
      real u2r(nyp),u2i(nyp),v2r(nyp),v2i(nyp)
      real u3r(nyp),u3i(nyp),v3r(nyp),v3i(nyp),w3r(nyp),w3i(nyp),eps2,wr
      real eps3,wi,alpha,betta,x2a,z2b,x3p,x3n,reos,dnx,dnz,dnxz
      logical os,physos
      integer nypin
      real m1(scalar),pr(scalar)
      integer poi_scal
      logical poi_zero
c
c     Perturbation equations
c
      logical pert,baseflow

      logical specm,pertfromfile
      integer n_specm,i_spec
      integer nalfa_s(10),nbeta_s(10)
      real    y0_s(10),yscale1_s(10),amp_s(10)

c
c     Perturbation from file  
c
      real u_real(nyp),u_imag(nyp),v_real(nyp),v_imag(nyp)
      real    w_real(nyp),w_imag(nyp),amp_f(10),phi_f(10)
      real temp_real(nyp),temp_imag(nyp)
      integer n_pertfromfile,nalfa_f(10),nbeta_f(10)
      complex init_cond(nyp,4),im
      character*80 initcond(10)
      integer scal_pos

c
c     Jose: Additional variables for channel flow with varying viscosity
c
      real vel(nyp),eta_file(nyp)
      real gr(scalar)
      integer ivarvisc
c
c     Jose: For base flow noise
c
      real vel_nois(nyp)
      real ranb, ran2, sumran, ranav
      integer :: bas_nois = -10254
      real :: bas_ed = 0.01
      logical :: ibas_nois = .false.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc Ritabrata: defining scalar flag. all changes made by me start with "rt_"
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c    nyp is the number of y points (ny). nx, ny, nz are defined in par.f which this file uses

      real rt_temp(nyp)
      real rt_salt(nyp)
      real rt_temperature_noise(nyp) ! I do not know temperature noise is really required

      real rt_pr(scalar)   ! is it really required? or has it been alreay defined          

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc Flags for switching on scalars, similar to ivarvisc. Sees if
ccccccc I am reading any scalar data from files 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer rt_tempflag
      integer rt_saltflag
      


      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*                             Simson        '//
     &     '                      *'
      write(*,*) '*                         bls $Rev: 975 $   '//
     &     '                      *'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*i*******************************************'//
     &     '***********************'
      write(*,*)

      if (nproc.gt.1) then
         write(*,*) 'Compile with nproc=1'
      end if
c
c     Default settings
c
      bstart=0.
      blength=0.
      rlam=0.0
      spanv=0.0
      asbl=.false.
      filein=.false.
      spat=.false.
      im=cmplx(0.,1.)
c
c     Reads parameters from bls.i
c
      open(unit=10,status='old',file='bls.i')
c
c     Stick in the input velocity field first if there is one
c
      read(10,9000) namnut
9000  format(a80)
      read(10,*,err=1010) reb
c      write(*,*) 'Here without the error in the bls.i file '
      write(*,*) 'Initial field (namnut)                : ',trim(namnut)
      write(*,*) 'Reynolds number (re)                  :',reb
      goto 1020
c
c     Two file names read namnut
c
c     The following lines until the 1020 line won't be read unless
c     there is an error
1010  continue
      backspace(10)
      filein=.true.
      write(*,*) 'Here with the error in the bls.i file '
      namnin=trim(namnut)
      read(10,9000) namnut
      write(*,*) 'Read field (namnin)                   : ',trim(namnin)

      read(10,*) reb
      write(*,*) 'Reynolds number (re)                  :',reb

1020  read(10,*) xlb
      write(*,*) 'Domain length (xlb)                   :',xlb

      read(10,*) h2
      write(*,*) 'Domain height (h2)                    :',h2

      read(10,*) zlb
      write(*,*) 'Domain width (zlb)                    :',zlb

      read(10,9000) namflo
      write(*,*) 'Similarity data file (namflo)         : ',trim(namflo)

      read(10,*) fltype
      write(*,*) 'Flow type (fltype)                    :',fltype

      if (fltype.eq.-1.or.fltype.eq.-2.or.
     &    fltype.eq. 7.or.fltype.eq. 8) then
         read(10,*) rlam
         write(*,*) 'Power law exp (rlam)                  :',rlam
      end if

      if (fltype.eq.-2.or.fltype.eq.8) then
         read(10,*) spanv
         write(*,*) 'Spanwise vel (spanv)                  :',spanv
      end if

      if (fltype.ge.4) then
         spat=.true.
         read(10,*) bstart
         write(*,*) 'Blending start (bstart)               :',bstart

         read(10,*) blength
         write(*,*) 'Blending length (blength)             :',blength
      end if

      if (fltype.eq.3) then
         read(10,*) asbl
         if (asbl) then
            write(*,*) 'Asymptotic suction boundary layer (asbl=T)'
         else
            write(*,*) 'Blasius boundary layer'
         end if
      end if
c
c     Perturbation mode
c
      read(10,*) pert
      write(*,*) 'Perturbation mode (pert)              :',pert
      baseflow=.false.
      if (filein) then
         baseflow=.true.
      else
         if (pert) baseflow=.true.
      end if
c
c     Shift velocity
c
      read(10,*) ushift
      write(*,*) 'Shift velocity (ushift)               :',ushift
c
c     Localized disturbances
c
      read(10,*) locdi
      write(*,*) 'Localized perturbation (locdi)        :',locdi
      if (locdi) then
         read(10,*) ditype
         write(*,*) '  Disturbance type (ditype)           :',ditype
         read(10,*) amp
         write(*,*) '  Amplitude (amp)                     :',amp
         read(10,*) theta
         write(*,*) '  (theta)                             :',theta
         read(10,*) xscale
         write(*,*) '  (xscale)                            :',xscale
         read(10,*) xloc0
         write(*,*) '  (xloc0)                             :',xloc0
         read(10,*) yscale
         write(*,*) '  (yscale)                            :',yscale
         if ( fltype.eq.1.or.fltype.eq.2.or.
     &        fltype.eq.4.or.fltype.eq.5 ) then
            if ( yscale.lt.1.or.yscale.gt.1. ) then
               write(*,*) 'Warning! The polynomial is not defined ',
     &              'between -1 and 1 anymore.'
            end if
         end if
         read(10,*) zscale
         write(*,*) '  (zscale)                            :',zscale
         read(10,*) ipoly
         write(*,*) '  (ipoly)                             :',ipoly
      end if
c
c     Gaussian disturbance
c
      read(10,*) gaussian
      write(*,*) 'Gaussian perturbation (gaussian)      :',gaussian
      if (gaussian) then
         read(10,*) amp
         read(10,*) y0
         read(10,*) yscale1
         read(10,*) walfa
         read(10,*) wbeta
         k2=walfa**2+wbeta**2
      end if
c
c     Wave perturbation
c
      read(10,*) waves
      write(*,*) 'Wave perturbation (waves)             :',waves
      if (waves) then
         read(10,*) energy
         read(10,*) ystart
         read(10,*) yend
         read(10,*) yrise
         read(10,*) yfall
         read(10,*) walfa
         if (ystart.gt.yend) then
            write(*,*) 'ystart greater than yend'
            stop
         end if
         if (ystart.lt.0) write(*,*) 'ystart < 0 '
         if (abs(2.*pi/xlb-walfa).gt.1e-8) then
            write(*,*) 'walfa dont fit the box walfa=',walfa
         end if
         read(10,*) wbeta
         if (abs(2.*pi/zlb-wbeta).gt.1e-8) then
            write(*,*) 'wbeta dont fit the box wbeta=',wbeta
         end if
         k2=walfa**2+wbeta**2
         amp=sqrt(energy*h2*2/((yrise+yfall)*0.4057052525+
     &        (yend-ystart-yrise-yfall)+(1/yrise+1/yfall)*
     &        1.638270581/k2))
         if (walfa.eq.0.or.wbeta.eq.0) amp=amp/sqrt(2.)
      end if
c
c     Read Orr-Sommerfeld modes
c
      read(10,*) os
      write(*,*) 'OS modes (os)                         :',os
      if (os) then
         open(77,status='old',file='os.i')
         read(77,*) reos,alpha,betta,wr,wi
         read(77,*) eps2,eps3,physos,nypin
         if (nypin.ne.nyp.or.reb.ne.reos) then
            write(*,*) 'Error: Re, Ny must be same in os.i and bls.i.'
            write(*,*) 'os.i: ',nypin,reos,'bls.i: ',reb,nyp
            stop
         end if
         if (nz.eq.1) eps3=0.
         if (abs(eps2).gt.0.) then
            do i=1,nyp
               read(77,*) u2r(i),u2i(i),v2r(i),v2i(i)
            end do
         end if
         if (abs(eps3).gt.0.) then
            do i=1,nyp
               read(77,*) u3r(i),u3i(i),v3r(i),v3i(i),w3r(i),w3i(i)
            end do
         end if
         close(77)
      end if
c
c     Spectral space mode
c
      read(10,*) specm
      write(*,*) 'Spectral space mode (specm)           :',specm
      if (specm) then
         read(10,*) n_specm
         if (n_specm.gt.10) then
            write(*,*)'n_specm is larger than 10:
     &           change in bls.f to allow this'
            stop
         end if
         do i=1,n_specm
            read(10,*) amp_s(i)
            read(10,*) y0_s(i)
            read(10,*) yscale1_s(i)
            read(10,*) nalfa_s(i)
            read(10,*) nbeta_s(i)
            write(*,*)'Adds to a selected pair of wavenumbers a'//
     &           ' prescribed'
            write(*,*)'perturbation that has a gaussian shape,'//
     &           ' modified for B.C.'//
     &           ' The chosen wavenumbers are:'
            write(*,*)'alpha=',nalfa_s(i),'and beta=',nbeta_s(i)

         end do
      end if
c
c     Read perturbation from file
c
      read(10,*) pertfromfile
      write(*,*) 'Initial cond. from file (pertfromfile):'
     &     ,pertfromfile
      if (pertfromfile) then
         read(10,*) n_pertfromfile
         if (n_pertfromfile.gt.10) then
            write(*,*)'n_pertfromfile is larger than 10:
     &           change in bls.f to allow this'
            stop
         end if
         do i=1,n_pertfromfile
            read(10,*) initcond(i)
            read(10,*) amp_f(i)
            read(10,*) phi_f(i)
            read(10,*) nalfa_f(i)
            read(10,*) nbeta_f(i)
            write(*,*)'Adds to a selected pair of wavenumbers a'//
     &           ' prescribed'
            write(*,*)'perturbation that is read from text file. The'//
     &           ' chosen wavenumbers are:'
            write(*,*)'alpha=',nalfa_f(i),'and beta=',nbeta_f(i)
         end do
      end if
c
c     Noise perturbation
c
      noise=.false.
      read(10,*) noise
      write(*,*) 'Noise (noise)                         :',noise
      if (noise) then
         read(10,*) ed
         write(*,*) '  Energy density (ed)            :',ed
         read(10,*) nxn
         write(*,*) '  # x-modes (nxn)                :',nxn
         read(10,*) nyn
         write(*,*) '  # y-modes (nyn)                :',nyn
         read(10,*) nzn
         write(*,*) '  # z-modes (nzn)                :',nzn
         read(10,*) seed
         write(*,*) '  Random seed (seed)             :',seed
      end if

c
c     Varying Viscosity
c
      read(10,*) ivarvisc
      do ith = 1,scalar
c     rt_scal_fl(ith) - flag for reading initial profile of ith scalar
c     Keep this as 1 always
         read(10,*) rt_scal_fl(ith)
c     if rt_scal_fl.ne. 0, then read the name of the file providing the
c     initial scalar profile.
         if (rt_scal_fl(ith).ne.0) then            
            read(10,*) rt_scal_prof(ith)
         end if
c     rt_pr(ith) - Prandtl number of ith scalar         
         read(10,*) rt_pr(ith)
c     rt_ri(ith) - Richardson number of ith scalar
         read(10,*) rt_ri(ith)
      end do
c      if(ivarvisc.ne.0) then
c         do x = 1,scalar
c            read(10,*) pr(x)
c            read(10,*) gr(x)
c         end do
c      end if
      write(*,*)'  Varying viscosity               :',ivarvisc


      ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccc Ritabrata scalar reading
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! reading temperature from file "rt_temp.dat"
      do ith = 1,scalar
         if (rt_tempflag.ne.0)then
            open(54, status='old', file=rt_scal_prof(ith)) ! what is 'old'?
            do y = 1,nyp
               read(54,*)eta_scal_file(y,ith),scal_from_file(y,ith)
            end do
            close(54)
         end if
      end do

! reading salinity from file "rt_sal.dat"
         
c      if (rt_saltflag.ne.0)then
c         open(64, status='old', file='rt_sal.dat') 
c         do y = 1,nyp
c            read(64,*)eta_s_file(y),s_from_file(y)
c         end do
c         close(64)
c      end if
      
                
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



c
c     Poiseuille flow
c
c Jose: Reading in the ivarvisc variable from bls.i

c         ivarvisc = 1
c
c Jose: varying viscosity
c
      if(ivarvisc.ne.0) then
         open(36,status='old',file='vel.dat')
         do y=1,nyp
            read(36,*)eta_file(y),vel(y)
c     write(*,*)'simson y ', eta(y), ' input_file y',
c     &              eta_file(y), ' Velocity value', vel(y)
         end do
         close(36)
      end if
      
      t=0.
      xs=0.
      u0low=-ushift
      dstar=2./h2
      re=reb/dstar
      xl=xlb*dstar
      zl=zlb*dstar
c
c     Read input file namnin if needed
c
      if (filein) then
         call rdiscp(ur,rein,xlin,zlin,t,xs,
     &        dstarin,flin,bstin,bslin,rlamin,spanin,spat,namnin,3,pxy,
     &        urx)
         if (dstarin.ne.dstar) then
            write(*,*) 'dstarin not equal to dstar'
            stop
         end if
         if (rein.ne.re) then
            write(*,*) 'rein not equal to re'
            stop
         end if
         if (xlin.ne.xl) then
            write(*,*) 'xlin not equal to xl'
            stop
         end if
         if (zlin.ne.zl) then
            write(*,*) 'zlin not equal to zl'
            stop
         end if
      end if
c
c     Blending function
c
      bstart=bstart*dstar
      blength=blength*dstar
c
c     Initialize coordinates
c
      do x=1,nx/2
         alfa(x)=2.*pi/xl*real(x-1)
      end do
      beta(1)=0.
      do z=2,nz/2+1
         beta(z)=2.*pi/zl*real(z-1)
         beta(nz+2-z)=-2.*pi/zl*real(z-1)
      end do
      do y=1,nyp
         eta(y)=cos(pi*real(y-1)/real(nyp-1))
      end do
      deta(1)=(1.-eta(2))*.5
      deta(nyp)=deta(1)
      do y=2,nyp-1
         deta(y)=(eta(y-1)-eta(y+1))*.5
      end do
      if ( fltype.eq.-2.or.fltype.eq.-1.or.fltype.eq.3.or.
     &     fltype.eq. 6.or.fltype.eq. 7.or.fltype.eq.8.or.
     &     fltype.eq. 9) then
         do y=1,nyp
            etab(y)=(eta(y)+1.)/dstar
         end do
      else
         do y=1,nyp
            etab(y)=eta(y)
         end do
      end if
      fbla = 0.
c
c     Initialize x, y and z-transforms
c
      call vrffti(nx,prex,0)
      call vcosti(nyp,prey,0)
      call vcffti(nz,prez,0)
c
c     Read velocity profile and displacement thickness coefficient
c
      if ( fltype.eq.-2.or.fltype.eq.-1.or.fltype.eq.3.or.
     &     fltype.eq. 6.or.fltype.eq. 7.or.fltype.eq.8.or.
     &     fltype.eq. 9) then
         if (.not.asbl) then
            call rbla(fbla,dfbla,d2fbla,d3fbla,gbl,dgbl,d2gbl,th,dth,
     &           d2th,nbla,dybla,rlam,dstar2,mbla,
     &           namflo,fltype,m1,pr,scalar,thgrad)
         end if
c
c     Compute base flow bu
c
         if (spat.and.fltype.ne.9) then
c
c     Spatially growing boundary layer
c
            x0=dstar**2*re*(rlam+1.)/(dstar2**2*2.)
            call blvel(bu,xs,x0,eta,dybla,dfbla,gbl,th,nbla,rlam,spanv,
     &           re,xl,bstart,blength,u0low,prex,prey,nx,nyp,wb,
     &           xbl1,xbl2,bst,bc,m1,scalar,thgrad,pr)
         else
c
c     Parallel boundary layer
c     If asymptotic suction boundary layer flow (2)
c     compute analytical solution
c
            if (asbl) then
               do y=1,nyp
                  ysuc=(1.-cos(pi*float(y-1)/float(nyp-1)))*h2/2.
                  usuc(y)=1.-exp(-ysuc)
               end do
               do y=1,nyp
                  do x=1,nx
                     bu(x,y,1)=usuc(nyp+1-y)+u0low
                     bu(x,y,2)=0.
                     bu(x,y,3)=0.
                  end do
               end do
            else
               do y=1,nyp
                  bu(1,y,1)=cubip(etab(y)*dstar2,dfbla,dybla,nbla)+u0low
                  bu(1,y,2)=0.
                  bu(1,y,3)=spanv*cubip(etab(y)*dstar2,gbl,dybla,nbla)
                  do j=1,scalar
                     bu(1,y,j+scalar)=
     &                    cubip(etab(y)*dstar2,th(1,j),dybla,nbla)
                  end do
               end do
               bu(1,nyp,1)=u0low
               do y=1,nyp
                  do x=1,nx
                     bu(x,y,1)=bu(1,y,1)
                     bu(x,y,2)=0.
                     bu(x,y,3)=bu(1,y,3)
                     do j=1,scalar
                        bu(x,y,3+j)=bu(1,y,3+j)
                     end do
                  end do
               end do
            end if
         end if
      else if ( fltype.eq.1.or.fltype.eq.4 ) then

c     Jose: Some random perturbation in the base flow itself.
c
         vel_nois = 0.0
         if(ibas_nois) then
            ranb = 0.0
            sumran = 0.0
            do y=2,nyp-1
               ranb = ran2(bas_nois)
               sumran = sumran + ranb
c     write(*,*) ranb
               vel_nois(y) = ranb
            end do
            ranav = sumran/real(nyp)

            do y = 2,nyp-1
               vel_nois(y) = vel_nois(y) - ranav
c            write(*,*) vel_nois(y)
            end do
            
         end if

         do y=1,nyp

            do x=1,nx
c
c     Parabolic velocity profile
c
               if(ivarvisc.eq.1) then
                  bu(x,y,1) = vel(y) + bas_ed*vel_nois(y)
               else
                  bu(x,y,1)=1.0 - eta(y)**2 + u0low + bas_ed*vel_nois(y)
c                  bu(x,y,1)=1.0-eta(y)**2+u0low
c               bu(x,y,1)=1.0-eta(y)**2-.375
               end if
               bu(x,y,2)=0.0
               bu(x,y,3)=0.0
            end do
c            write(*,*) bu(16,y,1), vel(y)
         end do

      else if ( fltype.eq.2.or.fltype.eq.5 ) then
c
c     Couette flow
c
         do y=1,nyp
            do x=1,nx
               bu(x,y,1)=eta(y)+u0low
               bu(x,y,2)=0.0
               bu(x,y,3)=0.0
            end do
         end do
      end if

c
c     Scalar profiles and parameters
c
      poi_zero = .false.
c     Jose: Check below for definition
      poi_scal = 1
      if (scalar.gt.1) then
         write(*,*) '***********************************'
         write(*,*) 'Scalar profiles for Poiseuille flow'
         write(*,*) 'poi_zero = ',poi_zero
         write(*,*) 'poi_scal = ',poi_scal
         write(*,*) 'Please change directly in bls.f'
         write(*,*) '***********************************'
      end if
      do ith=1,scalar
c     Jose: To read the Prandtl number from the bls.i file
c Jose: I had earlier defined Pran in par.f
c     if(ivarvisc.ne.0) then
c     read(10,*) pr(ith)
c     read(10,*) gr(ith)
c     else
c     pr(1) = 0.01
c     endif
         
c     if (scalar.eq.4) then
c     pr(2) = 0.1
c     pr(3) = 0.02
c     pr(4) = 0.01
c     end if
         pr(ith) = rt_pr(ith)
         gr(ith) = rt_ri(ith)
         
         m1(ith) = 0.
         do y=1,nyp
            do x=1,nx
               if (poi_zero) then
c     
c     Velocity zero
c     
                  bu(x,y,1) = 0.
               end if
               if (poi_scal.eq.0) then
c     
c     Scalar zero
c     
                  bu(x,y,3+ith) = 0.
               else if (poi_scal.eq.1) then
c     
c     Linear scalar profile, lower wall 0, upper wall 1
c     
c     Jose: enforce the same boundary conditions
                  bu(x,y,3+ith) = 1.-(1.-eta(y))/2.
               else if (poi_scal.eq.2) then
c
c     Linear scalar profile, lower wall 1, upper wall 0
                  bu(x,y,3+ith) = (1.-eta(y))/2.
c     else if(poi_scal.eq.3) then
c
c     Jose:    Linear scalar profile, lower (-1), upper (0)
c     
c     bu(x,y,3+ith) = -(eta(y)-1.)/2.
               else
                  write(*,*) 'poi_scal not between 0 and 3. Stop.'
                  stop
               end if
            end do
         end do
      end do

      close(unit=10)

c
c     Initial velocity field:
c
c     For Blasius flow :        u_0 =  Blasius profile+u0low, delta star =1
c                               v_0 = -dpsi/dz
c                               w_0 =  dpsi/dy
c
c                               For zscale=0 the disturbance is two dimensional:
c                               u_0 =  Blasius profile+u0low+dpsi/dy
c                               v_0 = -dpsi/dx
c                               w_0 =  0.0
c                               where psi=amp*p(y/yscale)*exp(-(x/xscale)^2)
c
c     For Poiseuille flow :     u_0 =  1 - y*y
c                               v_0 =  dpsi/dz
c                               w_0 = -dpsi/dy
c
c     For Couette flow :        u_0 =  y
c                               v_0 =  dpsi/dz
c                               w_0 = -dpsi/dy
c
c     where
c                               psi =  amp*x/xscale*z/zscale*p(y/yscale)^2*
c                                      *exp(-(x/xscale)^2-(z/zscale)^2)
c
c     is one of the available disturbance types (ditype)
c
c     ipoly is the choice of polynom :
c                          1 : p(y) = (1-y*y)^2 => sym v
c                          2 : p(y) = y(1-y*y)^2 => antisym v
c                          3 : p(y) = (1+4y)(1-y*y)^2 => nonsym v
c                          4 : p(y) = y^3*exp(-y^2)
c                          5 : p(y) = (1+y)^2(1-y)^5 => focused toward wall y=0
c                          6 : p(y) = (1+y)^2(1-y)^2
c
c     The disturbance is rotated the angle theta about the y-axis ie :
c                                x  =  x'*cos(theta) - z'*sin(theta)
c                                z  =  x'*sin(theta) + z'*cos(theta)
c                                u' =  w*sin(theta)+bu(y)
c                                w' =  w*cos(theta)
c
c     (To enable larger timesteps the walls can be moved backwards
c     at u=-.375 for Poiseuille flow)
c

      if (filein) then
         call getxzp(u2(1,1,1),nyp,1,ur,sym)
         u0lowin=u2(1,1,1)
      end if
c

cc
cc Ritabrata:
cc


c

c     Generate initial bla field
c
      do y=1,nyp
         if (y.eq.nyp/100+2) then
            write(*,*)
            write(*,*) 'Generated 1%'
         end if
         if (y.eq.nyp/10+2) write(*,*) 'Generated 10%'
         if (y.eq.nyp/2+2) write(*,*) 'Generated 50%'
         if (y.eq.nyp*9/10+2) then
            write(*,*) 'Generated 90%'
            write(*,*)
         end if
         if (filein) then
            do i=1,3
               sym=1.
               if (i.eq.3) sym=-1.
               call getxzp(u2(1,1,i),y,i,ur,sym)
               if (i.eq.1) u2(1,1,1)=u2(1,1,1)-u0lowin+u0low
               call vcfftb(u2(1,1,i),u2(2,1,i),
     &              w,w(2,1),nz,nx/2,nx+2,2,prez)
               call vrfftb(u2(1,1,i),u2(2,1,i),
     &              w,w(2,1),nx,nz,2,nx+2,prex)
            end do
         end if
c
c     Clean u2 variable
c
         if (pert) then
c$$$            do i=1,3
c$$$               do j=1,nz
c$$$                  do k=1,nx+2 
c$$$                     u2(k,j,i)=0.0
c$$$                  end do
c$$$               end do
c$$$            end do
            u2 = 0.0
         end if

         if (locdi) then
            y1=etab(y)/yscale
            scale1 = sqrt(2.0)*exp(0.5)
            scalef = (2.0/3.0)**(1.5)*exp(1.5)
            if (ipoly.eq.1) poly=(1.0-y1*y1)**2
            if (ipoly.eq.1) dpoly=-4.0*y1*(1.0-y1*y1)/yscale
            if (ipoly.eq.2) poly=y1*(1.0-y1*y1)**2
            if (ipoly.eq.2) dpoly=1./yscale*(1.0-y1*y1)*(1.0-5.0*y1*y1)
            if (ipoly.eq.3) poly=(1.0+4.0*y1)*(1.0-y1*y1)**2
            if (ipoly.eq.3) dpoly=1./yscale*(4.0-4.0*y1-20.0*y1*y1)*
     &                            (1.0-y1*y1)
            if (ipoly.eq.4) poly=y1**3*exp(-y1*y1)
            if (ipoly.eq.4) dpoly=1./yscale*(3.*y1*y1-2.*y1**4)*
     &                            exp(-y1*y1)
            if (ipoly.eq.5) poly=(1.0+y1)**2*(1.0-y1)**5
            if (ipoly.eq.5) dpoly=1./yscale*2.*(1.0+y1)*(1.0-y1)**5
     &           -5.*(1.0+y1)**2*(1.0-y1)**4
            if (ipoly.eq.6) poly=(1.0+y1)**2*(1.0-y1)**2
            if (ipoly.eq.6) dpoly=1./yscale*2.*(1.0+y1)*(1.0-y1)**2
     &           -2.*(1.0+y1)**2*(1.0-y1)
c
c     Vortex pair
c
            if (ditype.eq.1) then
               do z=1,nz
                  if (zscale.ne.0) then
                     zp=zlb*real(z-nz/2-1)/real(nz)/zscale
                  end if
                  do x=1,nx
                     xp=xlb*real(x-nx/2-1)/real(nx)/xscale
                     if (.not.baseflow) then
                        do j=1,3+scalar
                           u2(x,z,j)=bu(x,y,j)
                        end do
                     end if
                     if (zscale.gt.0) then
                        x1=(xp*cos(theta)-zp*sin(theta))-xloc0/xscale
                        z1=xp*sin(theta)+zp*cos(theta)
                        w1=amp*x1*z1*dpoly*exp(-x1**2-z1**2)*
     &                       scale1*scalef
                        u2(x,z,1)=u2(x,z,1)+w1*sin(theta)*xscale

                        u2(x,z,2)=u2(x,z,2)-
     &                       amp*x1*poly*(1.0-2.0*(z1)**2)*
     &                       exp(-x1**2-z1**2)*scale1*scalef
                        u2(x,z,3)=u2(x,z,3)+w1*cos(theta)*zscale
                     else
                        xp=xp-xloc0/xscale
                        u2(x,z,1)=u2(x,z,1)+amp*dpoly*exp(-xp*xp)
     &                       *scale1*scalef
                        u2(x,z,2)=u2(x,z,2)+amp*poly*2.*xp/xscale*
     &                       exp(-xp*xp)*scale1*scalef
                     end if
                  end do
               end do
            end if
c
c     Axi-symmetric
c
            if (ditype.eq.2) then
               zscale = xscale
               do z=1,nz
                  zp=zlb*real(z-nz/2-1)/real(nz)/zscale
                  do x=1,nx
                     xp=xlb*real(x-nx/2-1)/real(nx)/xscale
                     if (.not.baseflow) then
                        do j=1,3+scalar
                           u2(x,z,j)=bu(x,y,j)
                        end do
                     end if
                     x1=xp-xloc0/xscale
                     z1=zp
                     rsquare = x1**2 + z1**2
                     r1 = sqrt(rsquare)
                     w1 = -amp*dpoly*exp(-rsquare)/2.0*scalef
                     u2(x,z,1)=u2(x,z,1)+w1*x1
                     u2(x,z,2)=u2(x,z,2)+amp*poly*(1.0-(r1)**2)*
     &                    exp(-rsquare)*scalef
                     u2(x,z,3)=u2(x,z,3)+w1*z1
                  end do
               end do
            end if
c
c     Wave packet
c
            if (ditype.eq.3) then
               do z=1,nz
                  zp=zlb*real(z-nz/2-1)/real(nz)/zscale
                  do x=1,nx
                     xp=xlb*real(x-nx/2-1)/real(nx)/xscale
                     if (.not.baseflow) then
                        do j=1,3+scalar
                           u2(x,z,j)=bu(x,y,j)
                        end do
                     end if
                     x1=xp-xloc0/xscale
                     z1=zp
                     w1 = amp*exp(-x1**2-z1**2)*scalef
                     u2(x,z,1)=u2(x,z,1)-w1*dpoly*x1
                     u2(x,z,2)=u2(x,z,2)+w1*poly*(1.0 -
     &                    2.0*x1**2)
                     u2(x,z,3)=u2(x,z,3)+0.0
                  end do
               end do
            end if
c
c     Perturbation type imported from sta.
c     ipoly=4 in sta now corresponds to ipoly=5 in bls
c
            if (ditype.eq.4) then
               do z=1,nz
                  zp=zlb*real(z-nz/2-1)/real(nz)
                  do x=1,nx
                     xp=xlb*real(x-nx/2-1)/real(nx)
                     if (.not.baseflow) then
                        do j=1,3+scalar
                           u2(x,z,j)=bu(x,y,j)
                        end do
                     end if
                     x1=xp*cos(theta)-zp*sin(theta)-xloc0/xscale
                     z1=xp*sin(theta)+zp*cos(theta)
                     if (ipoly.eq.5) then
                        w1=-amp*exp(-(x1**2+z1**2)/(xscale**2))
                     else
                        w1=-amp*x1/xscale*z1*dpoly*
     &                       exp(-(x1/xscale)**2-(z1/zscale)**2)
                     end if
                     if (ipoly.eq.5) then
                        u2(x,z,1)=u2(x,z,1)+w1*0.5*x1*dpoly
                        u2(x,z,2)=u2(x,z,2)-w1*poly*
     &                       (1.0-(x1**2+z1**2)/xscale**2)
                        u2(x,z,3)=u2(x,z,3)+w1*0.5*z1*dpoly
                     else
                        u2(x,z,1)=u2(x,z,1)+w1*sin(theta)
                        u2(x,z,2)=u2(x,z,2)+amp*x1/xscale*poly*(1.0-2.0*
     &                       (z1/zscale)**2)*exp(-(x1/xscale)**2
     &                       -(z1/zscale)**2)
                        u2(x,z,3)=u2(x,z,3)+w1*cos(theta)
                     end if
                  end do
               end do
            end if
         else if (waves) then
            y1=etab(y)
            yfunc=step((y1-ystart)/yrise)
     &           -step((y1-yend)/yfall+1)
            dyfunc=dstep((y1-ystart)/yrise)/yrise
     &           -dstep((y1-yend)/yfall+1)/yfall
            do z=1,nz
               zp=zlb*real(z-nz/2-1)/real(nz)
               do x=1,nx
                  xp=xlb*real(x-nx/2-1)/real(nx)
                  if (.not.baseflow) then
                     do j=1,3+scalar
                        u2(x,z,j)=bu(x,y,j)
                     end do
                  end if
                  u2(x,z,1)=u2(x,z,1)-amp*walfa*2*sin(walfa*xp)*
     &                 cos(wbeta*zp)*dyfunc/k2
                  u2(x,z,2)=u2(x,z,2)+amp*2*cos(walfa*xp)*cos(wbeta*zp)*
     &                 yfunc
                  u2(x,z,3)=u2(x,z,3)-amp*2*wbeta*cos(walfa*xp)*
     &                 sin(wbeta*zp)*dyfunc/k2
               end do
            end do
         else if (gaussian) then
            y1=etab(y)
            gaussy=(1./(sqrt(2.*pi)*yscale1))*
     &           exp(-(y1-y0)**2/(2.*(yscale1)**2))
            Dgaussy=-((y1-y0)/(sqrt(2.*pi)*yscale1**3))*
     &           exp(-(y1-y0)**2/(2.*(yscale1)**2))

            do z=1,nz
               zp=zlb*float(z-nz/2-1)/float(nz)
               do x=1,nx
                  xp=xlb*float(x-nx/2-1)/float(nx)
                  if (.not. baseflow) then
                     do j=1,3+scalar
                        u2(x,z,j)=bu(x,y,j)
                     end do
                  end if
                  u2(x,z,1)=u2(x,z,1)-amp*walfa*4*sin(walfa*xp)*
     &                 cos(wbeta*zp)*Dgaussy/k2
                  u2(x,z,2)=u2(x,z,2)+amp*4*cos(walfa*xp)*cos(wbeta*zp)*
     &                 gaussy
                  u2(x,z,3)=u2(x,z,3)-amp*4*wbeta*cos(walfa*xp)*
     &                 sin(wbeta*zp)*Dgaussy/k2
               end do
            end do
         else if (os) then
            dnx=xlb/real(nx)
            dnz=zlb/real(nz)
            dnxz=1./(real(nx)*real(nz))
            if (physos) then
               do z=1,nz
c     zp=dnz*real(z-1)                      !     0<z<lz
                  zp=dnz*real(z-nz/2-1) ! -lz/2<z<lz/2
                  do x=1,nx
c     xp=dnx*real(x-1)                    !     0<x<lx
                     xp=dnx*real(x-nx/2-1) ! -lx/2<x<lx/2
                     if (.not.baseflow) then
                        do j=1,3+scalar
                           u2(x,z,j)=bu(x,y,j)
                        end do
                     end if
                     x2a=alpha*xp
                     z2b=betta*zp
                     x3p=alpha*xp+betta*zp
                     x3n=alpha*xp-betta*zp
c
c     Initial velocity field:blasius flow +2-D TS+3-D pair of oblique waves
c
                     u2(x,z,1)=u2(x,z,1)
     &                    +eps2*(u2r(y)*cos(x2a)-u2i(y)*sin(x2a))
     &                    +2.*eps3*(u3r(y)*cos(x2a)-u3i(y)*sin(x2a))
     &                    *cos(z2b)
                     u2(x,z,2)=u2(x,z,2)
     &                    +eps2*(v2r(y)*cos(x2a)-v2i(y)*sin(x2a))
     &                    +2.*eps3*(v3r(y)*cos(x2a)-v3i(y)*sin(x2a))
     &                    *cos(z2b)
                     u2(x,z,3)=u2(x,z,3)
     &                    -2.*eps3*(w3r(y)*sin(x2a)+w3i(y)*cos(x2a))
     &                    *sin(z2b)
                  end do
               end do
            else
c
c     If (.not.physos), put initial conditions in the wave space
c
               if (filein) then
                  write(*,*)
     &                 'The combination os=.true., physos=.false. and'
                  write(*,*) 'filein is not implemented '
                  stop
               end if
               if (spat) then
                  write(*,*)
     &                 'The combination os=.true., physos=.false. and'
                  write(*,*) 'spat=.true. is not implemented '
                  stop
               end if
               do z=1,nz
                  do x=1,nx
                     u2(x,z,1)=0.
                     u2(x,z,2)=0.
                     u2(x,z,3)=0.
                  end do
               end do
c
c     Mean part(kx=kz=0)
c
               u2(1,1,1)=bu(1,y,1)
               u2(1,1,3)=bu(1,y,3)
c
c     Nonzero 2-D parts(kx=1,kz=0)
c
               u2(3,1,1)=eps2*u2r(y)
               u2(4,1,1)=eps2*u2i(y)
               u2(3,1,2)=eps2*v2r(y)
               u2(4,1,2)=eps2*v2i(y)
c
c     Nonzero 3-D parts
c     kx=kz=1 only
c
               if (abs(eps3).gt.0..and.nz.gt.1) then
                  u2(3,2,1)=eps3*u3r(y)
                  u2(4,2,1)=eps3*u3i(y)
                  u2(3,2,2)=eps3*v3r(y)
                  u2(4,2,2)=eps3*v3i(y)
                  u2(3,2,3)=eps3*w3r(y)
                  u2(4,2,3)=eps3*w3i(y)
c
c     kx=-kz=1 only
c
                  u2(3,nz,1)=eps3*u3r(y)
                  u2(4,nz,1)=eps3*u3i(y)
                  u2(3,nz,2)=eps3*v3r(y)
                  u2(4,nz,2)=eps3*v3i(y)
                  u2(3,nz,3)=-eps3*w3r(y)
                  u2(4,nz,3)=-eps3*w3i(y)
               end if
               call putxzp(u2(1,1,1),y,1,ur)
               call putxzp(u2(1,1,2),y,2,ur)
               call putxzp(u2(1,1,3),y,3,ur)
            end if
         else
c
c     Go here if no disturbances
c
            if (.not.baseflow) then
               do j=1,3+scalar
                  do z=1,nz
                     do x=1,nx
                        u2(x,z,j)=bu(x,y,j)
                     end do
                  end do
               end do
            end if
         end if
c
c     Do the wall-parallel Fourier transform and put into main storage
c
         if (.not.(os.and..not.physos)) then
            do i=1,3+scalar
c     Transform in x
               call vrfftf(u2(1,1,i),u2(2,1,i),
     &              w(1,1),w(2,1),nx,nz,2,nx+2,prex)
c     Transform in z
               call vcfftf(u2(1,1,i),u2(2,1,i),
     &              w(1,1),w(2,1),nz,nx/2,nx+2,2,prez)
c
c     Normalize
c
               do z=1,nz
                  do x=1,nx
                     u2(x,z,i)=u2(x,z,i)*(1./real(nx)/real(nz))
                  end do
               end do
               call putxzp(u2(1,1,i),y,i,ur)
            end do
         end if
      end do
c
c     Spectral space mode
c
      if (specm) then
         write(*,*)'Specm initial condition in (alfa,beta)'

         do i_spec=1,n_specm
            write(*,*)alfa(nalfa_s(i_spec)+1)*dstar,
     &           beta(nbeta_s(i_spec)+1)*dstar
c$$$  write(*,*)alfa(nalfa+1)*dstar,beta(nz-nbeta+1)*dstar
            do y=1,nyp
               y1=etab(y)
               if (fltype.eq.1.or.fltype.eq.2.or.fltype.eq.4.
     &              or.fltype.eq.5) then
                  gauss=(1-y1**2)**2
                  dgauss=-4*y1*(1-y1**2)
               else
                  gauss=(y1**2)
                  dgauss=2*y1
               end if
               gaussy=gauss*(1./(sqrt(2.*pi)*yscale1_s(i_spec)))*
     &              exp(-(y1-y0_s(i_spec))**2/(2.*
     &              (yscale1_s(i_spec))**2))
              
               Dgaussy=-gauss*((y1-y0_s(i_spec))/(sqrt(2.*pi)*
     &              yscale1_s(i_spec)**3))*
     &              exp(-(y1-y0_s(i_spec))**2/(2.*
     &              (yscale1_s(i_spec))**2))+
     &              dgauss*(1./(sqrt(2.*pi)*yscale1_s(i_spec)))*
     &              exp(-(y1-y0_s(i_spec))**2/(2.*
     &              (yscale1_s(i_spec))**2))

               ur(nalfa_s(i_spec)+1,y,nbeta_s(i_spec)+1,1)=
     &              ur(nalfa_s(i)+1,y,nbeta_s(i_spec)+1,1)
     &              +amp_s(i_spec)*Dgaussy*cmplx(0.,1.)*
     &           alfa(nalfa_s(i_spec)+1)*dstar

               ur(nalfa_s(i_spec)+1,y,nbeta_s(i_spec)+1,2)=
     &              ur(nalfa_s(i_spec)+1,y,nbeta_s(i_spec)+1,2)
     &              +amp_s(i_spec)*(alfa(nalfa_s(i_spec)+1)**2
     &              +beta(nbeta_s(i_spec)+1)**2)*
     &              gaussy*dstar**2

               ur(nalfa_s(i_spec)+1,y,nbeta_s(i_spec)+1,3)=
     &              ur(nalfa_s(i_spec)+1,y,nbeta_s(i_spec)+1,3)
     &              +amp_s(i_spec)*Dgaussy*cmplx(0.,1.)*
     &              beta(nbeta_s(i_spec)+1)*dstar

c$$$            nbeta=nz-nbeta
c$$$
c$$$            ur(nalfa+1,y,nbeta+1,1)=ur(nalfa+1,y,nbeta+1,1)
c$$$     &           +amp*Dgaussy*cmplx(0.,1.)*alfa(nalfa+1)*dstar
c$$$
c$$$            ur(nalfa+1,y,nbeta+1,2)=ur(nalfa+1,y,nbeta+1,2)
c$$$     &           +amp*(alfa(nalfa+1)**2+beta(nbeta+1)**2)*
c$$$     &           gaussy*dstar**2
c$$$
c$$$            ur(nalfa+1,y,nbeta+1,3)=ur(nalfa+1,y,nbeta+1,3)
c$$$     &           +amp*Dgaussy*cmplx(0.,1.)*beta(nbeta+1)*dstar
            end do
         end do
      end if
c
c     Reading perturbations from file
c
!     Both the pertfromfile and nois are added directly to the strengths of 
!     the x and z fourier transformed variables. 

      if (pertfromfile) then
!     Note that the phase difference has to be in terms of number of 2*pi. 
!     See definition below
         scal_pos = 8 + pressure + (scalar -1)
         do i=1,n_pertfromfile
            open(unit=17,file=initcond(i))
            do y=1,nyp
               read (17,*) u_real(y),u_imag(y),v_real(y),v_imag(y),
     &              w_real(y),w_imag(y),temp_real(y),temp_imag(y)
               init_cond(y,1) = cmplx(u_real(y),u_imag(y))
     &              *exp(im*phi_f(i)*2*pi)
               init_cond(y,2) = cmplx(v_real(y),v_imag(y))
     &              *exp(im*phi_f(i)*2*pi)
               init_cond(y,3) = cmplx(w_real(y),w_imag(y))
     &              *exp(im*phi_f(i)*2*pi)
               init_cond(y,4) = cmplx(temp_real(y),temp_imag(y))
     &              *exp(im*phi_f(i)*2*pi)
            end do
            do y=1,nyp
               ur(nalfa_f(i)+1,y,nbeta_f(i)+1,1)=
     &              ur(nalfa_f(i)+1,y,nbeta_f(i)+1,1)
     &              +amp_f(i)*init_cond(y,1)
               ur(nalfa_f(i)+1,y,nbeta_f(i)+1,2)=
     &              ur(nalfa_f(i)+1,y,nbeta_f(i)+1,2)
     &              +amp_f(i)*init_cond(y,2)
               ur(nalfa_f(i)+1,y,nbeta_f(i)+1,3)=
     &              ur(nalfa_f(i)+1,y,nbeta_f(i)+1,3)
     &              +amp_f(i)*init_cond(y,3)
               if(scalar.eq.1) then
                  ur(nalfa_f(i)+1,y,nbeta_f(i)+1,4)=
     &              ur(nalfa_f(i)+1,y,nbeta_f(i)+1,4)
     &              +amp_f(i)*init_cond(y,4)             
               end if
           
               if(nalfa_f(i).eq.0) then
                  
                  ur(nalfa_f(i)+1,y,nz-nbeta_f(i)+1,1)=
     &                 ur(nalfa_f(i)+1,y,nz-nbeta_f(i)+1,1)
     &                 +amp_f(i)*conjg(init_cond(y,1))
                  ur(nalfa_f(i)+1,y,nz-nbeta_f(i)+1,2)=
     &                 ur(nalfa_f(i)+1,y,nz-nbeta_f(i)+1,2)
     &                 +amp_f(i)*conjg(init_cond(y,2))
                  ur(nalfa_f(i)+1,y,nz-nbeta_f(i)+1,3)=
     &                 ur(nalfa_f(i)+1,y,nz-nbeta_f(i)+1,3)
     &                 +amp_f(i)*conjg(init_cond(y,3))
                  if(scalar.eq.1) then
                     ur(nalfa_f(i)+1,y,nz-nbeta_f(i)+1,4)=
     &                 ur(nalfa_f(i)+1,y,nz-nbeta_f(i)+1,4)
     &                 +amp_f(i)*conjg(init_cond(y,4))
                  end if
               end if
            end do
         end do
      end if
c
c     If the noise flag is true, noise in the form of random Stokes'
c     modes is added to the output
c     the noise has a total energy density given by ed
c     this is distributed evenly among the nxn*nyn*nzn*2 modes
c     note that nyn should not exceed approx 2/3*ny to allow the
c     Stokes' modes to be resolved
c     nzn should be odd to give the same number of modes for positive and
c     negative beta
c
      if (noise) call stnois(ur,seed,nxn,nyn,nzn,ed,
     &     alfa,beta,eta,uxy)
c
c     Write to file namnut
c
      write(*,'(a,i3,a,a)') ' Writing output file with',3+scalar,
     &     ' components : ',trim(namnut)
      do i=1,scalar
         write(*,*) '  scalar ',i,':',pr(i),m1(i),gr(i)
      end do
      write(*,*)'filename: ',namnut
      call wdiscpbl(ur,re,pr,gr,xl,zl,t,xs,dstar,fltype,
     &    bstart,blength,rlam,m1,spanv,namnut,3+scalar,pxy,urx,ivarvisc)

      end program bls
