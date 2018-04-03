c************************************************************************
c
c $HeadURL:  $
c $LastChangedDate:  $
c $LastChangedBy: deusebio@MECH.KTH.SE $
c $LastChangedRevision: 1449 $
c
c ***********************************************************************

      subroutine varvisc(taur,taui,ivarvisc,yb,alfa,beta,u0low,w0low,
     &     prexn,prezn,presn,prean,prex,prez,my_node,realg1,realg2,ur,
     &     ui,re,eta,it,scalr,scali,nur,nui,tau3r,tau3i,s2r,s2i,l2r,l2i,
     &     m2r,m2i,om2r,om2i,du2r,du2i,u2r,u2i,delT,Tl,Tref,xl,zl)

c     subroutine which add a term taur and taui which will act as a 
c     varying viscosity term. Be careful...this will restrict the CFL condition 
c     as well. Is it worth it? 


c     I saw that gur and gui are used in prhs.f and then on les.f so i should check 
c     if it will make any differences...


      implicit none

      include 'par.f'

      real pi
      parameter(pi = 3.1415926535897932385)


c
c     Declaration (Local)
c
      integer yb
      integer npl,nxy
      integer x,y,xy,z,i,ll,j,zz
      integer yy
c      real u2si(nx/2,nyp), u2sr(nx/2,nyp)

      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real wr(nxp/2+1,mby,nzd),wi(nxp/2+1,mby,nzd)

      real u0low,w0low 
      real eta(nyp)
      real u2r((nxp/2+1)*mby,nzd,3),u2i((nxp/2+1)*mby,nzd,3)
      real du2r((nxp/2+1)*mby,nzd,3,3),du2i((nxp/2+1)*mby,nzd,3,3)
      real om2r((nxp/2+1)*mby,nzd,3),om2i((nxp/2+1)*mby,nzd,3)
      integer it
c
c     Jose
c   
      real scalr((nxp/2+1)*mby,nzd,1),scali((nxp/2+1)*mby,nzd,1)
      real c1, temp, Tref, Tl, delT
      integer scalf,ivarvisc
      real gd2, re
      real nur((nxp/2+1)*mby,nzd),nui((nxp/2+1)*mby,nzd)
      real temp2, nu1, nu2

      real randsr((nxp/2+1)*mby,nzd,1),randsi((nxp/2+1)*mby,nzd,1)
      real randnur((nxp/2+1)*mby,nzd),randnui((nxp/2+1)*mby,nzd)
      real randnur_tmp((nxp/2+1)*mby,nzd)
      real randnui_tmp((nxp/2+1)*mby,nzd)

      real prexn(nx+15),prezn(nz*2+15),presn(nz+2+15),prean(nz*3/4+15)
      real prex(nxp+15),prez(nzp*2+15)
      real alfa(nx/2*mbz),beta(nz)

      real s2r   ((nxp/2+1)*mby,nzd,6)
      real s2i   ((nxp/2+1)*mby,nzd,6)
      real l2r   ((nxp/2+1)*mby,nzd,6)
      real l2i   ((nxp/2+1)*mby,nzd,6)
      real m2r   ((nxp/2+1)*mby,nzd,6)
      real m2i   ((nxp/2+1)*mby,nzd,6)
c      real test   ((nxp/2+1)*mby,nzd,6)

      real a(6),b(6)
      real sum_temp1, sum_temp2
      real xl,zl
      
      real taur (memnx,memny,memnz,6+3*scalar)
      real taui (memnx,memny,memnz,6+3*scalar)
      real tau3r((nxp/2+1),nzp,6) ,tau3i((nxp/2+1),nzp,6)
c
c     MPI
c
      integer my_node
      integer realg1,realg2 
      
      npl = min(mby,nyp-yb+1)
      nxy = (nxp/2+1)*npl-(nxp/2+1-nx/2)
      scalf = 8 + pressure + (scalar -1)

c
c     Initialising variables
c
      
      u2r = 0.
      u2i = 0.
      om2r = 0.
      om2i = 0.
      s2r = 0.
      s2i = 0.
      m2r = 0.
      m2i = 0.
      l2r = 0.
      l2i = 0.

c     --------------------------------------
c     computing Sij
c     --------------------------------------

c
c     Get the velocity and the vorticity omega_x and omega_z so we 
c     can avoid to use dcheb.f
c
c      write(*,*)'herr'
      if (nproc.eq.1) then
         do i=1,3
            call getxz(u2r(1,1,i),u2i(1,1,i),yb,i,0,ur,ui)
         end do
         call getxz(om2r(1,1,1),om2i(1,1,1),yb,4,0,ur,ui)
         call getxz(om2r(1,1,3),om2i(1,1,3),yb,5,0,ur,ui)

      else
#ifdef MPI
         do i=1,3
            call getpxz(u2r(1,1,i),u2i(1,1,i),yb,i,0,ur,ui,
     &           realg1,realg2,my_node)
         end do
         call getpxz(om2r(1,1,1),om2i(1,1,1),yb,4,0,ur,ui,
     &        realg1,realg2,my_node)
         call getpxz(om2r(1,1,3),om2i(1,1,3),yb,5,0,ur,ui,
     &        realg1,realg2,my_node)

#endif
      end if


c
c     Subtract lower wall shifting velocity
c

      do y=1,npl
         xy=1+(y-1)*(nxp/2+1)
         u2r(xy,1,1)=u2r(xy,1,1)-u0low
         u2r(xy,1,3)=u2r(xy,1,3)-w0low
      end do

     
c     
c     Calculate all velocity derivatives
c     horizontal derivatives by multiplication by alfa or beta
c
      do i=1,3
         do z=1,nz/2
            do y=1,npl
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)
                  du2r(xy,z,i,1)=-alfa(x)*u2i(xy,z,i)
                  du2i(xy,z,i,1)= alfa(x)*u2r(xy,z,i)
                  du2r(xy,z,i,3)=-beta(z)*u2i(xy,z,i)
                  du2i(xy,z,i,3)= beta(z)*u2r(xy,z,i)
               end do
            end do
         end do
      end do
c
c     The following loop accounts for dealiasing grid
c

      do i=1,3
         do z=nz/2+1,nz
            do y=1,npl
               do x=1,nx/2
                  zz = nzp - nz + z
                  xy=x+(y-1)*(nxp/2+1)
                  du2r(xy,zz,i,1)=-alfa(x)*u2i(xy,z,i)
                  du2i(xy,zz,i,1)= alfa(x)*u2r(xy,z,i)
                  du2r(xy,zz,i,3)=-beta(z)*u2i(xy,z,i)
                  du2i(xy,zz,i,3)= beta(z)*u2r(xy,z,i)
               end do
            end do
         end do
      end do

c     
c     The wall normal derivatives can be found from :
c     omx  = dwdy - dvdz
c     omz  = dvdx - dudy
c     i.e. :
c     dudy = dvdx - omz
c     dvdy =-dudx - dwdz   (from continuity)
c     dwdy = omx  + dvdz
c     
      do z=1,nz/2
         do y=1,npl
            do x=1,nx/2
               xy=x+(y-1)*(nxp/2+1)
               du2r(xy,z,1,2)=-alfa(x)*u2i(xy,z,2) - om2r(xy,z,3)
               du2i(xy,z,1,2)= alfa(x)*u2r(xy,z,2) - om2i(xy,z,3)
               du2r(xy,z,2,2)= alfa(x)*u2i(xy,z,1) +
     &              beta(z)*u2i(xy,z,3)
               du2i(xy,z,2,2)=-alfa(x)*u2r(xy,z,1) -
     &              beta(z)*u2r(xy,z,3)
               du2r(xy,z,3,2)= om2r(xy,z,1)        -
     &              beta(z)*u2i(xy,z,2)
               du2i(xy,z,3,2)= om2i(xy,z,1)        +
     &              beta(z)*u2r(xy,z,2)
            end do
         end do
      end do

      do z=nz/2+1, nz
         do y=1,npl
            do x=1,nx/2
               xy=x+(y-1)*(nxp/2+1)
               zz = nzp - nz + z
               du2r(xy,zz,1,2)=-alfa(x)*u2i(xy,z,2) - om2r(xy,z,3)
               du2i(xy,zz,1,2)= alfa(x)*u2r(xy,z,2) - om2i(xy,z,3)
               du2r(xy,zz,2,2)= alfa(x)*u2i(xy,z,1) +
     &              beta(z)*u2i(xy,z,3)
               du2i(xy,zz,2,2)=-alfa(x)*u2r(xy,z,1) -
     &              beta(z)*u2r(xy,z,3)
               du2r(xy,zz,3,2)= om2r(xy,z,1)        -
     &              beta(z)*u2i(xy,z,2)
               du2i(xy,zz,3,2)= om2i(xy,z,1)        +
     &              beta(z)*u2r(xy,z,2)
            end do
         end do
      end do

c     
c     Transform to physical space
c     
      do i=1,3
c     
c     No shifting for LES quantities
c     
         if (nfzsym.eq.0) then
            do j=1,3
               call vcfftb(du2r(1,1,i,j),du2i(1,1,i,j),wr,wi,nzp,
     &              nxy,(nxp/2+1)*mby,1,prez)
	    end do
         else
c     Understand this later. Right now not relevant as nfzsym = 0.
c
c
            if (i.le.2) then
               call vcffts(u2r(1,1,i),u2i(1,1,i),wr,
     &              nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
               call vcftab(om2r(1,2,i),om2i(1,2,i),
     &              wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
               call vcffts(du2r(1,1,i,1),du2i(1,1,i,1),wr,
     &              nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
               call vcffts(du2r(1,1,i,2),du2i(1,1,i,2),wr,
     &              nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
               call vcftab(du2r(1,2,i,3),du2i(1,2,i,3),
     &              wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
            else
               call vcftab(u2r(1,2,i),u2i(1,2,i),
     &              wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
               call vcffts(om2r(1,1,i),om2i(1,1,i),wr,
     &              nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
               call vcftab(du2r(1,2,i,1),du2i(1,2,i,1),
     &              wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
               call vcftab(du2r(1,2,i,2),du2i(1,2,i,2),
     &              wr,wi,nz/2-1,nxy,(nxp/2+1)*mby,1,prean)
               call vcffts(du2r(1,1,i,3),du2i(1,1,i,3),wr,
     &              nz/2+1,nxy,(nxp/2+1)*mby,1,presn)
            end if
         end if
         do j=1,3
            call vrfftb(du2r(1,1,i,j),du2i(1,1,i,j),wr,wi,
     &           nxp,nzpc*mby,1,nxp/2+1,prex)
         end do
      end do

c     
c     The velocity derivatives are in physical space 
c     

c     
c     Computation of the strain rate tensor Sij
c     

      do y=1,npl
         do z=1,nzp
            do x=1,nxp/2
               xy=x+(y-1)*(nxp/2+1)
               
               s2r(xy,z,1) =     du2r(xy,z,1,1)
               s2r(xy,z,2) = .5*(du2r(xy,z,2,1)+du2r(xy,z,1,2))
               s2r(xy,z,3) = .5*(du2r(xy,z,3,1)+du2r(xy,z,1,3))
               s2r(xy,z,4) =     du2r(xy,z,2,2)
               s2r(xy,z,5) = .5*(du2r(xy,z,2,3)+du2r(xy,z,3,2))
               s2r(xy,z,6) =     du2r(xy,z,3,3)
               
               s2i(xy,z,1) =     du2i(xy,z,1,1)
               s2i(xy,z,2) = .5*(du2i(xy,z,2,1)+du2i(xy,z,1,2))
               s2i(xy,z,3) = .5*(du2i(xy,z,3,1)+du2i(xy,z,1,3))
               s2i(xy,z,4) =     du2i(xy,z,2,2)
               s2i(xy,z,5) = .5*(du2i(xy,z,2,3)+du2i(xy,z,3,2))
               s2i(xy,z,6) =     du2i(xy,z,3,3)

            end do
         end do
      end do
      
c     
c     Computation of the viscosity
c     

      scalr = 0.
      scali = 0.
c      write(*,*) 'scalf:   ', scalf

c
c     For the scalar, first on N physical points. 
c
      if (nproc.eq.1) then
         call getxz(scalr(1,1,1),scali(1,1,1),yb,scalf,0,ur,ui)
      else
#ifdef MPI
         call getpxz(scalr(1,1,1),scali(1,1,1),yb,scalf,0,ur,ui,
     &        realg1,realg2,my_node)
#endif
      end if

      call vcfftb(scalr,scali,wr,wi,nz,nx/2,nxp/2+1,1,prezn)
      call vrfftb(scalr,scali,wr,wi,nx,nz,1,nxp/2+1,prexn)
c      do z = 1, nzp
c         do x = 1, nxp/2 + 1
c            write(20,*) x, z, scalr(x,z,1), scali(x,z,1)
c         end do
c      end do      

c
c     We have normalised scalar field now in physical space.
c

      randsr = scalr
      randsi = scali
      nur = 0.
      nui = 0.
      randnur = 0.
      randnui = 0.

      if (ivarvisc.ne.0) then            
c     
c     Jose: Viscosity field
c     Viscosity model - Arrhenius model. See Sameen, Govindarajan - 2007 JFM.
c     Tref = For defining reference viscosity. Tl = Lower wall temperature
c     delT = T_up - T_low
c     if (delT.lt.0.) then
c     Unstable stratification. 
c     Tref = 295 - delT
c     Tl = Tref
c     else
c     Stable stratification
c     Tref = 295 + delT
c     Tl = 295
c     end if
c     Remember only deviations of stress from that in uniform viscosity is required         
         c1 = exp(-ct/Tref)
         do z=1,nz
            do x=1,nx/2
		temp = delT*randsr(x,z,1) + Tl
c		temp2 = delT*scalr(x,z,1) + Tl

            	randnur(x,z) = c1*exp(ct/temp) - 1.0
c                write(21,*) x, z, temp, temp2

		temp = delT*randsi(x,z,1) + Tl
c		temp2 = delT*scali(x,z,1) + Tl
c                write(22,*) x, z, temp, temp2

                randnui(x,z) = c1*exp(ct/temp) - 1.0
	    enddo
         end do

      end if

      randnur_tmp = randnur
      randnui_tmp = randnui

c      do z = 1,1
c         do x = 1, nx/2
c            yy = 2*x - 1
c            gd2 = xl/real(nx)*real(yy-1)
c            write(22,*) gd2, randnur(x,z)
c            yy = 2*x
c            gd2 = xl/real(nx)*real(yy-1)
c            write(22,*) gd2, randnui(x,z)            
c         end do
c      end do
c
c     First forward Fourier transform to get nur, nui in Fourier space over N points.
c
      call vrfftf(randnur(1,1),randnui(1,1),wr,wi,
     &        nx,nz,1,nxp/2+1,prexn) 
      call vcfftf(randnur(1,1),randnui(1,1),wr,wi,nz,
     &        nxy,(nxp/2+1)*mby,1,prezn)

c
c     Normalise the Fourier coefficients.
c     If the normal backward Fourier transform were to be done now,
c     a factor of nx*nz will be multiplied in the result. 
c     Hence below
c
      randnur = randnur/real(nx*nz)
      randnui = randnui/real(nx*nz)

c
c     Remove the oddball in x
c
      x = nx/2 + 1
      do z = 1,nz
         randnur(x,z) = 0.
      end do
c
c     Remove the oddball in z
c
      z = nz/2 + 1
      do x = 1,nx/2
         randnur(x,z) = 0.
         randnui(x,z) = 0.
      end do

c      do z = 1, nzp
c         do x = 1, nxp/2 + 1
c            write(23,*) x, z, randnur(x,z), randnui(x,z)
c         end do
c     end do

c
c     Pad to get obtain it on the dealiased grid.
c     In z

      do z = nz/2+1,nz
         zz = nzp - nz + z
         do x = 1,nx/2
            randnur(x,zz) = randnur(x,z)
            randnur(x,z) = 0.
            randnui(x,zz) = randnui(x,z)
            randnui(x,z) = 0.
         end do
      end do

c
c     Backward Fourier transform onto to 3*N/2 points
c
      call vcfftb(randnur,randnui,wr,wi,nzp,nx/2,nxp/2+1,1,prez)
      call vrfftb(randnur,randnui,wr,wi,nxp,nzp,1,nxp/2+1,prex)

c      do z = 1,1
c         do x = 1, nxp/2
c            yy = 2*x - 1
c            gd2 = xl/real(nxp)*real(yy-1)
c            write(26,*) gd2, randnur(x,z)
c            yy = 2*x
c            gd2 = xl/real(nxp)*real(yy-1)
c            write(26,*) gd2, randnui(x,z)
c         end do
c      end do
      
c     
c     Calculation of the taur and taui terms
c     

      do z=1,nzp
         do x=1,nxp/2
            do i=1,6
               m2r(x,z,i)=-2.*randnur(x,z)*s2r(x,z,i)
               m2i(x,z,i)=-2.*randnui(x,z)*s2i(x,z,i)
            end do
         end do
      end do
c     
c     We have the stresses in physical space, let's go back to the 
c     spectral space
c     
      do i=1,6
         call vrfftf(m2r(1,1,i),m2i(1,1,i),wr,wi,
     &        nxp,nzpc*mby,1,nxp/2+1,prex)
         call vcfftf(m2r(1,1,i),m2i(1,1,i),wr,wi,nzp,
     &        nxy,(nxp/2+1)*mby,1,prez)
      end do 
c
c     Check here what the outputs are like. What are the frequencies 
c     that come out.
c
c      do i = 1,6
c         yy = 23 + i - 1
c         do z = 1, nzp
c            do x = 1, nxp/2 + 1
c               write(yy,*) x, z, m2r(x,z,2), m2i(x,z,2)
c            end do
c         end do
c      end do
c
c     It is seen that several frequencies are excited for wavenumbers
c     greater than nx/2, nz. All this must be discarded. They are however
c     fairly low.
c

c
c     Zero oddballs (The original frequencies)
c

      do i=1,6
         if(nz.gt.1) then
            z=nz/2+1
            do x=1,nx/2
               xy=x+(y-1)*(nxp/2+1)
               m2r(xy,z,i) = 0.
               m2i(xy,z,i) = 0.
            end do
            
            z = nzp - nz + nz/2 + 1
            do x=1,nx/2
               xy=x+(y-1)*(nxp/2+1)
               m2r(xy,z,i) = 0.
               m2i(xy,z,i) = 0.
            end do
         end if
         x=nx/2+1
         xy=x+(y-1)*(nxp/2+1)
         do z=1,nz/2
            m2r(xy,z,i) = 0.
            m2i(xy,z,i) = 0.
         end do
         do z=nz/2+1,nz
            zz = nzp + z - nz
            m2r(xy,zz,i) = 0.
            m2i(xy,zz,i) = 0.
         end do
      end do


c
c     Blow up array (could be avoided by changing putxz) and normalize
c
c     m is the one without padding while l has got 0s in the middle
c
      l2r=0.
      l2i=0.
      gd2=1./real(nxp*nzp)

      do ll=1,6
         do y=1,npl
            do z=1,nz/2
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)
                  l2r(xy,z,ll) = m2r(xy,z,ll)*gd2/re
                  l2i(xy,z,ll) = m2i(xy,z,ll)*gd2/re
               end do
            end do
            do z=nz/2+1,nz
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)
                  zz = z+nzp-nz
                  l2r(xy,zz,ll) = m2r(xy,zz,ll)*gd2/re
                  l2i(xy,zz,ll) = m2i(xy,zz,ll)*gd2/re
               end do
            end do
         end do
      end do

c      do i = 1,6
c         yy = 23 + i - 1
c         do z = 1, nzp
c            do x = 1, nxp/2 + 1
c               write(yy,*) x, z, l2r(x,z,2), l2i(x,z,2)
c            end do
c         end do
c      end do

c
c     Finally let's add this contribution to the previous tau
c

c
c     get, add and store back again
c
c      if(iles.eq.0) then
      do ll=1,6
         do y=1,npl
            yy= yb+(y-1)
c
c     ATTENTION: The 'getxz' procedure should be done in parallel
c
            if (nproc.eq.1) then
               call getxz(tau3r(1,1,ll),tau3i(1,1,ll),yy,ll,1,
     &              taur,taui)
            else
#ifdef MPI
               call getpxz(tau3r(1,1,ll),tau3i(1,1,ll),yy,ll,1,
     &              taur,taui,realg1,realg2,my_mode)
#endif
            end if

            do z=1,nz/2
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)
                  l2r(xy,z,ll) = l2r(xy,z,ll) + tau3r(xy,z,ll)
                  l2i(xy,z,ll) = l2i(xy,z,ll) + tau3i(xy,z,ll)
               end do
            end do
            do z=nz/2+1,nz
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)
                  zz = z+nzp-nz
                  l2r(xy,zz,ll) = l2r(xy,zz,ll)+tau3r(xy,zz,ll)
                  l2i(xy,zz,ll) = l2i(xy,zz,ll)+tau3i(xy,zz,ll)
               end do
            end do
         end do
      end do

c      do i = 1,6
c         yy = 30 + i - 1
c         do z = 1, nzp
c            do x = 1, nxp/2 + 1
c               write(yy,*) x, z, l2r(x,z,2), l2i(x,z,2)
c            end do
c         end do
c      end do
c     end if

c      write(12,314) m2r(1,1,2)

c
c     Store back again
c

c
c     Save the output
c
           
c      do z=1,nzp
c         do x=1,nxp/2+1
c            do i=1,6
c               a=l2r(x,z,2)
c               b=l2i(x,z,2)
c            end do
c            write(44,'(12F)') a,b
c
c        end do
c     end do
c     
c     Put tau_ij onto global field
c
314   format(f12.6)
      if (nproc.eq.1) then
         do i=1,6
            call putxz(l2r(1,1,i),l2i(1,1,i),yb,i,taur,taui)
         end do
      else
#ifdef MPI
         do i=1,6
            call putpxz(l2r(1,1,i),l2i(1,1,i),yb,i,taur,taui,
     &           realg1,realg2,my_node)
         end do
#endif
      end if

      end subroutine varvisc
