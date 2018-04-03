c************************************************************************
c
c $HeadURL:  $
c $LastChangedDate:  $
c $LastChangedBy: deusebio@MECH.KTH.SE $
c $LastChangedRevision: 1449 $
c
c ***********************************************************************

      subroutine varvisc2(taurvis,tauivis,ivarvisc,yb,alfa,beta,u0low,
     &     w0low,prexn,prezn,presn,prean,my_node,realg1,
     &     realg2,ur,ui,re,eta,cr,ci,nur,nui,tauplr,taupli,s2r,s2i,
     &     l2r,l2i,m2r,m2i,om2r,om2i,du2r,du2i,u2r,u2i,delT,Tl,Tref,
     &     xl,zl,alp_co,beta_co,alp_fl,alp_wn,beta_fl,beta_wn,beta_loc)

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

      integer alp_co,beta_co
      integer alp_wn(nx/2),beta_wn(nz)
      integer beta_loc(2*nz/3 + 1)
      real alp_fl(nx/2),beta_fl(nz)
c      real u2si(nx/2,nyp), u2sr(nx/2,nyp)

      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real wr(nxp/2+1,mby,nzd),wi(nxp/2+1,mby,nzd)

      real u0low,w0low 
      real eta(nyp)
      real u2r((nxp/2+1)*mby,nzd,3),u2i((nxp/2+1)*mby,nzd,3)
      real du2r((nxp/2+1)*mby,nzd,3,3),du2i((nxp/2+1)*mby,nzd,3,3)
      real om2r((nxp/2+1)*mby,nzd,3),om2i((nxp/2+1)*mby,nzd,3)

c
c     Jose
c
      real taurvis(memnx,memny,memnz,6), tauivis(memnx,memny,memnz,6)   
      real cr((nxp/2+1)*mby,nzd),ci((nxp/2+1)*mby,nzd)
      real c1, temp, Tref, Tl, delT
      integer scalf,ivarvisc
      real gd2, re
      real nur((nxp/2+1)*mby,nzd),nui((nxp/2+1)*mby,nzd)
c      real temp2, nu1, nu2

      real randsr((nxp/2+1)*mby,nzd),randsi((nxp/2+1)*mby,nzd)
c      real randnur((nxp/2+1)*mby,nzd),randnui((nxp/2+1)*mby,nzd)
c     real randnur_tmp((nxp/2+1)*mby,nzd)
c      real randnui_tmp((nxp/2+1)*mby,nzd)

      real prexn(nx+15),prezn(nz*2+15),presn(nz+2+15),prean(nz*3/4+15)
      real alfa(nx/2*mbz),beta(nz)
      real xl,zl

      real s2r   ((nxp/2+1)*mby,nzd,6)
      real s2i   ((nxp/2+1)*mby,nzd,6)
      real l2r   ((nxp/2+1)*mby,nzd,6)
      real l2i   ((nxp/2+1)*mby,nzd,6)
      real m2r   ((nxp/2+1)*mby,nzd,6)
      real m2i   ((nxp/2+1)*mby,nzd,6)
c      real test   ((nxp/2+1)*mby,nzd,6)

c      real a(6),b(6)
c      real sum_temp1, sum_temp2

      
c      real taur (memnx,memny,memnz,memnxyz)
c      real taui (memnx,memny,memnz,memnxyz)
      real tauplr((nxp/2+1),nzp) ,taupli((nxp/2+1),nzp)
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

      cr = 0.
      ci = 0.
      nur = 0.
      nui = 0.
      s2r = 0.
      s2i = 0.
      m2r = 0.
      m2i = 0.
      l2r = 0.
      l2i = 0.
      u2r = 0.
      u2i = 0.
      om2r = 0.
      om2i = 0.
      du2r = 0.
      du2i = 0.
      
c     --------------------------------------
c     computing Sij
c     --------------------------------------
c      if (yb.eq.nyp) then
c         write(*,*)'Here start.my node = ', my_node
c         write(*,*)'yb =  ',yb
c      end if
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
         call getxz(cr(1,1),ci(1,1),yb,scalf,0,ur,ui)
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
         call getpxz(cr(1,1),ci(1,1),yb,scalf,0,
     &        ur,ui,realg1,realg2,my_node)         
#endif
      end if

c
c     Random check
c

c      if(yb.eq.15) then
c         do 
c      endif
      


c
c     Zero oddball: Done in getpxz.f
c

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
         do z=1,nz
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
c     The wall normal derivatives can be found from :
c     omx  = dwdy - dvdz
c     omz  = dvdx - dudy
c     i.e. :
c     dudy = dvdx - omz
c     dvdy =-dudx - dwdz   (from continuity)
c     dwdy = omx  + dvdz
c     
      do z=1,nz
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



      do i=1,3
c     
c     No shifting for LES quantities
c     
         if (nfzsym.eq.0) then
            do j=1,3
               call vcfftb(du2r(1,1,i,j),du2i(1,1,i,j),wr,wi,nz,
     &              nxy,(nxp/2+1)*mby,1,prezn)
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
     &           nx,nzc*mby,1,nxp/2+1,prexn)
         end do
      end do

c     
c     The velocity derivatives are in physical space in nx x nz grid.
c     

c     
c     Computation of the strain rate tensor Sij
c     

      do y=1,npl
         do z=1,nz
            do x=1,nx/2
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
c     Bringing back scalar to real space
      call vcfftb(cr,ci,wr,wi,nz,nx/2,nxp/2+1,1,prezn)
      call vrfftb(cr,ci,wr,wi,nx,nz,1,nxp/2+1,prexn)

c
c     We have normalised scalar field now in physical space.
c

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
            temp = delT*cr(x,z) + Tl
            nur(x,z) = c1*exp(ct/temp) - 1.0

            temp = delT*ci(x,z) + Tl
            nui(x,z) = c1*exp(ct/temp) - 1.0

         enddo
      end do

c
c     First forward Fourier transform to get nur, nui in Fourier space over N points.
c
      call vrfftf(nur(1,1),nui(1,1),wr,wi,
     &        nx,nz,1,nxp/2+1,prexn) 
      call vcfftf(nur(1,1),nui(1,1),wr,wi,nz,
     &        nxy,(nxp/2+1)*mby,1,prezn)

c
c     Normalise the Fourier coefficients.
c     If the normal backward Fourier transform were to be done now,
c     a factor of nx*nz will be multiplied in the result. 
c     Hence below
c
      nur = nur/real(nx*nz)
      nui = nui/real(nx*nz)


c
c     Jose: Dealiased Fourier coefficients of stress.
c     Without moving to physical space and back. 
c     Attempt to save time per iteration.

      do i = 1,6
         do z = 1,nz/2
            wn_b = beta_wn(z)
            zz = nz + 2 - z 
            do x = 1,nx/2
               wn_a = alp_wn(x)
               sumr_p = 0.
               sumi_p = 0.
               sumr_n = 0.
               sumi_n = 0.
               loc_a1 = 1
               do tmpa = 1,(wn_a  + 1)
                  loc_a2 = wn_a - loc_a1 + 2
                  fact_a = alp_fl(loc_a1)*alp_fl(loc_a2)
                  if (fact_a.ne.0) then
                     bt_st = wn_b - beta_co
                     do tmpb1 = bt_st,beta_co,1
                        if(tmpb1.ge.0) then
                           loc_b1 = tmpb1  + 1
                        elseif 
                           loc_b1 = nz + 1 + tmpb1
                        endif
c     loc_b1 = nz + 1 + tmpb1
                        tmpb2 = wn_b - tmpb1
                        
c     if(abs(tmpb2).le.beta_co) then
c     For the positive set of spanwise wavenumbers
c     
                        if(tmpb2.ge.0) then
                           loc_b2 = tmpb2  + 1
                        elseif 
                           loc_b2 = nz + 1 + tmpb2
                        endif
c     factb = 
c     if(tmpb2.eq.tmpb1) then
                        tr_p = nur(loc_a1,loc_b1)*s2r(loc_a2,loc_b2) -
     &                       nui(loc_a1,loc_b1)*s2i(loc_a2,loc_b2)
                        ti_p = nur(loc_a1,loc_b1)*s2r(loc_a2,loc_b2) +
     &                       nui(loc_a1,loc_b1)*s2i(loc_a2,loc_b2)
c     else
c     tr = nur(loc_a1,loc_b1)*s2r(loc_a2,loc_b2) -
c     &                       nui(loc_a1,loc_b1)*s2i(loc_a2,loc_b2) +
c     &                       nur(loc_a1,loc_b2)*s2r(loc_a2,loc_b1) -
c     &                       nui(loc_a1,loc_b2)*s2i(loc_a2,loc_b1)
                        
c     ti = nur(loc_a1,loc_b1)*s2r(loc_a2,loc_b2) +
c     &                       nui(loc_a1,loc_b1)*s2i(loc_a2,loc_b2) +
c     &                       nur(loc_a1,loc_b2)*s2r(loc_a2,loc_b1) +
c     &                       nui(loc_a1,loc_b2)*s2i(loc_a2,loc_b1)
c     endif
                        
                        sumr_p = sumr_p + 2*tr_p
                        sumi_p = sumi_p + 2*ti_p
c     endif
                        if(wn_b.ne.0) then
                           tmpb1a = -1*tmpb1
                           tmpb2a = -1*tmpb2
                           if(tmpb1a.ge.0) then
                              loc_b1 = tmpb1  + 1
                           elseif 
                              loc_b1 = nz + 1 + tmpb1a
                           endif
                           if(tmpb2a.ge.0) then
                              loc_b2 = tmpb2  + 1
                           elseif 
                              loc_b2 = nz + 1 + tmpb2a
                           endif
                           tr_n=nur(loc_a1,loc_b1)*s2r(loc_a2,loc_b2,i)-
     &                          nui(loc_a1,loc_b1)*s2i(loc_a2,loc_b2,i)
                           ti_n=nur(loc_a1,loc_b1)*s2r(loc_a2,loc_b2,i)+
     &                          nui(loc_a1,loc_b1)*s2i(loc_a2,loc_b2,i)
                           sumr_n = sumr_n + 2*tr_n
                           sumi_n = sumi_n + 2*ti_n
                        endif
                     enddo
                  endif
                  
                  loc_a1 = loc_a1 + 1
               enddo
               m2r(x,z,i) = -1.*sumr_p
               m2i(x,z,i) = -1.*sumi_p
               if(wn_b.ne.0)then 
                  m2r(x,zz,i) = -1.*sumr_p
                  m2i(x,zz,i) = -1.*sumi_p
               endif
            enddo
         enddo
      enddo

c
c     The oddball wavenumber is never evaluated.
c     Hence it is taken care of.
c


      do z=1,nz
         do x=1,nx/2
            do i=1,6
               m2r(x,z,i)=-2.*nur(x,z)*s2r(x,z,i)
               m2i(x,z,i)=-2.*nui(x,z)*s2i(x,z,i)
            end do
         end do
      end do

c     
c     We have the stresses in physical space, let's go back to the 
c     spectral space
c     
      do i=1,6
         call vrfftf(m2r(1,1,i),m2i(1,1,i),wr,wi,
     &        nx,nzc*mby,1,nxp/2+1,prexn)
         call vcfftf(m2r(1,1,i),m2i(1,1,i),wr,wi,nz,
     &        nxy,(nxp/2+1)*mby,1,prezn)
      end do 
c
c     Zero oddballs
c     2/3 dealiasing as well.

      do i=1,6

c     Dealiasing and oddball
c     xy = starting index
         xy = nx/3 + 2
         do z=1,nz/2
            do x = xy,nx/2
               m2r(xy,z,i) = 0.
               m2i(xy,z,i) = 0.
            end do
         end do

         if(nz.gt.1) then
c     Oddball
            z=nz/2+1
            do x=1,nx/2
               xy=x+(y-1)*(nxp/2+1)
               m2r(xy,z,i) = 0.
               m2i(xy,z,i) = 0.
            end do 

c     Dealiasing
            zz = nz/3 + 2
c     zz = starting index
            do z = zz,nz/2
               do x = 1,nx/2
                  xy=x+(y-1)*(nxp/2+1)
                  m2r(xy,z,i) = 0.
                  m2i(xy,z,i) = 0.
                  j = nz + 2 - z;
                  m2r(xy,j,i) = 0.
                  m2i(xy,j,i) = 0.
               end do
            end do
         end if

      end do

c
c     Blow up array (could be avoided by changing putxz) and normalize
c
c     m is the one without padding while l has got 0s in the middle
c
      gd2=1./real(nx*nz)

      do ll=1,6
         do y=1,npl
            do z=1,nz
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)
                  l2r(xy,z,ll) = m2r(xy,z,ll)*gd2/re
                  l2i(xy,z,ll) = m2i(xy,z,ll)*gd2/re
               end do
            end do
         end do
      end do
           
c     
c     Put tau_ij onto global field
c
314   format(f12.6)

      if (nproc.eq.1) then
         do i=1,6
            call putxz(l2r(1,1,i),l2i(1,1,i),yb,i,taurvis,tauivis)
         end do
      else
#ifdef MPI
         do i=1,6
            call putpxz(l2r(1,1,i),l2i(1,1,i),yb,i,taurvis,tauivis,
     &           realg1,realg2,my_node)
         end do
#endif
      end if

      end subroutine varvisc2
