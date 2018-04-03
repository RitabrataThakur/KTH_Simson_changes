C************************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/bla/bla.f $
c $LastChangedDate: 2013-03-15 11:44:18 +0100 (Fri, 15 Mar 2013) $
c $LastChangedBy: pschlatt@MECH.KTH.SE $
c $LastChangedRevision: 1821 $
c
c ***********************************************************************
#ifdef WRAPPER
      subroutine bla(ur,ui,u0low,u0upp,w0low,w0upp,du0upp,re,pr,
     &     xlold,zl,t,xs,dstar,fltype,bstart,blength,rlam,m1,spanv,
     &     spat,varsiz,w3,urx,nxtmp,my_node,all_nodes,
     &     tstatus,tmin,alam,aturb)
#else
      program bla
#endif

      implicit none

#include "par.f"
#ifdef MPI
      include 'mpif.h'
#endif
c
c     Flags to turn on/off certain features
c
      integer,parameter :: iiles = 1
      integer,parameter :: iiles_smag = iiles*1

      integer,parameter :: iimhd = 0

      integer,parameter :: iisfd  = 0

      integer,parameter :: iiwall = 1
      integer,parameter :: iibase = 0
      integer,parameter :: iilen  = 0
      integer,parameter :: iiwave = 0
      integer,parameter :: iipert = 0
      integer,parameter :: iibf3  = 0

!     What to fill for these parameters? Old values 500, 30, 36 
      integer,parameter :: nsave = 100
      integer,parameter :: mpl   = 10
      integer,parameter :: nwave = 10

      integer,parameter :: arnoldi = 0

#ifdef WRAPPER
      integer,parameter :: wrapper = 1
#else
      integer,parameter :: wrapper = 0
#endif
c
c     Main storage (distributed in z among the processors)
c
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
c
c     In case of memory related problems, the following lines
c     can be tried out:
c
c      common /MAIN/ ur,ui
c
c      static ur,ui
c
c      real,allocatable :: ur(:,:,:,:),ui(:,:,:,:)
c      allocate(ur(memnx,memny,memnz,memnxyz)
c      allocate(ui(memnx,memny,memnz,memnxyz)
c
c
c     Step 2 storage (nonlinbl)
c

      real u2r ((nxp/2+1)*mby*nzd*3       ,nthread)
      real u2i ((nxp/2+1)*mby*nzd*3       ,nthread)
      real om2r((nxp/2+1)*mby*nzd*3       ,nthread)
      real om2i((nxp/2+1)*mby*nzd*3       ,nthread)
      real wr  ((nxp/2+1)*mby*nzd         ,nthread)
      real wi  ((nxp/2+1)*mby*nzd         ,nthread)
      real du2r((nxp/2+1)*mby,nzd,3,3     ,nthread)
      real du2i((nxp/2+1)*mby,nzd,3,3     ,nthread)
      real th2r((nxp/2+1)*mby*nzd*4*scalar,nthread)
      real th2i((nxp/2+1)*mby*nzd*4*scalar,nthread)
c
c     Particles storage
c
      real part     (10,nppart)
      real pos      (3,nppart)
      real pint     (3,npart,3)
      real taupp    (nppart)
      integer p,pseed
      real pdummy,pdummy1
c
c     Step 3 storage (linearbl,prhs)
c
      real u3r   (nx/2*mbz*nyp*3,nthread),u3i  (nx/2*mbz*nyp*3,nthread)
      real om3r  (nx/2*mbz*nyp*3,nthread),om3i (nx/2*mbz*nyp*3,nthread)
      real h3r   (nx/2*mbz*nyp*3,nthread),h3i  (nx/2*mbz*nyp*3,nthread)
      real domyr (nx/2*mbz*nyp  ,nthread),domyi(nx/2*mbz*nyp  ,nthread)
      real domyh (nx/2*mbz*nyp  ,nthread),omyh (nx/2*mbz*nyp  ,nthread)
      real homyr (nx/2*mbz*nyp  ,nthread),homyi(nx/2*mbz*nyp  ,nthread)
      real dvr   (nx/2*mbz*nyp  ,nthread),dvi  (nx/2*mbz*nyp  ,nthread)
      real dvh   (nx/2*mbz*nyp  ,nthread),dvh2 (nx/2*mbz*nyp  ,nthread)
      real vr    (nx/2*mbz*nyp  ,nthread),vi   (nx/2*mbz*nyp  ,nthread)
      real vh    (nx/2*mbz*nyp  ,nthread),vh2  (nx/2*mbz*nyp  ,nthread)
      real hvr   (nx/2*mbz*nyp  ,nthread),hvi  (nx/2*mbz*nyp  ,nthread)
      real w3    (nx/2*mbz*nyp  ,nthread)
      real d2omyh(nx/2*mbz*nyp  ,nthread)
      real th3r  (nx/2*mbz*nyp*scalar,nthread)
      real th3i  (nx/2*mbz*nyp*scalar,nthread)
      real hthr  (nx/2*mbz*nyp*scalar,nthread)
      real hthi  (nx/2*mbz*nyp*scalar,nthread)
      real pthr  (nx/2*mbz*nyp*scalar,nthread)
      real pthi  (nx/2*mbz*nyp*scalar,nthread)
      real dthr  (nx/2*mbz*nyp*scalar,nthread)
      real dthi  (nx/2*mbz*nyp*scalar,nthread)
      real c     (nx/2*mbz*ny    ,nthread)
      real d     (nx/2*mbz*ny    ,nthread)
      real e     (nx/2*mbz*ny    ,nthread)
      real q     (nx/2*mbz*(ny+2),nthread)
      real puw(ny,2+scalar)
c
c     Step 3 collection variables
c
      real f  (nx/2*mbz*nyp,5+2*scalar,nthread)
      real bis(nx/2*mbz*nyp,6+2*scalar,nthread)
      real d2v(nx/2*mbz*nyp,4         ,nthread)
c
c     Variables for tbc>0
c
      real fth_bc   (nx/2*mbz*nyp*2*scalar,nthread)
      real d2th_bc  (nx/2*mbz*nyp*4*scalar,nthread)
      real dth3r_bc (nx/2*mbz*nyp  *scalar,nthread)
      real dth3i_bc (nx/2*mbz*nyp  *scalar,nthread)
      real dth3r_bc2(nx/2*mbz*nyp  *scalar,nthread)
      real dth3i_bc2(nx/2*mbz*nyp  *scalar,nthread)
c
c     Partial right-hand sides (prhs)
c
      real pomyr (nx/2,mbz,nyp),pomyi (nx/2,mbz,nyp)
      real pvr   (nx/2,mbz,nyp),pvi   (nx/2,mbz,nyp)
c
c     Temporary storage
c
      real boxr  (nx/2,mbz,nyp),boxi  (nx/2,mbz,nyp)
c
c     Pressure
c
      real h2r((nxp/2+1)*mby*nzd*3*pressure,nthread)
      real h2i((nxp/2+1)*mby*nzd*3*pressure,nthread)
      real rv1r(nx/2,nzc,2*pressure)
      real rv1i(nx/2,nzc,2*pressure)
c
c     Time stepping parameters
c
      logical vart
      real t,tc,tmax,tleft,dtmax,dtn,dtnp1
      real an,bn,anp1,bnp1,cnp1,anrk(4,2),bnrk(4,2),cnrk(4,2)
      integer it,maxit,nst
      integer rksubstep
      integer stop_now
c
c     Numerics
c
      logical icorr,cim
c
c     File communication
c
      logical write_inter
      logical varsiz
      character*80 namnin,namnut,nmsave(nsave)
      character*80 namnut_inter,namnut_sfd,dumpname
      real tsave(nsave)
      integer isave
      real urx(nx)
c
c     Runtime plane saving
c
      real wxy(nx,nyp),wyz(nyp,nz)
      real plxy(nx,nyp/nproc+1)
      character*80 nampl(mpl)
      integer ipl,npl,tpl(mpl,3)
      real cpl(mpl)
      integer wplblodd(mpl)
c
c     Runtime statistics
c
      logical do_stat,do_press,do_press_nl
      logical longli,fileurms
      integer icfl,iamp,mwave,ixys,ixyss
      character*80 namamp,namwav,namext,namxys
      real cfl,cflmax,cflp((nyp/nproc+1)*nproc),sumw,txys
      integer kx(nwave),kz(nwave)
      real amp(nyp,20)
      complex campw(nyp,4,nwave)
      integer iext
      real vext(nyp,2,6),cext(nyp,2,6,2)
      real xys     (nx,nyp/nproc+1,nxys)
      real xysth   (nx,nyp/nproc+1,nxysth,scalar)
c
c     Flow
c
      character*80 namfre,nambfl
      logical spat,gall,tabfre,rbfl
      integer fltype,ibc
      real re,px,rot,cdev,pr(scalar),m1(scalar),retau,thgrad(scalar)
      real gr(scalar)
      real mflux
      logical cflux
      real u0low,u0upp,w0low,w0upp,du0upp
      real spanv,rlam
c
c     Disturbance volume force
c
      integer loctyp
      real fp1,fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7
      real fpds8,fpdds1,fpdds2,fpdds3,fpdds4,fpdds5,fpdds6
c
c     Fringe function
c
      real fmax,fstart,fend,frise,ffall
      real cphas1(osnf,nxp/2,nzd),cphas2(osnf,nxp/2,nzd)
      real sphas1(osnf,nxp/2,nzd),sphas2(osnf,nxp/2,nzd)
      real fring1(nxp/2),fring2(nxp/2)
      real xc1(nxp/2),xc2(nxp/2)
c
c     Trip forcing
c
      logical tripf
      integer seed,nzt,ntdt
      real tdt,ttzc(nzp+2,4),tamps,tampt,tysc,txsc,tx0
c
c     Blowing and suction
c
      integer wbci
      real wp1,wp2
      real wpds1,wpds2,wpds3,wpds4,wpds5
      real wpds6,wpds7,wpds8,wpds9,wpds10,wpds11,wpds12
      real wpdds1,wpdds2
      real wallvr(nxp/2+1,nzd*iiwall),wallvi(nxp/2+1,nzd*iiwall)
c
c     Boundary conditions for scalar fiels
c
      integer tbc(scalar)
      real dtheta0_upp(scalar),dtheta0_low(scalar)
      real theta0_low(scalar),theta0_upp(scalar)
c
c     Wall roughness (wbci==-1)
c
      integer nxm,nxn,nzn,lenx
      parameter(nxm=nxp/2+1,nxn=nx/2+1,nzn=nz/2+1)
      parameter (lenx=nx/8)
      integer every,nzr,rseed,hst,hen,hlen
      real delh,h_rough,hstart,hend,hrise,hfall,psi,tpsi,switch
      real wallur(nxm,memnz*iiwall),wallui(nxm,memnz*iiwall)
      real wlvr(nxm,memnz*iiwall),wlvi(nxm,memnz*iiwall)
      real wallwr(nxm,memnz*iiwall),wallwi(nxm,memnz*iiwall)
      real rzd(nz+2)
      character(32) wallroughinput,roughfile
      logical updat,v_wall,taylor4,zrand,rghfil,hrms
c
c     Jet in crossflow (wbci==-2)
c
      real jet_diam,xjet,zjet,jetmag
      integer jet_prof
c
c     Wave generation
c
      real evr(nyp),evi(nyp),eur(nyp),eui(nyp),ewr(nyp),ewi(nyp)
      real afw,alfaw,betaw,ampob,amp2d
      real ev2r(nyp),ev2i(nyp),eu2r(nyp),eu2i(nyp),afw2d,alf2d
c
c     Orr-Sommerfeld modes
c
      logical osmod,osdamp
      real osamp
      integer osn,osnumb
      real osymax,osbeta(osnf), osre, osomega(osnf)
      real osalr(osnf), osali(osnf)
      real osur(osnf,nyp), osvr(osnf,nyp), oswr(osnf,nyp)
      real osui(osnf,nyp), osvi(osnf,nyp), oswi(osnf,nyp)
      character*80 osfil
c
c     Asymptotic suction layer
c
      logical suction,asbl
      real vsuc
c
c     Streak generation
c
      logical streak
      character*80 str_nam1,str_nam2
      real iampst,ampst(2),tsmoo(4),tsmst(2),tsmend(2)
      real betast(2),omegast(2),phist
      real uust_r(nyp,180,2*iiwave),vvst_r(nyp,180,2*iiwave)
      real wwst_r(nyp,180,2*iiwave)
      real uust_i(nyp,180,2*iiwave),vvst_i(nyp,180,2*iiwave)
      real wwst_i(nyp,180,2*iiwave)
      integer ndxst
c
c     Wave in fringe generation parameters
c
      logical waves
      real waamp,wamoo,wamst,waiamp
      real omega
      real ucwa(nyp,180*iiwave),uswa(nyp,180*iiwave)
      real vcwa(nyp,180*iiwave),vswa(nyp,180*iiwave)
      integer ndxwa
c
c     Box extension
c
      real xlold
      integer nxtmp
c
c     Geometrics
c
      real alfa(nx/2,mbz),beta(nz)
      real eta(nyp),deta(nyp),wint(nyp),wd1(nyp,4),xl,zl
      real gridx(nx),gridy(nyp),gridz(nz)
      real deltaxyz2(nyp)
      real zs,zsc,xs,xsc,xs0
      real x0,dstar,dstar2,xblc,xbl
c
c     Base flow
c
      integer nbla
      real dybla
      real fbla(mbla,7+3*scalar)
      real wbr(nxp/2+1,nyp),wbi(nxp/2+1,nyp)
      real bstart,blength
c
c     Base flow and derivatives
c
      real bu1(nxp/2+1,nyp,3+scalar),bu2(nxp/2+1,nyp,3+scalar)
      real bu3(nxp/2+1,nyp*iipert,3+scalar)
      real bu4(nxp/2+1,nyp*iipert,3+scalar)
      real bu1jr(nxp/2+1,3+scalar,3),bu1ji(nxp/2+1,3+scalar,3)
c
c     Fixed base flow and derivatives for blshift, blfou
c
      real bu1f(nxp/2+1,nyp*iibase,3+scalar)
      real bu2f(nxp/2+1,nyp*iibase,3+scalar)
      real bu1jr0(nxp/2+1,3+scalar,3*iibase)
      real bu1ji0(nxp/2+1,3+scalar,3*iibase)
c
c     Base flow for lenbox
c
      real bucorr(nxp,nyp,3*iilen)
      real pxz(nxp+2,nzd*iilen),wxz(nxp+2,nzd)
c
c     FFT preprocessing data
c
      real prex(nxp+15),prey(nyp*2+15)
      real prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)
      real prexn(nx+15),prezn(nz*2+15),presn(nz+2+15),prean(nz*3/4+15)
      real prezr(nzp+15)
c
c     Two-point correlation
c
      character*80 corrnam,corrnam_x
      real corrx(mcorr), corry(mcorr)
      real corrz(mcorr), corry_x(mcorr)
      integer corrt(mcorr),corrt_x(mcorr)
      real corrxx(mcorr),corryy(mcorr)
      real corrzz(mcorr),corryy_x(mcorr)
      real corr(nzc+2,mcorr),corr_x(nx+2,mcorr)
      real totcorr(nzc+2),totcorr_x(nx+2)
      real corr_ms(4,mcorr),totcorr_ms(4),totcorr_ms_x(4)
      real corr_ms_x(4,mcorr)
      integer corrxi(mcorr),corrxodd(mcorr),corryi(mcorr)
      integer corrzi(mcorr),corryi_x(mcorr)
      logical corrf,corrf_x
      integer ncorr,ncorr_x
      real corrwsave(nzc+15)
      real corrwsave_x(nx+15)
c
c     Time series
c
      real series(msamp,0:mser),totseries(msamp,0:mser)
      character(len=80) namser
      real serc(3,mser)
      real sercc(3,mser)
      integer serci(4,mser)
      integer sert(mser)
      logical serf
      integer nser,sercount,nsamp
c
c     Constant
c
      real pi
      parameter (pi = 3.1415926535897932385)
c
c     Loop indices, pointers
c
      integer yb,zb,myb,i,j
c
c     Indices for parallelization
c
      integer ybp,zbp,ip
c
c     Timing
c
      character*20 sttimec,timdat
      real sttimet,delta_wall_time
      real,parameter :: delta_wall_time_int = 10.
      integer kode
      ! Start time
      real wtime0,ctime0
      ! Timings for one time step
      real wtime1,ctime1
      real wtime2,ctime2
      ! Partial timing
      real wtime3,ctime3
      real wtime4,ctime4
      ! Different timings
      real wtime_stat,ctime_stat
      real wtime_nonlin,ctime_nonlin
      real wtime_linear,ctime_linear
      real wtime_io,ctime_io
      real wtime_les,ctime_les
      real wtime_tau,ctime_tau
      real wtime_startup,ctime_startup
      real wtime_end,ctime_end
      real wtime_rest,ctime_rest
      real wtime_trip,ctime_trip
      real wtime_ffun,ctime_ffun
      real wtime_cfl,ctime_cfl
      real wtime_planes,ctime_planes
      real wtime_press,ctime_press
      real wtime_part,ctime_part
      real cpumax,wallmax
      real tcom,tser
c
c     MHD
c
      real fmhdr(iimhd*memnx*memny*memnz*2)
      real fmhdi(iimhd*memnx*memny*memnz*2)
      real j2r((nxp/2+1)*mby*nzd*4*iimhd,nthread)
      real j2i((nxp/2+1)*mby*nzd*4*iimhd,nthread)
      integer imhd
      real b0(3),mhd_n
c
c     Jose: Varying viscosity variables
c
      real delT, Tl, Tu, Tref
      integer ivarvisc
c
c     Jose: Unstable mode excitation check
c
      integer,parameter :: j_unch = 5
      logical j_un_chk
      integer n_un, j_un_cnt
      character*80 nam_un_chk, fil_nam
      character fil_no
      integer un_alp(j_unch),un_bet(j_unch)
      real vel_unchk_r(nyp,3,j_unch), vel_unchk_i(nyp,3,j_unch)
      real vel_unchk_pr(nyp,3,j_unch), vel_unchk_pi(nyp,3,j_unch)
      real vel_unchk_lr(nyp,3,j_unch), vel_unchk_li(nyp,3,j_unch)
      real en(j_unch)
      real proj1(4,j_unch), proj2(4,j_unch)
c
c     Jose: Entropy measure for chaos.
c
      logical ent_stat, j_en_flag
      real om_entj(nyp,nzd,(nxp/2+1))
      real om_prevr(nyp,nzd,(nxp/2+1)),om_previ(nyp,nzd,(nxp/2+1))
      real ent(nyp), dz, t_je
      integer en_ty_fl
      integer ent_cal, i_jen
      character*80 nam_entj

c     Jose: ent is ma multipurpose variable.
c
c     Jose: For time averaged mean flow
c
      real tav_u(nyp)
      integer tav_cnt,itag
      logical t_av

      integer req1,npget,k
#ifdef MPI
      integer status1(mpi_status_size), status2(mpi_status_size)
#endif
c     Large-eddy simulation (LES)
c     iiles can be 0 or 1 to determine whether the LES arrays are allocated.
c     iiles_smag (0 or 1) is to allocate arrays for the dynamic procedure.
c
      logical sgs
      ! Filtered velocity field / SGS force
      real gur(iiles*memnx*memny*memnz*5)
      real gui(iiles*memnx*memny*memnz*5)
      real gthr(iiles*memnx*memny*memnz*2*scalar)
      real gthi(iiles*memnx*memny*memnz*2*scalar)
      real gsr(iiles*memnx*memny*memnz*scalar)
      real gsi(iiles*memnx*memny*memnz*scalar)
      ! SGS stresses for velocity and scalars
      real taur(iiles*iiles_smag*memnx*memny*memnz*(6+3*scalar),nthread)
      real taui(iiles*iiles_smag*memnx*memny*memnz*(6+3*scalar),nthread)
      ! Temporary fields
      real gu3r  (iiles*nx/2*mbz*nyp  ,nthread)
      real gu3i  (iiles*nx/2*mbz*nyp  ,nthread)
      real ggu3r (iiles*nx/2*mbz*nyp  ,nthread)
      real ggu3i (iiles*nx/2*mbz*nyp  ,nthread)
      real gth3r  (iiles*nx/2*mbz*nyp  ,nthread)
      real gth3i  (iiles*nx/2*mbz*nyp  ,nthread)
      real ggth3r (iiles*nx/2*mbz*nyp  ,nthread)
      real ggth3i (iiles*nx/2*mbz*nyp  ,nthread)
      real gplane(iiles*nx/2*mbz*nyp  ,nthread)
      real tau3r (iiles*nx/2*mbz*nyp*6,nthread)
      real tau3i (iiles*nx/2*mbz*nyp*6,nthread)
      real tau2r (iiles*(nxp/2+1)*mby*nzd*6,nthread)
      real tau2i (iiles*(nxp/2+1)*mby*nzd*6,nthread)
      real s2r   (iiles*(nxp/2+1)*mby*nzd*6,nthread)
      real s2i   (iiles*(nxp/2+1)*mby*nzd*6,nthread)
      real l2r   (iiles*(nxp/2+1)*mby*nzd*6,nthread)
      real l2i   (iiles*(nxp/2+1)*mby*nzd*6,nthread)
      real ls2r   (iiles*(nxp/2+1)*mby*nzd*3*scalar,nthread)
      real ls2i   (iiles*(nxp/2+1)*mby*nzd*3*scalar,nthread)
      real m2r   (iiles*(nxp/2+1)*mby*nzd*6,nthread)
      real m2i   (iiles*(nxp/2+1)*mby*nzd*6,nthread)
      real b2r   (iiles*(nxp/2+1)*mby*nzd*6,nthread)
      real b2i   (iiles*(nxp/2+1)*mby*nzd*6,nthread)
      real gu2r  (iiles*(nxp/2+1)*mby*nzd*3,nthread)
      real gu2i  (iiles*(nxp/2+1)*mby*nzd*3,nthread)
      real gom2r (iiles*(nxp/2+1)*mby*nzd*3,nthread)
      real gom2i (iiles*(nxp/2+1)*mby*nzd*3,nthread)
      real gdu2r (iiles*(nxp/2+1)*mby*nzd*9,nthread)
      real gdu2i (iiles*(nxp/2+1)*mby*nzd*9,nthread)
      real sabsr (iiles*(nxp/2+1)*mby*nzd,nthread)
      real sabsi (iiles*(nxp/2+1)*mby*nzd,nthread)
      real sdr   (iiles*(nxp/2+1)*mby*nzd,nthread)
      real sdi   (iiles*(nxp/2+1)*mby*nzd,nthread)
      real cr    (iiles*(nxp/2+1)*mby*nzd,nthread)
      real ci    (iiles*(nxp/2+1)*mby*nzd,nthread)
      ! Jose: Varying viscosity additional working arrays
      ! Should I use more of these?
      real scalr   (iiles*(nxp/2+1)*mby*nzd,nthread)
      real scali   (iiles*(nxp/2+1)*mby*nzd,nthread)
      real nur   (iiles*(nxp/2+1)*mby,nthread)
      real nui   (iiles*(nxp/2+1)*mby,nthread)

      ! Filter
      real gewy1(nyp,5),gewy2(nyp,5)
      real gewspecx(nx/2+2)
      real gewspecz(nz/2+2)
      real diags(nyp,5)
      real filtxz(iiles*nx/2,nz)
      real lpfxz(iiles*nx/2,nz)
      ! Secondary highpass filter
      real gewy12(nyp,5),gewy22(nyp,5)
      real diags2(nyp,5)
      real filtxz2(iiles*nx/2,nz)
      real lpfxz2(iiles*nx/2,nz)
      integer iord,ihighorder
      real cutoff
      ! Diverse
      integer iles,ineg
      real chi,chi_scaled,cs,prt
c
c     Selective frequency damping (SFD)
c
      integer isfd
      logical sfdzero
      real sfd_chi,sfd_delta
      real sfdur(memnx,memny,memnz,iisfd*6)
      real sfdui(memnx,memny,memnz,iisfd*6)
c
c     Perturbation equations and linearised version
c
      logical lin,pert
      real bom1(nxp/2+1,nyp*iipert,3+2*scalar)
      real bom2(nxp/2+1,nyp*iipert,3+2*scalar)
      real bom3(nxp/2+1,nyp*iipert,3)
      real bom4(nxp/2+1,nyp*iipert,3)
      real appr((nxp/2+1)*iipert),appi((nxp/2+1)*iipert)
      real app2r((nxp/2+1)*iipert),app2i((nxp/2+1)*iipert)
      real wlin(nxp/2+1,nyp*iipert)
c
c     3D base flow
c
      logical bf3
      character*80 bf3nam
      real bf3ur(memnx,memny,memnz,6*iibf3)
      real bf3ui(memnx,memny,memnz,6*iibf3)
      real bf3u2r(nxp/2+1,nzd,nyp/nproc+1,6*iibf3)
      real bf3u2i(nxp/2+1,nzd,nyp/nproc+1,6*iibf3)
      real bf3tempr(nxp/2+1,nzd,3*iibf3)
      real bf3tempi(nxp/2+1,nzd,3*iibf3)
c
c     Arnoldi
c
      character*80 arnoldi_namnin,arnoldi_namnut
      real arnoldi_tend, arnoldi_dt, arnoldi_tstart
      real arpckin(nx*ny*nz/nproc*3*arnoldi)
      real arpckut(nx*ny*nz/nproc*3*arnoldi)
      integer maxn
c
c     Wrapper
c
      real curamp
      real tmin
      real alam,aturb
      logical tstatus
c
c     OpenMP
c
      integer ompnproc,ompnthreads,ompmaxthreads
#ifdef OPENMP
      integer omp_get_num_procs, omp_get_num_threads
      integer omp_get_thread_num, omp_get_max_threads
#endif
c
c     MPI
c
      integer all_nodes,my_node
      integer realg1,realg2
      integer, parameter :: nypp = nyp/nproc+(min(nproc,2)-1)
#ifdef MPI
#ifdef ALLTOALL
      integer realg1_vec,realg2_vec
      integer block_lengths(2),displacements(2),types(2)
      integer(kind=mpi_address_kind) displacements2(2)
#endif
      integer my_procs, my_threads
      character(len=mpi_max_processor_name) my_name
      character(len=mpi_max_processor_name) all_names(nproc)
      integer all_threads(nproc),all_procs(nproc)
      integer ierror
c     For MPI-2
c      integer irequired,iprovided
#endif

c-----------------------------------
c
c     Equivalence statements
c     **   mbox2 = (nxp/2+1)*nzd*mby*nthread
c     **   mbox3 = nx/2*nyp*mbz*nthread
c
c     These statements alias the work arrays in steps 2 and 3
c     onto each other. Therefore, a substantial memory saving
c     can be achieved.
c
c     The equivalence statements can safely be commented out
c     if not needed (or unsure).
c
c     Remaining big variables: plxy,wxy,wyz,wbr,wbi
c
c-----------------------------------
#ifdef MEMORY_SAVE
      real alloc2(mbox2,32+8*scalar+6*pressure)
      real alloc3(mbox3,47+21*scalar)

      equivalence (alloc2,alloc3)

      equivalence (alloc2(1,1),u2r)
      equivalence (alloc2(1,4),u2i)
      equivalence (alloc2(1,7),om2r)
      equivalence (alloc2(1,10),om2i)
      equivalence (alloc2(1,13),wr)
      equivalence (alloc2(1,14),wi)
      equivalence (alloc2(1,15),du2r)
      equivalence (alloc2(1,24),du2i)
      equivalence (alloc2(1,33),th2r)
      equivalence (alloc2(1,33+4*scalar),th2i)
      equivalence (alloc2(1,33+8*scalar),h2r)
      equivalence (alloc2(1,33+8*scalar+3*pressure),h2i)

      equivalence (alloc3(1,1),u3r,h3r)
      equivalence (alloc3(1,4),u3i,h3i)
      equivalence (alloc3(1,7),om3r)
      equivalence (alloc3(1,10),om3i)
      equivalence (alloc3(1,13),domyr,pomyr)
      equivalence (alloc3(1,14),domyi,pomyi)
      equivalence (alloc3(1,15),domyh,pvr)
      equivalence (alloc3(1,16),omyh,pvi)
      equivalence (alloc3(1,17),homyr,boxr)
      equivalence (alloc3(1,18),homyi,boxi)
      equivalence (alloc3(1,19),dvr)
      equivalence (alloc3(1,20),dvi)
      equivalence (alloc3(1,21),dvh)
      equivalence (alloc3(1,22),dvh2)
      equivalence (alloc3(1,23),vr)
      equivalence (alloc3(1,24),vi)
      equivalence (alloc3(1,25),hvr)
      equivalence (alloc3(1,26),hvi)
      equivalence (alloc3(1,27),w3)
      equivalence (alloc3(1,28),d2omyh)
      equivalence (alloc3(1,29),d2v)
      equivalence (alloc3(1,33),f)
      equivalence (alloc3(1,38+2*scalar),bis)
      equivalence (alloc3(1,44+4*scalar),th3r)
      equivalence (alloc3(1,44+5*scalar),th3i)
      equivalence (alloc3(1,44+6*scalar),hthr)
      equivalence (alloc3(1,44+7*scalar),hthi)
      equivalence (alloc3(1,44+8*scalar),pthr)
      equivalence (alloc3(1,44+9*scalar),pthi)
      equivalence (alloc3(1,44+10*scalar),dthr)
      equivalence (alloc3(1,44+11*scalar),dthi)
      equivalence (alloc3(1,44+12*scalar),c)
      equivalence (alloc3(1,45+12*scalar),d)
      equivalence (alloc3(1,46+12*scalar),e)
      equivalence (alloc3(1,47+12*scalar),fth_bc)
      equivalence (alloc3(1,47+14*scalar),d2th_bc)
      equivalence (alloc3(1,47+18*scalar),dth3r_bc)
      equivalence (alloc3(1,47+19*scalar),dth3i_bc)
      equivalence (alloc3(1,47+20*scalar),dth3r_bc2)
      equivalence (alloc3(1,47+21*scalar),dth3i_bc2)
#endif
c      common /fields/ ur,ui,u2r,u2i,h2r,h2i

c-----------------------------------------------------------------------
c     BEGIN OF PROGRAM
c-----------------------------------------------------------------------
c
c     Initialize timers
c
      call time_string(sttimec)
      call ctim(ctime0,wtime0)
      ctime3=ctime0
      wtime3=wtime0
c
c     Initialize array d2v
c
      d2v = 0.

      if (wrapper.eq.0) then
#ifdef MPI
c
c     Startup MPI
c
      call mpi_init(ierror)
c
c     For MPI-2
c      irequired = MP_THREAD_MULTIPLE
c      call mpi_init_thread(irequired,iprovided,ierror)
c      write(ioe,*) irequired,iprovided
c      stop

      call mpi_comm_rank(mpi_comm_world,my_node,ierror)
      call mpi_comm_size(mpi_comm_world,all_nodes,ierror)
#else
c
c     Set MPI variables to serial settings
c
      my_node = 0
      all_nodes = 1
#endif
      endif

      if (my_node.eq.0) then
         write(ios,*) '********************************************'//
     &              '***********************'
         write(ios,*) '*                                           '//
     &              '                      *'
         write(ios,*) '*                             Simson        '//
     &              '                      *'
         write(ios,*) '*                        bla $Rev: 1821 $   '//
     &              '                      *'
         write(ios,*) '*                                           '//
     &              '                      *'
         write(ios,*) '********************************************'//
     &              '***********************'

         write(ios,*)
         write(ios,*)
         write(ios,*) '>>>>>>  General information  <<<<<<<'
         write(ios,*) '--------------------------------------------'//
     &        '-----------------------'
         write(ios,*) 'Started              : ',sttimec

      end if

#ifdef MPI
c#ifndef ALLTOALL
c
c     Define vector types for normal getpxy/putpxy
c     
      call mpi_type_vector(memnz,memnx,memnx*memny,
     &     mpi_double_precision,realg1,ierror)
      call mpi_type_vector(memnz,memnx,nxp/2+1,
     &     mpi_double_precision,realg2,ierror)
      call mpi_type_commit(realg1,ierror)
      call mpi_type_commit(realg2,ierror)
c
c     MPI environment
c
      call mpi_findout(mpi_comm_world,all_names,my_name,
     &     all_procs,my_procs,all_threads,my_threads)

#endif

c
c     Set OpenMP variables to serial settings
c
      ompnproc = 1
      ompnthreads = 1
c
c     Step 1
c
c
c     Print compile-time parameters
c
      call ppar(my_node,all_nodes,ompnthreads,ompnproc)
c
c     Runge-Kutta/CN coefficients
c
      call rkcoeff(anrk,bnrk,cnrk)
c
c     Read run-time parameters
c
      call rparambl(namnin,namnut,
     &     tmax,dtn,dtmax,
     &     vart,nst,varsiz,rot,cdev,
     &     loctyp,fp1,fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7,
     &     fpds8,fpdds1,fpdds2,fpdds3,fpdds4,fpdds5,
     &     wbci,wp1,wp2,
     &     wpds1,wpds2,wpds3,wpds4,wpds5,
     &     wpds6,wpds7,wpds8,wpds9,wpds10,wpds11,wpds12,
     &     wpdds1,wpdds2,wallroughinput,
     &     tripf,tamps,tampt,txsc,tysc,nzt,tdt,seed,tx0,
     &     tsave,nmsave,nsave,cflmax,icfl,iamp,namamp,longli,ibc,cim,
     &     icorr,spat,tabfre,namfre,rbfl,nambfl,gall,
     &     fmax,fstart,fend,frise,ffall,
     &     ampob,amp2d,osmod,osdamp,osamp,osfil,iext,namext,ixys,ixyss,
     &     namxys,txys,maxit,kx,kz,mwave,nwave,namwav,
     &     ipl,nampl,npl,tpl,cpl,mpl,xl,cpumax,wallmax,
     &     corrf,corrf_x,corrnam,corrnam_x,
     &     corrx,corry,corrz,corry_x,corrt,corrt_x,ncorr,ncorr_x,
     &     write_inter,namnut_inter,my_node,suction,asbl,vsuc,
     &     streak,iampst,ampst,tsmoo,tsmst,tsmend,phist,
     &     str_nam1,str_nam2,isfd,sfd_delta,sfd_chi,sfdzero,namnut_sfd,
     &     sgs,waves,waamp,wamoo,wamst,waiamp,pert,lin,
     &     namser,serc,sert,serf,nser,tbc,
     &     dtheta0_upp,dtheta0_low,theta0_low,theta0_upp,cflux,retau,
     &     iimhd,imhd,b0,mhd_n,fileurms,bf3,bf3nam,
     &     jet_prof,jet_diam,xjet,zjet,jetmag,delT)

      if (my_node.eq.0) then
         write(ios,*)
         write(ios,*)
         write(ios,*) '>>>>>>  Initialization  <<<<<<<'
         write(ios,*) '--------------------------------------------'//
     &        '-----------------------'
      end if

      call rdiscbl(ur,ui,u0low,u0upp,w0low,w0upp,du0upp,re,pr,
     &           xlold,zl,t,xs,dstar,fltype,bstart,blength,rlam,m1,
     &           spanv,spat,namnin,varsiz,3+scalar,boxr,boxi,w3,urx,
     &           nxtmp,my_node,gr,ivarvisc)

c     Set to zero for spatial cases
c     Will however be overwritten by the values computed from the
c     base flow (bflow.f)
c
      if (abs(u0low).lt.1e-8) u0low = 0.
      if (abs(w0low).lt.1e-8) w0low = 0.
      if (abs(u0upp).lt.1e-8) u0upp = 0.
      if (abs(w0upp).lt.1e-8) w0upp = 0.
c
c     Introduce boundary-layer scaling
c
      call rescale(tmax,dtn,dtmax,rot,tsave,nsave,dstar,xl,
     &     fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7,
     &     fpds8,fpdds1,fpdds2,fpdds3,fpdds4,fpdds5,fpdds6,
     &     wpds1,wpds2,wpds3,wpds4,wpds5,
     &     wpds6,wpds7,wpds8,wpds9,wpds10,wpds11,wpds12,
     &     wpdds1,wpdds2,
     &     tamps,tampt,txsc,tysc,tdt,txys,tx0,
     &     fmax,fstart,fend,frise,ffall,tsmst,tsmend,tsmoo,
     &     sfd_delta,sfd_chi,dtheta0_low,dtheta0_upp)

      xl=xlold
c
c     Compute some constant values for the simulation
c
      call preprbl(alfa,beta,eta,deta,xl,zl,prex,prey,prez,pres,prea,
     &     prexn,prezn,presn,prean,prezr,wint,wd1,
     &     gridx,gridy,gridz,dstar,deltaxyz2)

      x0 = 0.

c
c     Jose: fft tests. For the sake of understanding
c
c      call fft_test(gridx,gridz,prexn,prezn,prex,prez,alfa,beta,xl,zl)

c      call stopnow(9801)


c
c     Jose: Entropy initial file
c
      j_en_flag=.false.
      i_jen = 0
      t_je = 0.
      call r_jent(j_en_flag,en_ty_fl,nam_entj,i_jen,t_je)
      if(j_en_flag) then 
         dz = abs(gridz(1)-gridz(2))
         if(my_node.eq.0) then
            open(unit=42,file=nam_entj)
         end if
      end if

c
c     Jose: Unstable mode check. Reading parameters
c
      j_un_chk=.false.
      call r_un_chk(j_unch,j_un_chk,n_un,nam_un_chk,un_alp,un_bet,
     &     vel_unchk_r,vel_unchk_i,wint)

c
c     Jose: Reference temperature for varying viscosity case.
c
      if (ivarvisc.eq.1) then
c     Reference viscosity defined by a specific value of T
c     corresponding to hot or cold wall whichever is greater
c     The lowest temperature is always considered to 295 K
         if (delT.lt.0.) then
c Unstable stratification. 
            Tref = 295 - delT
            Tl = Tref
         else
c Stable stratification
            Tref = 295 + delT
            Tl = 295
         end if
      elseif (ivarvisc.eq.2) then
c     Reference viscosity defined by a mean value of T
c     i.e Tref = 0.5*(Tu + Tl)
c     The lowest temperature is always considered to 295 K
         if (delT.lt.0.) then
c Unstable stratification. 
            Tref = 295 - delT
            Tl = Tref
            Tu = 295
         else
c Stable stratification
            Tref = 295 + delT
            Tl = 295
            Tu = Tref
         end if
         Tref = 0.5*(Tu + Tl)
      end if

c      call stopnow(31813)

c$$$      write(*,*) 'j_en_flag: ', j_en_flag
c      write(*,*) 'en_ty_fl: ', en_ty_fl
c$$$      write(*,*) 'nam_entj: ', nam_entj
c$$$      write(*,*) 't_je: ', t_je
c$$$      write(*,*) 'i_jen: ', i_jen
c$$$
c      call stopnow(88)

c      write(*,*) 'the adjusted probe points.' 

c
c     Initialize time series
c
      if (serf) then
         do i=1, nser
c
c     x-coordinate
c
            sercc(1,i) = serc(1,i)*dstar
            serci(1,i) = modulo(nint(sercc(1,i)/xl*real(nx))+nx/2,nx)+1
            serci(4,i) = mod(serci(1,i),2)
            serci(1,i) = (serci(1,i)+1)/2
            sercc(1,i) = nint(serc(1,i)/xl*real(nx))*xl/real(nx)
c
c     y-coordinate
c     input should be from 0-2.
c
            sercc(2,i) = serc(2,i)*dstar-1.
            sercc(2,i) = max(sercc(2,i),-1.)
            sercc(2,i) = min(sercc(2,i),1.)
            serci(2,i) = int(acos(sercc(2,i))/pi*real(nyp-1)+1.5)
            sercc(2,i) = (cos(pi*real(serci(2,i)-1)/
     &           real(nyp-1))+1.)/dstar
c
c     z-coordinate
c
            sercc(3,i) = serc(3,i)*dstar
            serci(3,i) = modulo(nint(sercc(3,i)/zl*real(nz))+nz/2,nz)+1
            sercc(3,i) = nint(serc(3,i)/zl*real(nz))*zl/real(nz)

c            write(ios,*) i,serci(2,i),sercc(2,i)
         end do
c
c     Check array size
c
         nsamp = ixyss/ixys+min(1,mod(ixyss,ixys))
         if (nsamp.gt.msamp) then
            if (my_node.eq.0) then
               write(ioe,*) 'increase msamp to at least ',nsamp
            end if
            call stopnow(56455435)
         end if
         series = 0.
      end if
      sercount = 0
      
c      call stopnow(42)
c
c     xblc: current position for parallel flows in terms of reference
c           velocity and current time
c     xs:   shifting of the box since time zero
c     xs0:  shifting of the box since time zero (at beginning of simulation)
c
      xblc = x0+cdev*t
      xs0  = xs
      if (abs(u0low).gt.1.e-12) xs0=0.
c      write(*,*)'x0 = ', x0
c      write(*,*)'cdev = ', cdev, ', xs = ', xs

      if (my_node.eq.0) then
         write(ios,*)
         write(ios,*) 'Boundary conditions (including shifts)'
         write(ios,*) '  u0low,w0low : ',u0low,w0low
         write(ios,*) '  u0upp,w0upp : ',u0upp,w0upp
         write(ios,*) '  du0upp      : ',du0upp*dstar
      end if

c
c     Get time step dt for first iteration
c
      if (vart) then
         call getdt(dtn,cflmax,rot,deta,xl,zl,prex,prez,pres,prea,
     &        ur,ui,u2r,u2i,bu1,bu2,pert,lin,wr,wi,
     &        realg1,realg2,wbci,bf3,bf3u2r,bf3u2i,my_node,
     &        bf3tempr,bf3tempi)
         if (my_node.eq.0) then
            write(ios,*)
            write(ios,*) 'First time step      : ',dtn/dstar
         end if
      end if
      
c      write(*,*) 'Jose here. Line 751. '
c      call stopnow(18)

      if (nst.gt.2) then
         an   = dtn*anrk(1,nst-2)
         cnp1 = 0.
      else
         an   = dtn
         cnp1 = an*.5
      end if
      dtnp1 = dtn
      anp1  = an
      bnp1  = 0.
c
c     Large-Eddy Simulation
c
      iles = 0

c     Build partial rhs for first Euler step (i.e. next step)
c     filtering in wall-normal direction in gur,gui
c
      call prhs(ur,ui,puw,alfa,beta,re,pr,anp1,prey,u3r,u3i,om3r,om3i,
     &     th3r,th3i,
     &     pomyr,pomyi,pvr,pvi,pthr,pthi,boxr,boxi,w3,dthr,dthi,
     &     iles,gur,gui,taur,taui,gu3r,gu3i,ggu3r,ggu3i,
     &     gthr,gthi,gsr,gsi,gth3r,gth3i,ggth3r,ggth3i,
     &     gewy1,gewy2,gewy12,gewy22,iord,
     &     diags,diags2,gplane,
     &     filtxz,filtxz2,lpfxz,lpfxz2,my_node,ihighorder,cs)
c
c     Initialize some variables
c
      sumw=  0.
      xys = 0.
      xysth = 0.
      stop_now = 0
      it = 1
      zs = 0.
      ctime_stat = 0.
      wtime_stat = 0.
      ctime_les = 0.
      wtime_les = 0.
      ctime_nonlin = 0.
      wtime_nonlin = 0.
      ctime_linear = 0.
      wtime_linear = 0.
      ctime_io = 0.
      wtime_io = 0.
      ctime_tau = 0.
      wtime_tau = 0.
      ctime_rest = 0.
      wtime_rest = 0.
      ctime_trip = 0.
      wtime_trip = 0.
      ctime_ffun = 0.
      wtime_ffun = 0.
      ctime_cfl = 0.
      wtime_cfl = 0.
      ctime_planes = 0.
      wtime_planes = 0.
      ctime_press = 0.
      wtime_press = 0.
      ctime_part = 0.
      wtime_part = 0.
      tcom = 0.
      tser = 0.
      delta_wall_time = 0.

      if(j_en_flag) then 
         ent_cal = 0
         om_entj = 0.
         om_prevr = 0.
         om_previ = 0.
      end if

      t_av = .false.
      tav_cnt = 0
      tav_u = 0.

      if(j_un_chk) then
         j_un_cnt = 0
         if(my_node.eq.0) then
            do ybp = 1,n_un
               yb = 48 + ybp
               fil_no = achar(yb)
               fil_nam = fil_no//'_'//trim(nam_un_chk)
c               write(*,*) fil_nam
c               open(unit=yb,file=fil_nam)
            end do
         end if

         vel_unchk_pr = 0.
         vel_unchk_pi = 0.
      end if

c      call stopnow(21297)

      do_stat = .false.
      do_press = .false.
      do_press_nl = .false.
      ent_stat=.false.
c
c     Pressure gradient for temporal channel flow (fltype=1)
c     For all other cases (e.g. boundary-layer flows) it is set to zero.
c
c     AM
c     If pertubation mode then set mflux and px to zero
c     :: No pressure gradient :: the flow is driven by the base flow
c
      if (fltype.eq.1.and..not.pert) then

         if (cflux) then
c
c     Set mflux for fixed mass flux
c     px is only used as a first guess
c
            mflux = 4./3.
            px = -2./re
            if (my_node.eq.0) then
               write(ios,*)
               write(ios,*) 'Fixed mass flux for channel flow : '
               write(ios,*) '   mflux  : ',mflux
            end if
         else
c
c     Set px for fixed pressure gradient
c
c     The pressure gradient in laminar channel flows (with Re based
c     on the centre-line velocity), the skin friction is given as
c            retau = sqrt(2.*re)
c
            mflux = 0.
            px = -(retau/re)**2
            retau = sqrt(-px)*re
            if (my_node.eq.0) then
               write(ios,*)
               write(ios,*) 'Fixed pressure gradient ' //
     &              'for channel flow : '
               write(ios,*) '   px     : ',px
               write(ios,*) '   Re_tau : ',retau
            end if
         end if
      else
         px = 0.
         mflux = 0.
         cflux = .false.
         if (my_node.eq.0) then
            write(ios,*) 'No pressure gradient : ',px,mflux,cflux
         end if
      end if
c
c     Open files for amplitude, waves and extrema
c
      if (my_node.eq.0) then
         if (iamp.gt.0.and.wrapper.gt.0) then
            open(unit=15,file=namamp,position='append')
         else if (iamp.gt.0) then
            open(unit=15,file=namamp)
         end if
         if (mwave.gt.0) open(unit=16,file=namwav)
         if (iext .gt.0) open(unit=18,file=namext)
      end if
c
c     Determine next full flow field to save
c
      isave = 1
      j=0
      do i=1,nsave
         if (t.ge.tsave(i).and.j.eq.0) then
            isave=isave+1
         else
            j=1
         end if
      end do
      if (my_node.eq.0) then
         write(ios,*) 'Next field to save   : isave =',isave,
     &        ' at t=',tsave(isave)/dstar
      end if

#ifdef MPI
      if (nproc.gt.1) call mpi_barrier(mpi_comm_world,ierror)
#endif

      sttimet = t

      if (my_node.eq.0) then

         call ctim(ctime4,wtime4)
         wtime_startup = wtime4-wtime3
         ctime_startup = ctime4-ctime3
         ctime3=ctime4
         wtime3=wtime4
         ctime1=ctime4
         wtime1=wtime4

         write(ios,'(a,f10.3,a,f10.3,a)') ' CPU time             : ',
     &        ctime4-ctime0,
     &        ' sec'
         write(ios,'(a,f10.3,a,f10.3,a)') ' Wall time            : ',
     &        wtime4-wtime0,' sec'
         write(ios,*)
         write(ios,*)
         write(ios,*) '>>>>>>  Starting simulation  <<<<<<'
         write(ios,*) '--------------------------------------------'//
     &              '-----------------------'
         write(ios,*) 'Starting iteration at t=',t/dstar
         write(ios,*)
c
c     Write out first step file
c
         open(file='step.out',unit=87)
         write(87,'(a20,f18.9,a,f18.9,a)') 'Time now: ',t/dstar,
     &        ' (max. ',tmax/dstar,')'
         write(87,'(a20,f18.9)') 'Started with: ',sttimet/dstar
         write(87,'(a20,f18.9)') 'Delta: ',(t-sttimet)/dstar
         write(87,'(a20,i18,a,i18,a)') 'Iterations: ',it-1,
     &        ' (max. ',maxit,')'
         call time_string(timdat)
         write(87,'(a20,a20)') 'Time now: ',timdat
         write(87,'(a20,a20)') 'Started: ',sttimec
         write(87,'(a,f18.9,a,f18.9,a)') 'CPU time       ',
     &        ctime4-ctime0,
     &        ' seconds (max. ',cpumax,')'
         write(87,'(a,f18.9,a,f18.9,a)') 'Wallclock time ',
     &        wtime4-wtime0,
     &        ' seconds (max. ',wallmax,')'
         close(87)
      end if

c------------------------------------
c>>>>>>  START OF INTEGRATION  <<<<<<
c------------------------------------

c
c     Timestep loop, starting with it=1
c
 1    continue
c
c     Step 2
c
c     Time of the substep
c
      rksubstep = mod(it-1,nst)+1
      tc=t+cnp1

c      write(ios,*) '>>>>>> RK STEP <<<<<<',rksubstep,it
c      write(ios,*) '      step: ',t/dstar,'-',(t+dtn)/dstar
c      write(ios,*) '   substep: t=',tc/dstar
c
c     Move centre of the box according to u0low, w0low
c
      xsc=xs-cnp1*u0low
      zsc=zs-cnp1*w0low
c
c     Only needed for temporal boundary-layer simulations:
c     Streamwise position of the boundary layer
c
      xbl =x0+cdev*t
      xblc=x0+cdev*tc
c
c     Decide whether statistics are taken
c
      if (ixys.gt.0.and.mod(it-1-nst,ixys).eq.0..and.t.ge.txys) then

         if (t.ge.txys) then
            do_stat = .true.
         end if
      else
         do_stat = .false.
      end if

      call ctim(ctime4,wtime4)
      ctime_rest = ctime_rest + (ctime4-ctime3)
      wtime_rest = wtime_rest + (wtime4-wtime3)
      ctime3=ctime4
      wtime3=wtime4
c
c     Write planes
c
      if (npl.gt.0.and.mod(it-1,ipl).eq.0) then

         if (nproc.eq.1) then
c
c     Serial plane writing
c
            call wplbl_serial(npl,nampl,tpl,cpl,
     &           it,u2r,u2i,wxz,wxy,wyz,wr,wi,ur,ui,
     &           mpl,re,xl,zl,t,xs,dstar,fltype,
     &           prexn,prezn,presn,prean,u0low,w0low,alfa,zs,
     &           beta,wplblodd)
         else
c
c     Parallel plane writing
c
            call wplbl(npl,nampl,tpl,cpl,
     &           it,u2r,u2i,wxz,wxy,wyz,wr,wi,ur,ui,
     &           mpl,re,xl,zl,t,xs,dstar,fltype,
     &           prexn,prezn,presn,prean,u0low,w0low,alfa,zs,
     &           beta,wplblodd,my_node,realg1,realg2,plxy)
         end if

      end if

      call ctim(ctime4,wtime4)
      ctime_tau = ctime_tau + (ctime4-ctime3)
      wtime_tau = wtime_tau + (wtime4-wtime3)
      ctime3=ctime4
      wtime3=wtime4

      taur = 0.
      taui = 0.

      if (ivarvisc.gt.0.and.delT.ne.0) then
         do ybp=1,nypp

            ip = 1
#ifdef OPENMP
            ip = OMP_GET_THREAD_NUM()+1
#endif
            yb = (ybp-1)*nproc+my_node+1
            myb= (yb-1)/mby+1 
            call varvisc(taur(1,ip),taui(1,ip),ivarvisc,yb,alfa,beta,
     &           u0low,w0low,prexn,prezn,presn,prean,prex,prez,
     &           my_node,realg1,realg2,ur,ui,
     &           re,eta,it,cr(1,ip),ci(1,ip),nur(1,ip),nui(1,ip),
     &           tau3r(1,ip),tau3i(1,ip),s2r,s2i,l2r,l2i,m2r,m2i,
     &           om2r,om2i,du2r,du2i,u2r,u2i,delT,Tl,Tref,xl,zl)

         enddo
      end if

c      call stopnow(9801)

      if (iles.eq.2.or.iles.eq.3.or.ivarvisc.ne.0) then
c
c     Compute sigma_i=dtau_ij/dx_j
c
!$omp parallel do private(zbp,zb,ip)
         do zbp=1,memnz
#ifdef OPENMP
            ip = OMP_GET_THREAD_NUM()+1
#else
            ip = 1
#endif
            zb = my_node*memnz+zbp
            call sgs_stress(taur,taui,zb,zbp,alfa,beta,
     &           gur,gui,gu3r(1,ip),gu3i(1,ip),w3(1,ip),prey,
     &           u3r(1,ip),u3i(1,ip),tau3r(1,ip),tau3i(1,ip),
     &           gsr,gsi)
         end do

      end if
c-------------------------------------

      call ctim(ctime4,wtime4)
      ctime_tau = ctime_tau + (ctime4-ctime3)
      wtime_tau = wtime_tau + (wtime4-wtime3)
      ctime3=ctime4
      wtime3=wtime4

c-------------------------------------
c     STATISTICS
c-------------------------------------
c     independent xz planes in y
c     ybp    local index of xz plane
c     ip     thread number
c     yb     start index of box
c     mby    thickness of box
c     myb    box number
c-------------------------------------
      if (do_stat) then
         sercount = sercount + 1
         if (my_node.eq.0) then
            write(ios,*) '** Calculate statistics **', t/dstar,sercount
         end if

!$omp parallel do private(ybp,ip,yb,myb)
         do ybp=1,nypp
#ifdef OPENMP
            ip = OMP_GET_THREAD_NUM()+1
#else
            ip = 1
#endif
            yb = (ybp-1)*nproc+my_node+1
            myb= (yb-1)/mby+1
            call boxxys(xys,xysth,tc,xs,zs,yb,dtn,ur,ui,
     &           u0low,w0low,prexn,prezn,presn,prean,
     &           alfa,beta,u2r(1,ip),u2i(1,ip),
     &           h2r(1,ip),h2i(1,ip),
     &           om2r(1,ip),om2i(1,ip),
     &           wr(1,ip),wi(1,ip),
     &           du2r(1,1,1,1,ip),du2i(1,1,1,1,ip),
     &           th2r(1,ip),th2i(1,ip),
     &           corrf,corrf_x,corrxi,corryi,corrzi,corryi_x,corr,
     &           corr_x,corrwsave,corrwsave_x,ncorr,ncorr_x,corrxodd,
     &           corrt,corrt_x,corr_ms,corr_ms_x,my_node,realg1,realg2,
     &           ybp,taur,taui,tau2r(1,ip),tau2i(1,ip),iles,
     &           series,serci,sert,serf,nser,sercount,re,fltype,dstar,
     &           pert,bu3,bu4,bom3,bom4,
     &           fmhdr,fmhdi,j2r(1,ip),j2i(1,ip),b0,imhd,px,mflux,mhd_n)
         end do
         sumw=sumw+dtn
      end if
c-------------------------------------

      call ctim(ctime4,wtime4)
      ctime_stat = ctime_stat + (ctime4-ctime3)
      wtime_stat = wtime_stat + (wtime4-wtime3)
      ctime3=ctime4
      wtime3=wtime4

      if (ixys.gt.0.and.mod(it-1,ixyss).eq.0.and.t.ge.txys
     &     .and.it.gt.1.and.sumw.gt.0.) then
c
c     Write statistics and intermediate velocity field
c
         if (my_node.eq.0) then
            write(ios,*) '** Write statistics **', t/dstar,sumw/dstar
         end if
         call wxys(xys,xysth,
     &        totcorr,totcorr_x,corr,corr_x,ncorr,ncorr_x,corrf,corrf_x,
     &        totcorr_ms,totcorr_ms_x,corr_ms,corr_ms_x,my_node,
     &        series,serf,nser,totseries,
     &        sumw,re,xl,zl,t,dstar,fltype,
     &        bstart,blength,rlam,spanv,namxys,wxy,
     &        corrnam,corrnam_x,corrxx,corryy,corrzz,corryy_x,corrwsave,
     &        corrwsave_x,corrx,corry,corrz,corry_x,corrxi,corryi,
     &        corrzi,corryi_x,pr,m1,corrt,corrt_x,nsamp,namser,sercount,
     &        mhd_n,b0,it,ixyss)

         sercount = 0

         if (write_inter) then
            if (my_node.eq.0) then
               write(ios,*) '** Write intermediate velocity field **',
     &              t/dstar
            end if

            call wdiscbl(ur,ui,re,pr,m1,xl,zl,t,xs,dstar,fltype,
     &           bstart,blength,rlam,spanv,namnut_inter,3+scalar,gall,
     &           boxr,boxi,urx,alfa,zs,beta,my_node,gr)
         end if
      end if

      call ctim(ctime4,wtime4)
      ctime_io = ctime_io + (ctime4-ctime3)
      wtime_io = wtime_io + (wtime4-wtime3)
      ctime3=ctime4
      wtime3=wtime4

c######################################

c      tcom = 0.
c      tser = 0.

c-------------------------------------
c     NONLINEAR PART (explicit)
c-------------------------------------
c     independent xz planes in y
c     ybp    local index of xz plane
c     ip     thread number
c     yb     start index of box
c     mby    thickness of box
c     myb    box number
c-------------------------------------
!$omp parallel do private(ybp,ip,yb,myb)

      if(j_en_flag.and.(t.gt.t_je).and.i_jen.gt.0.and.
     &     mod(it-1-nst,i_jen).eq.0) then 
         ent_stat=.true.
      end if

      do ybp=1,nypp
#ifdef OPENMP
         ip = OMP_GET_THREAD_NUM()+1
#else
         ip = 1
#endif
         yb =(ybp-1)*nproc+my_node+1
         myb=(yb-1)/mby+1
         call nonlinbl(amp,campw,kx,kz,nwave,cflp(myb),vext,cext,
     &        iext,ur,ui,tc,xsc,zsc,xs,yb,rot,spat,fltype,
     &        fstart,fend,nst,
     &        bu1,bu2,
     &        loctyp,fp1,fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,
     &        fpds7,fpds8,fpdds1,fpdds2,fpdds3,fpdds4,fpdds5,
     &        tripf,txsc,tysc,ttzc,tx0,
     &        osmod,osnumb,
     &        osur,osui,osvr,osvi,oswr,oswi,
     &        xc1,xc2,fring1,fring2,cphas1,cphas2,sphas1,sphas2,
     &        evr,evi,eur,eui,ewr,ewi,afw,alfaw,betaw,ampob,amp2d,
     &        ev2r,ev2i,eu2r,eu2i,afw2d,alf2d,
     &        fbla,nbla,dybla,rlam,cdev,re,xblc,alfa,
     &        beta,eta,deta,xl,zl,prex,prez,pres,prea,it,icfl,
     &        iamp,u2r(1,ip),u2i(1,ip),om2r(1,ip),om2i(1,ip),
     &        th2r(1,ip),th2i(1,ip),
     &        wr(1,ip),wi(1,ip),
     &        iles,gur,gui,chi_scaled,prt,
     &        my_node,realg1,realg2,
     &        ampst,streak,betast,omegast,
     &        ndxst,uust_r,uust_i,vvst_r,vvst_i,wwst_r,wwst_i,
     &        tsmoo,tsmst,tsmend,iampst,phist,
     &        isfd,sfd_chi,sfd_delta,sfdur,sfdui,anp1,bnp1,wbci,
     &        pert,lin,bom1,bom2,x0,spanv,rv1r,rv1i,do_press_nl,
     &        tcom,tser,gsr,gsi,fmhdr,fmhdi,imhd,
     &        du2r(1,1,1,1,ip),du2i(1,1,1,1,ip),b0,mhd_n,
     &        bf3,bf3u2r,bf3u2i,ybp,bf3tempr,bf3tempi,pint,
     &        pr,gr,dstar,ivarvisc,ent_stat,ent_cal,en_ty_fl,om_entj,
     &        om_prevr,om_previ,ent,tav_u,t_av,
     &        j_un_chk,n_un,un_alp,un_bet,vel_unchk_pr,vel_unchk_pi)
      end do

c
c     To get the count for time averaged mean flow
c
      if(mod(it-1-nst,iamp).eq.0.and.iamp.gt.0.and.t_av) then
         tav_cnt = tav_cnt + 1
      end if

c      if(tav_cnt.gt.25) then
c         call stopnow(107)
c      end if
c      tav_u = 0.*tav_u

c
c     Communication and advection of particles
c
      call ctim(ctime4,wtime4)
      ctime_nonlin = ctime_nonlin + (ctime4-ctime3)
      wtime_nonlin = wtime_nonlin + (wtime4-wtime3)
      ctime3=ctime4
      wtime3=wtime4

c
c     Write amplitude file
c
      if (mod(it-1-nst,iamp).eq.0.and.iamp.gt.0) then
         call wampbl(t,it,amp,campw,kx,kz,
     &        mwave,nwave,wint,re,xl,zl,longli,dstar,
     &        rlam,xbl,my_node,eta,fbla,dybla,nbla,spanv,fileurms,
     &        fltype,px,gridy)
      end if
c
c     Jose: Evaluate "entropy" and write to file
c
c$$$      write(*,*) '#^*&^*#!^(*#!($^((^$(^!$'
c$$$      write(*,*) '#^*&^*#!^(*#!($^((^$(^!$'
c$$$      write(*,*) 'HERE I AM'
c$$$      write(*,*) '#^*&^*#!^(*#!($^((^$(^!$'
c$$$      write(*,*) '#^*&^*#!^(*#!($^((^$(^!$'
c$$$      write(*,*) 'ent_stat:  ', ent_stat
      if(ent_stat) then
         ent_cal = ent_cal + 1
         call entropy_j(om_entj,om_prevr,om_previ,ent_cal,ent,
     &        en_ty_fl,eta,dz,my_node,wint)
c         write(*,*) 'ent_cal: ', ent_cal
c         write(*,*) 'my_node: ', my_node
         
         if(ent_cal.gt.1.and.my_node.eq.0) then
c            write(*,*) 'ent_cal: ', ent_cal
c            write(*,*) 'en_ty_fl: ', en_ty_fl
            if(en_ty_fl.eq.1) then
               write(42,'(7(e25.16e3,tr1))') t/dstar,(ent(j),j=1,6)
            elseif(en_ty_fl.eq.2) then
c               write(*,*) '#^*&^*#!^(*#!($^((^$(^!$'
c               write(*,*) '#^*&^*#!^(*#!($^((^$(^!$'
c               write(*,*) 'HERE I AM'
c               write(*,*) '#^*&^*#!^(*#!($^((^$(^!$'
c               write(*,*) '#^*&^*#!^(*#!($^((^$(^!$'
               write(42,'(2(e25.16e3,tr1))') t/dstar,ent(1)
            end if
         end if
         ent_stat=.false.
      end if
c
c     Jose: Unstable mode checks
c
      if(mod(it-1-nst,iamp).eq.0.and.iamp.gt.0.and.j_un_chk) then
         en = 0.
         proj1 = 0.
         proj2 = 0.
         j_un_cnt = j_un_cnt + 1
         call j_unmod_cmp(j_un_cnt,n_un,vel_unchk_r,vel_unchk_i,
     &        vel_unchk_pr,vel_unchk_pi,vel_unchk_lr,vel_unchk_li,
     &        wint,en,proj1,proj2,my_node)
c
c     Write out the projection statistics
c
         if(my_node.eq.0) then
c            write(*,*) t/dstar,(en(j),proj1(1,j),proj1(2,j),proj2(1,j),
c     &           proj2(2,j),proj1(3,j),proj1(4,j),proj2(3,j),proj2(4,j),
c     &           j=1,n_un)            
            do ybp = 1,n_un
               yb = 48 + ybp
               write(*,'(10(e25.16e3,tr1))') t/dstar,en(ybp),
     &              proj1(1,ybp),proj1(2,ybp),proj2(1,ybp),proj2(2,ybp),
     &              proj1(3,ybp),proj1(4,ybp),proj2(3,ybp),proj2(4,ybp)               
            enddo

c            write(43,*) t/dstar,(en(j),proj1(1,j),proj1(2,j),proj2(1,j),
c     &           proj2(2,j),proj1(3,j),proj1(4,j),proj2(3,j),proj2(4,j),
c     &           j=1,n_un)
         end if
         vel_unchk_pr = 0.
         vel_unchk_pi = 0.
      end if

c      if(j_un_cnt.gt.1) call stopnow(138194)
c
c     Write extremum file
c
      if (mod(it-1,iext).eq.0.and.iext.gt.0) then
         call wextbl(t,vext,cext,re,xl,
     &        xbl,u0low,w0low,longli,dstar,eta,fltype,fend,
     &        my_node,fbla,dybla,nbla)
      end if

      call ctim(ctime4,wtime4)
      ctime_io = ctime_io + (ctime4-ctime3)
      wtime_io = wtime_io + (wtime4-wtime3)
      ctime3=ctime4
      wtime3=wtime4
c
c     Update cfl time-step restriction
c
      if (mod(it-1,icfl).eq.0.and.icfl.gt.0) then
         call wcflbl(dtn,cfl,cflp,rot,dstar,my_node)
      end if
c
c     Step 3
c
      if (vart.and.mod(it,nst).eq.0) then
         ! time step for NEXT iteration
         if(tmax.gt.0) then
            tleft = min(tmax,tsave(isave))-t-dtn
         else
            tleft = tsave(isave)-t-dtn
         end if
         if (cfl/cflmax.gt.10./9.) then
            dtnp1 = dtn* cflmax/cfl
         else
            dtnp1 = dtn*( 0.8 + 0.2*cflmax/cfl )
         end if
         dtnp1=min(dtnp1,dtmax)
         if (.5*dtnp1.lt.tleft.and.tleft.lt.4.*dtnp1) then
            if (my_node.eq.0) then
              if(tmax.gt.0) then
                 write(ios,*) 'Adjusting time step to ',
     &               min(tmax,tsave(isave))/dstar
              else
                 write(ios,*) 'Adjusting time step to ',
     &               tsave(isave)/dstar
              end if
            end if
            dtnp1=tleft/aint(tleft/dtnp1+1.)
         end if
      end if
c
c     RK integration weights for this iteration
c
      an=anp1
      bn=bnp1
c
c     RK integration weights for next iteration
c
      if (nst.eq.1) then
         bnp1=-dtnp1**2/(2.*dtn)
         anp1=dtnp1-bnp1
         cnp1=dtnp1*.5
      else
         anp1=dtnp1*anrk(mod(it,nst)+1,nst-2)
         bnp1=dtnp1*bnrk(mod(it,nst)+1,nst-2)
         cnp1=dtnp1*cnrk(mod(it,nst)+1,nst-2)
      end if

      call ctim(ctime4,wtime4)
      ctime_cfl = ctime_cfl + (ctime4-ctime3)
      wtime_cfl = wtime_cfl + (wtime4-wtime3)
      ctime3=ctime4
      wtime3=wtime4

c-------------------------------------
c     LINEAR PART
c-------------------------------------
c     independent xy planes in z
c     zbp    local index of xy plane
c     ip     thread number
c     zb     start index of box
c     mbz    thickness of box
c-------------------------------------

!$omp parallel do private(zbp,zb,ip)
      do zbp=1,memnz
#ifdef OPENMP
         ip = OMP_GET_THREAD_NUM()+1
#else
         ip = 1
#endif
         zb = my_node*memnz+zbp
         call linearbl(ur,ui,puw,zb,an,bn,anp1,bnp1,re,pr,alfa,beta,
     &        u0low,u0upp,w0low,w0upp,du0upp,px,ibc,cim,icorr,
     &        bu1jr,bu1ji,prey,h3r(1,ip),h3i(1,ip),
     &        hthr(1,ip),hthi(1,ip),pthr(1,ip),pthi(1,ip),
     &        dthr(1,ip),dthi(1,ip),
     &        u3r(1,ip),u3i(1,ip),om3r(1,ip),om3i(1,ip),
     &        th3r(1,ip),th3i(1,ip),f(1,1,ip),
     &        f(1,1,ip),f(1,2,ip),f(1,3,ip),f(1,4,ip),f(1,5,ip),
     &        f(1,4,ip),f(1,5,ip),f(1,6,ip),f(1,7,ip),
     &        bis(1,1,ip),bis(1,2,ip),
     &        bis(1,3,ip),bis(1,4,ip),bis(1,2,ip),
     &        bis(1,3,ip),bis(1,4,ip),bis(1,1,ip),bis(1,5,ip),
     &        bis(1,6,ip),bis(1,7,ip),bis(1,8,ip),
     &        domyr(1,ip),domyi(1,ip),domyh(1,ip),
     &        omyh(1,ip),d2omyh(1,ip),
     &        vr(1,ip),vi(1,ip),vh(1,ip),vh2(1,ip),
     &        dvr(1,ip),dvi(1,ip),dvh(1,ip),dvh2(1,ip),d2v(1,1,ip),
     &        d2v(1,2,ip),d2v(1,3,ip),d2v(1,4,ip),d2v(1,1,ip),
     &        d2v(1,2,ip),d2v(1,3,ip),homyr(1,ip),homyi(1,ip),
     &        hvr(1,ip),hvi(1,ip),q(1,ip),c(1,ip),d(1,ip),e(1,ip),
     &        w3(1,ip),
     &        wallvr,wallvi,wbci,
     &        wallur,wallui,wallwr,wallwi,v_wall,wlvr,wlvi,
     &        iles,gur,gui,taur,taui,gu3r(1,ip),gu3i(1,ip),
     &        ggu3r(1,ip),ggu3i(1,ip),
     &        gthr,gthi,gsr,gsi,gth3r(1,ip),gth3i(1,ip),
     &        ggth3r(1,ip),ggth3i(1,ip),
     &        gewy1,gewy2,gewy12,gewy22,iord,
     &        diags,diags2,gplane(1,ip),filtxz,filtxz2,
     &        lpfxz,lpfxz2,zbp,ihighorder,cs,my_node,tbc,
     &        dtheta0_upp,dtheta0_low,d2th_bc(1,ip),fth_bc(1,ip),
     &        dth3r_bc(1,ip),dth3i_bc(1,ip),
     &        dth3r_bc2(1,ip),dth3i_bc2(1,ip),theta0_low,theta0_upp,
     &        cflux,mflux,vsuc)
      end do
c-------------------------------------

      call ctim(ctime4,wtime4)
      ctime_linear = ctime_linear + (ctime4-ctime3)
      wtime_linear = wtime_linear + (wtime4-wtime3)
      ctime3=ctime4
      wtime3=wtime4

c---- RK SUB-SEP ITERATION FINISHED

      it=it+1

      if (mod(it-1,nst).eq.0) then
c
c     If at the end of a full RK step
c      write(ios,*) '>>>>>> FULL STEP <<<<<<'
c
c     Update time step and shift of lower wall
c
         t = t + dtn
         xs = xs - dtn*u0low
         zs = zs - dtn*w0low
         if (vart) dtn = dtnp1
c
c     Save velocity field?
c
         if (t+.5*dtn.gt.tsave(isave)) then
c
c     Check if file already exists
c
            open(unit=11,file=nmsave(isave),form='unformatted',
     &           status='old',iostat=kode)
#ifdef MPI
            if (nproc.gt.1) call mpi_barrier(mpi_comm_world,ierror)
#endif
c
c     If you want to replace files: uncomment next line
c
            kode = 1
            if (kode.eq.0) then
c
c     File already exists: stop
c
               close(11)
               write(ioe,*) 'Velocity field already exists: ',
     &              trim(nmsave(isave))
               call stopnow(57584)
            else
c
c     File does not exist: continue
c
               if (my_node.eq.0) then
                  write(ios,*) '** Write velocity field **',
     &                 t/dstar,isave
               end if

               call wdiscbl(ur,ui,re,pr,m1,xl,zl,t,xs,dstar,
     &              fltype,bstart,
     &              blength,rlam,spanv,nmsave(isave),3+scalar,gall,
     &              boxr,boxi,urx,alfa,zs,beta,my_node,gr)

            end if
            isave=isave+1
         end if

         if (my_node.eq.0) then
c
c     Writing out internal buffer for logfile
c
c            call cflush(6)
c
c     Writing out internal buffers for amplitude and plane files
c
            if (iamp.gt.0) call cflush(15)

            if(j_en_flag) then 
               call cflush(42)
            end if

            if(j_un_chk) then 
               do ybp = 1,n_un
                  yb = 48 + ybp
c                  call cflush(yb)
               end do
            end if   

            do i=1,npl
               call cflush(30+i)
            end do
         end if
      end if

      call ctim(ctime4,wtime4)
      ctime_io = ctime_io + (ctime4-ctime3)
      wtime_io = wtime_io + (wtime4-wtime3)
      ctime3=ctime4
      wtime3=wtime4

      if (mod(it-1,nst).eq.0) then
c
c     Start timing for next step...
c     
         if (my_node.eq.0) then
            call ctim(ctime2,wtime2)
            write(ios,'(a,f10.3,a,f10.3,a)') ' CPU time: ',
     &           ctime2-ctime1,' sec    Wall time: ',
     &           wtime2-wtime1,' sec'
            write(ios,*) 'Simulation time',t/dstar,
     &           'it: ',it-1
c     
c     Write out step-file
c     writing is done (approximately) every delta_wall_time_int
c     seconds. If positive, also history is written.
c     
            delta_wall_time = delta_wall_time + wtime2-wtime1
            
            if (delta_wall_time.gt.abs(delta_wall_time_int)) then
               delta_wall_time = 0.
               open(file='step.out',unit=87)
               write(87,'(a20,f18.9,a,f18.9,a)') 'Time now: ',t/dstar,
     &              ' (max. ',tmax/dstar,')'
               write(87,'(a20,f18.9)') 'Started with: ',sttimet/dstar
               write(87,'(a20,f18.9)') 'Delta: ',(t-sttimet)/dstar
               write(87,'(a20,i18,a,i18,a)') 'Iteration: ',it-1,
     &              ' (max. ',maxit,')'
               call time_string(timdat)
               write(87,'(a20,a20)') 'Time now: ',timdat
               write(87,'(a20,a20)') 'Started: ',sttimec
               write(87,'(a,f18.9,a,f18.9,a)') 'CPU time       ',
     &              ctime2-ctime0,
     &              ' seconds (max. ',cpumax,')'
               write(87,'(a,f18.9,a,f18.9,a)') 'Wallclock time ',
     &           wtime2-wtime0,
     &              ' seconds (max. ',wallmax,')'
               write(87,'(a,f18.9,a)') 'Last time step ',ctime2-ctime1,
     &              ' seconds (CPU time)'
               write(87,'(a,f18.9,a)') 'Last time step ',wtime2-wtime1,
     &              ' seconds (wallclock time)'
               close(87)
c
c     History file
c
c               if (delta_wall_time_int.gt.0) then
c                  open(file='history.out',position='append',unit=87)
c                  write(87,'(f18.9,a25,f30.3)') t/dstar,timdat,wtime2
c                  close(87)
c               end if
            end if
            
            ctime1=ctime2
            wtime1=wtime2
c
c     Decide iteration loop break criteria (only on node 0)
c
            stop_now = 0
c
c     For the jet in crossflow:
c     Terminate bla only when the statistics are computed
c
c            if (wbci.eq.-2) then
c               if (ixys.ne.0.and.mod(it-1,ixys).eq.0)  then
c            if (wbci.eq.-2 .and. ixys.gt.0) then


            if (ixys.gt.0.and.t.ge.txys.and.mod(it-1,ixys).ne.0) then
            else

c            if (mod(it-1,ixys).eq.0)  then

             ! Check stop.now file
               open(unit=11,file='stop.now',status='old',iostat=kode)
               if (kode.eq.0) then
                  close(11)
                  stop_now = 1
                  write(ios,*) 'STOP DUE TO stop.now FILE'
               end if

            ! Check wall time
               if (wallmax.gt.0) then
                  if ((wtime2-wtime0).gt.wallmax) then
                     stop_now = 1
                     write(ios,*) 'STOP DUE TO WALL TIME'
                  end if
               end if

            ! Check cpu time
               if (cpumax.gt.0) then
                  if ((ctime2-ctime0).gt.cpumax) then
                     stop_now = 1
                     write(ios,*) 'STOP DUE TO CPU TIME'
                  end if
               end if

            ! Check iterations
               if (maxit.gt.0) then
                  if (it-1.ge.maxit) then
                     stop_now = 1
                     write(ios,*) 'STOP DUE TO ITERATIONS'
                  end if
               end if

            ! Check integration time
               if(tmax.gt.0) then
                  if (t+.5*dtn.gt.tmax) then
                     stop_now = 1
                     write(ios,*) 'STOP DUE TO INTEGRATION TIME'
                  end if
               end if
            end if
         
         end if
#ifdef MPI
c
c     Communicate break criterion
c
         if (nproc.gt.1) then
            call mpi_bcast(stop_now,1,mpi_integer4,0,mpi_comm_world,
     &           ierror)
         end if
#endif

      end if
c
c     Continue integration for stop_now=0
c
      if (stop_now.eq.0) then
         goto 1
      end if

c------------------------------------
c>>>>>>   END OF INTEGRATION   <<<<<<
c------------------------------------

      if (my_node.eq.0) then
         write(ios,*) '** Stop iteration **',t/dstar
      end if

      call ctim(ctime4,wtime4)
      ctime_rest = ctime_rest + (ctime4-ctime3)
      wtime_rest = wtime_rest + (wtime4-wtime3)
      ctime3=ctime4
      wtime3=wtime4
c
c     Do all the work after time integration
c     (note that one could do some statistics)
c
c     Write statistics
c
      if (ixys.gt.0 .and.t.ge.txys.and.sumw.gt.0.) then
         if (my_node.eq.0) then
            write(ios,*) '** Write statistics **', t/dstar,sumw/dstar
         end if
         call wxys(xys,xysth,
     &        totcorr,totcorr_x,corr,corr_x,ncorr,ncorr_x,corrf,corrf_x,
     &        totcorr_ms,totcorr_ms_x,corr_ms,corr_ms_x,my_node,
     &        series,serf,nser,totseries,
     &        sumw,re,xl,zl,t,dstar,fltype,
     &        bstart,blength,rlam,spanv,namxys,wxy,
     &        corrnam,corrnam_x,corrxx,corryy,corrzz,corryy_x,corrwsave,
     &        corrwsave_x,corrx,corry,corrz,corry_x,corrxi,corryi,
     &        corrzi,corryi_x,pr,m1,corrt,corrt_x,nsamp,namser,sercount,
     &        mhd_n,b0,it,ixyss)
         sercount = 0
      end if

      if(t_av) then
c
c     Jose: Write out the time averaged mean velocity into a file
c
         call w_mean_u(tav_u,my_node,eta,tav_cnt)
c
      end if
c
c
c
      if (arnoldi.eq.1) then
         call bla2arpck(arpckut,ur,ui,maxn)
      else
         if (wrapper.eq.0) then
c
c     Write end velocity field
c
            if (my_node.eq.0) then
               write(ios,*) '** Write field **', t/dstar
            end if

            call wdiscbl(ur,ui,re,pr,m1,xl,zl,t,xs,dstar,fltype,
     &           bstart,blength,rlam,spanv,namnut,3+scalar,gall,
     &           boxr,boxi,urx,alfa,zs,beta,my_node,gr)
         end if
      end if
c
c     Close open files
c
      if (my_node.eq.0) then
         if (iamp.gt.0) close(unit=15)
         if (mwave.gt.0) close(unit=16)
         if (iext.gt.0) close(unit=18)
         do i=1,npl
            close(unit=30+i)
         end do
      end if
c
c     Jose: Closing files used for saving "entropy" calculations. 
c
      if(my_node.eq.0) then 
         if(j_en_flag) then
            close(unit=42)
         end if
         if(j_un_chk) then
            do ybp = 1,n_un
               yb = 48 + ybp
c               close(unit=yb)
            end do
         end if
      end if

c
c     Write summary
c
      if (my_node.eq.0) then
         call ctim(ctime4,wtime4)
         ctime_end = ctime4-ctime3
         wtime_end = wtime4-wtime3

#ifdef OPENMP
         ompnproc=OMP_GET_NUM_PROCS()
         ompmaxthreads=OMP_GET_MAX_THREADS()
#else
         ompnproc = 1
         ompmaxthreads = 1
#endif

         write(ios,*)
         write(ios,*)
         write(ios,*) '>>>>>>  Summary  <<<<<<'
         write(ios,*) '--------------------------------------------'//
     &              '-----------------------'
         write(ios,*)
         write(ios,'(a23,f18.9)') ' Time now             :',t/dstar
         write(ios,'(a23,f18.9)') ' Started with         :',
     &        sttimet/dstar
         write(ios,'(a23,f18.9)') ' Delta                :',
     &                          (t-sttimet)/dstar
         write(ios,'(a23,i6)') ' Iterations           :',it-1
         call time_string(timdat)
         write(ios,'(a24,a20)') ' Started              : ',sttimec
         write(ios,'(a24,a20)') ' Ended                : ',timdat
         write(ios,*)
         write(ios,*) 'Total CPU time       :',ctime4-ctime0,' seconds'
         write(ios,*) 'Total wallclock time :',wtime4-wtime0,' seconds'
         write(ios,*) 'Average time step    :',
     &              (wtime4-wtime0-wtime_startup-wtime_end)
     &        /real(it-1)*nst
         write(ios,*)

         write(ios,'(a,i4,a)') ' Compiled for         :',nproc,
     &        ' processors'
         write(ios,'(a,i4,a)') '                      :',nthread,
     &        ' threads'
         write(ios,'(a,i4,a)') ' Running on           :',ompnproc,
     &              ' OMP processors'
         write(ios,'(a,i4,a)') ' Running on           :',all_nodes,
     &              ' MPI processors'
         write(ios,'(a,i4,a)') ' Running with max.    :',ompmaxthreads,
     &              ' OMP threads'
         write(ios,*)
         write(ios,'(a,3a20)')   '             ','Wall time','CPU time',
     &        'Wall time per it.'
         write(ios,'(a,3f20.8)') ' LES        : ',
     &        wtime_les,ctime_les,wtime_les/real(it-1)*nst
         write(ios,'(a,3f20.8)') ' tau_ij     : ',
     &        wtime_tau,ctime_tau,wtime_tau/real(it-1)*nst
         write(ios,'(a,3f20.8)') ' statistics : ',
     &        wtime_stat,ctime_stat,wtime_stat/real(it-1)*nst
         write(ios,'(a,3f20.8)') ' nonlin     : ',
     &        wtime_nonlin,ctime_nonlin,wtime_nonlin/real(it-1)*nst
         write(ios,'(a,3f20.8)') ' particles  : ',
     &        wtime_part,ctime_part,wtime_part/real(it-1)*nst
         write(ios,'(a,3f20.8)') ' linear     : ',
     &        wtime_linear,ctime_linear,wtime_linear/real(it-1)*nst
         write(ios,'(a,3f20.8)') ' IO         : ',
     &        wtime_io,ctime_io,wtime_io/real(it-1)*nst
         write(ios,'(a,3f20.8)') ' trip       : ',
     &        wtime_trip,ctime_trip,wtime_trip/real(it-1)*nst
         write(ios,'(a,3f20.8)') ' ffun+mhd   : ',
     &        wtime_ffun,ctime_ffun,wtime_ffun/real(it-1)*nst
         write(ios,'(a,3f20.8)') ' CFL        : ',
     &        wtime_cfl,ctime_cfl,wtime_cfl/real(it-1)*nst
         write(ios,'(a,3f20.8)') ' planes     : ',
     &        wtime_planes,ctime_planes,wtime_planes/real(it-1)*nst
         write(ios,'(a,3f20.8)') ' press      : ',
     &        wtime_press,ctime_press,wtime_press/real(it-1)*nst
         write(ios,'(a,3f20.8)') ' rest       : ',
     &        wtime_rest,ctime_rest,wtime_rest/real(it-1)*nst
         write(ios,'(a,2f20.8)') ' startup    : ',
     &        wtime_startup,ctime_startup
         write(ios,'(a,2f20.8)') ' shutdown   : ',
     &        wtime_end,ctime_end
         write(ios,'(a,3f20.8)') ' sum        : ',
     &        wtime_stat+wtime_nonlin+wtime_part+wtime_linear+wtime_io+
     &        wtime_startup+wtime_end+wtime_les+wtime_tau+
     &        wtime_rest+wtime_trip+wtime_ffun+wtime_cfl+
     &        wtime_planes+wtime_press,
     &        ctime_stat+ctime_nonlin+ctime_part+ctime_linear+ctime_io+
     &        ctime_startup+ctime_end+ctime_les+ctime_tau+
     &        ctime_rest+ctime_trip+ctime_ffun+ctime_cfl+
     &        ctime_planes+ctime_press,
     &        (wtime_stat+wtime_nonlin+wtime_part+wtime_linear+wtime_io+
     &        wtime_les+wtime_tau+
     &        wtime_rest+wtime_trip+wtime_ffun+wtime_cfl+
     &        wtime_planes+wtime_press)/real(it-1)*nst
         write(ios,'(a,3f20.8)') ' nl(com,ser): ',tcom,tser,tcom+tser
c
c     Write final step file
c
         open(file='step.out',unit=87)
         write(87,*) '** Finished **'
         write(87,'(a22,f18.9)') ' Time now             : ',t/dstar
         write(87,'(a22,f18.9)') ' Started with         : ',
     &                           sttimet/dstar
         write(87,'(a22,f18.9)') ' Delta                : ',
     &                           (t-sttimet)/dstar
         write(87,'(a23,i6)') ' Iterations           : ',it-1
         write(87,'(a24,a20)') ' Started              : ',sttimec
         write(87,'(a24,a20)') ' Ended                : ',timdat
         write(87,*)
         write(87,*) 'Total CPU time       :',ctime4-ctime0,' seconds'
         write(87,*) 'Total wallclock time :',wtime4-wtime0,' seconds'
         write(87,*) 'Average time step    :',
     &               (wtime4-wtime0-wtime_startup-wtime_end)
     &        /real(it-1)*nst
         write(87,*)
         write(87,'(a,i4,a)') ' Compiled for         :',nproc,
     &        ' processors'
         write(87,'(a,i4,a)') '                      :',nthread,
     &        ' threads'
         write(87,'(a,i4,a)') ' Running on           :',ompnproc,
     &              ' OMP processors'
         write(87,'(a,i4,a)') ' Running on           :',all_nodes,
     &              ' MPI processors'
         write(87,'(a,i4,a)') ' Running with max.    :',ompmaxthreads,
     &              ' OMP threads'
         write(87,*)
         write(87,'(a,3a20)')   '             ','wall time','CPU time',
     &        'Wall time per it.'
         write(87,'(a,3f20.8)') ' LES        : ',
     &        wtime_les,ctime_les,wtime_les/real(it-1)*nst
         write(87,'(a,3f20.8)') ' tau_ij     : ',
     &        wtime_tau,ctime_tau,wtime_tau/real(it-1)*nst
         write(87,'(a,3f20.8)') ' statistics : ',
     &        wtime_stat,ctime_stat,wtime_stat/real(it-1)*nst
         write(87,'(a,3f20.8)') ' nonlin     : ',
     &        wtime_nonlin,ctime_nonlin,wtime_nonlin/real(it-1)*nst
         write(87,'(a,3f20.8)') ' particles  : ',
     &        wtime_part,ctime_part,wtime_part/real(it-1)*nst
         write(87,'(a,3f20.8)') ' linear     : ',
     &        wtime_linear,ctime_linear,wtime_linear/real(it-1)*nst
         write(87,'(a,3f20.8)') ' IO         : ',
     &        wtime_io,ctime_io,wtime_io/real(it-1)*nst
         write(87,'(a,3f20.8)') ' trip       : ',
     &        wtime_trip,ctime_trip,wtime_trip/real(it-1)*nst
         write(87,'(a,3f20.8)') ' ffun+mhd   : ',
     &        wtime_ffun,ctime_ffun,wtime_ffun/real(it-1)*nst
         write(87,'(a,3f20.8)') ' CFL        : ',
     &        wtime_cfl,ctime_cfl,wtime_cfl/real(it-1)*nst
         write(87,'(a,3f20.8)') ' planes     : ',
     &        wtime_planes,ctime_planes,wtime_planes/real(it-1)*nst
         write(87,'(a,3f20.8)') ' press      : ',
     &        wtime_press,ctime_press,wtime_press/real(it-1)*nst
         write(87,'(a,3f20.8)') ' rest       : ',
     &        wtime_rest,ctime_rest,wtime_rest/real(it-1)*nst
         write(87,'(a,2f20.8)') ' startup    : ',
     &        wtime_startup,ctime_startup
         write(87,'(a,2f20.8)') ' shutdown   : ',
     &        wtime_end,ctime_end
         write(87,'(a,3f20.8)') ' sum        : ',
     &        wtime_stat+wtime_nonlin+wtime_part+wtime_linear+wtime_io+
     &        wtime_startup+wtime_end+wtime_les+wtime_tau+
     &        wtime_rest+wtime_trip+wtime_ffun+wtime_cfl+
     &        wtime_planes+wtime_press,
     &        ctime_stat+ctime_nonlin+ctime_part+ctime_linear+ctime_io+
     &        ctime_startup+ctime_end+ctime_les+ctime_tau+
     &        ctime_rest+ctime_trip+ctime_ffun+ctime_cfl+
     &        ctime_planes+ctime_press,
     &        (wtime_stat+wtime_nonlin+wtime_part+wtime_linear+wtime_io+
     &        wtime_les+wtime_tau+
     &        wtime_rest+wtime_trip+wtime_ffun+wtime_cfl+
     &        wtime_planes+wtime_press)/real(it-1)*nst
         write(87,'(a,3f20.8)') ' nl(com,ser): ',tcom,tser,tcom+tser

         close(87)

         write(ios,*)
         write(ios,*) 'bla terminated normally.'
      end if
c
c     Finalize program
c
#ifdef MPI
      call mpi_barrier(mpi_comm_world,ierror)
      if(wrapper.eq.0) call mpi_finalize(ierror)
#endif

#ifdef WRAPPER
      end subroutine bla
#else
      end program bla
#endif
