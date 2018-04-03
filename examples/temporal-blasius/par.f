c ***********************************************************************
C
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/examples/temporal-blasius/par.f $
c $LastChangedDate: 2010-10-21 13:39:19 +0200 (Thu, 21 Oct 2010) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 1594 $
c
c ***********************************************************************
c
c     par.f contains size of problem
c
      integer nx,ny,nz
      integer nfxd,nfyd,nfzd,nfzsym
      integer nxp,nyp,nzp,nzc,nzat,nzst,nzd,nzpc
      integer nby,nbz,mby,mbz
      integer memnx,memny,memnz,memnxyz
      integer mbox2,mbox3
      integer nproc
      integer nthread
      integer nxys,nxysth
      integer osnf
      integer mcorr
      integer scalar
      integer pressure
      integer mbla
      integer mser,msamp
      integer ios,ioe
c
c     Adjustable parameters:
c     ======================
c
c     Number of spectral modes
c
      parameter (nx=80,ny=65,nz=80)
c
c     Dealiasing flags
c
      parameter (nfxd=1,nfyd=0,nfzd=1)
c
c     Symmetry flag
c
      parameter (nfzsym=0)
c
c     Boxsize
c
      parameter (mby=1,mbz=1)
c
c     Number of processors (MPI)
c
      parameter (nproc=1)
c
c     Number of threads (OpenMP)
c
      parameter (nthread=1)
c
c     Statistics
c
      parameter (nxys  =  96)
      parameter (nxysth=  35)
      parameter (mcorr =  30)
      parameter (mser  =  20)
      parameter (msamp =  16)
c
c     Pressure (0/1)
c
      parameter (pressure=1)
c
c     Passive scalar (0/1)
c
      parameter (scalar=0)
c
c     Number of waves for freestream OS-eigenmodes
c
      parameter (osnf=5)
c
c     Number of points in base flow
c
      parameter (mbla=20001)
c
c     Printing unit ios (output unit) and ioe (error unit)
c
      parameter (ios=6,ioe=6)
c
c     Computed parameters: (don't change them!)
c     =========================================
c
      parameter (nxp=nx+nfxd*nx/2,nyp=ny+nfyd*ny/2,nzp=nz+nfzd*nz/2)
      parameter (nzc=nz-nfzsym*nz/2)
      parameter (nzd=nzp-(nzp/2-2)*nfzsym)
      parameter (nzpc=nzp-(nzp/2-1)*nfzsym)
c
c     Symmetric transform lengths
c
      parameter (nzat=nzp/2-1,nzst=nzp/2+1)
c
c     Number of boxes
c
      parameter (nby=(nyp+mby-1)/mby,nbz=nzc/mbz)
c
c     Core storage size
c
      parameter (memnx=nx/2)
      parameter (memny=nyp)
      parameter (memnz=nzc/nproc)
c
c     Arrangement of main storage:
c
c     1-3:   velocities
c     4-5:   vorticities
c     6-7:   partial right-hand sides
c     8  :   pressure (if compiled accordingly)
c     9-11:  theta, dtheta/dy, prhs (if compiled accordingly)
c
      parameter (memnxyz=7+pressure+scalar*3)
c
c     Boxsize
c
      parameter (mbox2=(nxp/2+1)*nzd*mby*nthread)
      parameter (mbox3=nx/2*nyp*mbz*nthread)
c
c     Particles per processor
c      
      parameter (nppart=npart/nproc)
c
c     End of par.f
c
c ***********************************************************************
