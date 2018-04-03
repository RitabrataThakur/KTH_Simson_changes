c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/par.f $
c $LastChangedDate: 2011-04-14 10:38:01 +0200 (Thu, 14 Apr 2011) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 1659 $
c
c ***********************************************************************
c
c     par.f contains size of problem
c
      integer nx,ny,nz,nby,nbz,mby,mbz
      integer nfxd,nfyd,nfzd,nfzsym
      integer nxp,nyp,nzp,nzc,nzat,nzst,nzd,nzpc
      integer memnx,memny,memnz,memnxyz
      integer mbox2,mbox3,nproc,nthread,nxys,nxysth
      integer osnf,mcorr,mbla,mser,msamp,scalar,pressure
      integer ios,ioe
      integer npart,nppart
c     Jose variables
c      real delT, ct, Gras, Pran
      real ct
c
c     Adjustable parameters:
c     ======================
c
c     Number of spectral modes

      parameter (nx=32,ny=81,nz=32)
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
      parameter (nxysth=  44)
      parameter (mcorr =  30)
      parameter (mser  =  25)
      parameter (msamp =  32)
c
c     Pressure (0/1)
c
      parameter (pressure=0)
c
c     Passive scalars
c
      parameter (scalar=1)
c
c     Particles
c
      parameter (npart=0)
c
c     Dealiasing flags
c
      parameter (nfxd=1,nfyd=0,nfzd=1)
c
c     Number of waves for freestream OS-eigenmodes
c
      parameter (osnf = 0)
c
c     Number of points in base flow
c
      parameter (mbla = 20001)
c
c     Printing unit ios (output unit) and ioe (error unit)
c
      parameter (ios=6,ioe=6)
c
c     Temporally disabled parameters: (don't change them!)
c     ====================================================
c
c     Symmetry flag
c
      parameter (nfzsym=0)
c
c     Boxsize
c
      parameter (mby=1,mbz=1)
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
      parameter (memnz=nzc/nproc)
c
c     The next line has to be used if the mpi_alltoall implementation
c     in getpxz and putpxz is employed, i.e. ALLTOALL is set.
c      parameter (memny=(nyp/nproc+1)*nproc)
c
c     Otherwise, use the next line, i.e. ALLTOALL is not set.
      parameter (memny=nyp)
c
c     Arrangement of main storage:
c
c     1-3:   velocities
c     4-5:   vorticities
c     6-7:   partial right-hand sides
c     8  :   pressure (if compiled accordingly)
c     9-11:  theta, dtheta/dy, prhs (scalar 1, if compiled accordingly)
c     12-:   repeated entries 9-11 for additional scalars
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
c     For viscosity dependent on the temperature. See Sameen, Govindarajan (2007)
c     delT = T_up - T_low
      parameter(ct = 1879.9)
c ct: parameter constant for defining viscosity (Arrhenius model)
c      parameter(delT = 0)
c Gras: Richardson number 
c      parameter(Gras = 0.01)
c Pran: Prandtl number
c      parameter(Pran = 1)
c ***********************************************************************
