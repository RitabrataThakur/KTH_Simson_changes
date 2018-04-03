c************************************************************************
c     To check how basic fft transforms work in x-z plane
c     
c ***********************************************************************

      subroutine fft_test(gridx,gridz,prexn,prezn,prex,prez,alfa,beta,
     &     xl,zl)

      implicit none
      
      include 'par.f'

      
      real alfa(nx/2,mbz),beta(nz)
      real prex(nxp+15),prey(nyp*2+15)
      real prez(nzp*2+15)
      real prexn(nx+15),prezn(nz*2+15)
      real gridx(nx),gridz(nz)

      real rand1r((nxp/2+1),nzd),rand1i((nxp/2+1),nzd)
      real rand2r((nxp/2+1),nzd),rand2i((nxp/2+1),nzd)
      real wr((nxp/2+1),nzd), wi((nxp/2+1),nzd)
      real tmp1, tmp2, tmp3

c
c     As per the documentation, 
c     rand1r are the odd points and rand1i are even in physical space
c     In spectral space, real and imaginary. 
c

      integer x,z,xx,zz
      integer len1, len2
      real pi
      real xl,zl
      parameter (pi = 3.1415926535897932385)

      rand1r = 0.
      rand1i = 0.

c      do x = 1,nx/2
c         write(20,*) x, alfa(x,1)
c      end do
      
c      do z = 1,nz
c         write(19,*) z, beta(z)
c      end do
c
c     Define f = sin(x)*cos(z) over nx physical points
c
      do z = 1,nz
         do x = 1,nx/2
            xx = 2*x - 1
            rand1r(x,z) =sin(gridx(xx))*cos(gridz(z))
            xx = 2*x
            rand1i(x,z) = sin(gridx(xx))*cos(gridz(z))
c            write(21,*) x, z, rand1r(x,z), rand1i(x,z)
         end do
      end do

      rand2r = rand1r
      rand2i = rand1i
c
c     Forward Fourier transform in x
c     Real sequences to be transformed here. 
      call vrfftf(rand1r,rand1i,wr,wi,nx,nz,1,nxp/2+1,prexn)

c
c     Backward Fourier transform in x
c
c      call vrfftb(rand1r,rand1i,wr,wi,nx,nz,1,nxp/2+1,prexn)
      
c      write(*,*) "Output of the forward transform in x" 
c      do z = 1,nz
c        do x = 1,nx/2
c            write(21,*) x, z, rand1r(x,z), rand2r(x,z)
c         end do
c      end do

c     Works till here. 

c
c     Forward Fourier transform in z
c     Complex sequences 

      call vcfftf(rand1r,rand1i,wr,wi,nz,nx/2,nxp/2+1,
     &     1,prezn)


c      write(*,*) "Output of the forward transform in z" 
c      do z = 1,nzd
c         do x = 1,nxp/2+1
c            write(23,*) x, z, rand1r(x,z), rand1i(x,z)
c         end do
c      end do
c     Here, notice that nx/2 + 1 sequences have to be actually transformed.
c     However wavenumbers nx/2, nz/2 are oddball wavenumbers.
c     So on chosing to transform only nx/2 sequences, differences arise only
c     at the locations corresponding to nx/2 and nz/2 wavenumbers. 
c     An interesting observation perhaps is that the nx/2 wavenumber amplitudes
c     seem to be complex when nx/2 + 1 sequences are chosen.
c     Other important aspect: if nx/2 is chosen, the nx/2 + 1 row is the same
c     as the output of the x Fourier transform.

c
c     Remove oddball in z
c
      z = nz/2 + 1
      do x = 1,nx/2
         rand1r(x,z) = 0. 
         rand1i(x,z) = 0.
      end do
c
c     Remove oddball in x
c
      x = nx/2 + 1
      do z = 1,nz
         rand1r(x,z) = 0.
      end do

c
c     Normalisation
c      
      rand1r = rand1r/real(nx*nz)
      rand1i = rand1i/real(nx*nz)
c
c     Backward transforms: an extended grid
c

c
c     1. Shift the spectrum
c
      do z =nz/2+1,nz
         zz = nzp - nz + z
         do x = 1,nx/2
            rand1r(x,zz) = rand1r(x,z)
            rand1r(x,z) = 0.
            rand1i(x,zz) = rand1i(x,z)
            rand1i(x,z) = 0.
         end do
      end do

c
c     Backward Fourier transforms on extended grid
c
      call vcfftb(rand1r,rand1i,wr,wi,nzp,nx/2,nxp/2+1,1,prez)
      call vrfftb(rand1r,rand1i,wr,wi,nxp,nzp,1,nxp/2+1,prex)


      do z = 4,4
         tmp2 = gridz(z)
         tmp3 = cos(tmp2)
         do x = 1,nx/2
            xx = 2*x - 1
            write(21,*) gridx(xx), rand2r(x,z)/tmp3
            xx = 2*x
            write(21,*) gridx(xx), rand2i(x,z)/tmp3
         end do
      end do

      do z = 4,4
         tmp2 = zl/real(nzp)*real(z-1)
         tmp3 = cos(tmp2)
         do x = 1,nxp/2
            xx = 2*x - 1
            tmp1 = xl/real(nxp)*real(xx-1)
            write(22,*) tmp1, rand1r(x,z)/tmp3
            xx = 2*x
            tmp1 = xl/real(nxp)*real(xx-1)
            write(22,*) tmp1, rand1i(x,z)/tmp3
         end do
      end do

c
c     For defining a nonlinear term.
c

c      write(*,*) 'To see what is stored in the final spot'
c      write(*,*) rand1r(nxp/2+1,z),  rand1i(nxp/2+1,z)

c
c     Comparing the results before and after FFT
c
c      do z = 1, nz
c         do x = 1,nx/2
c            write(*,*) x, z
c            write(25,*) x,z,rand1r(x,z),rand2r(x,z)
c            write(26,*) x,z,rand1i(x,z),rand2i(x,z)
c rand1i(x,z),
c     &           rand2i(x,z)
            
c         end do
c      end do




      end subroutine
