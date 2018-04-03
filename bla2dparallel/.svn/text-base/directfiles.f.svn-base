c     To compile on BG/L
c     /opt/ibmcmp/xlf/bg/11.1/bin/xlf -qextname -qrealsize=8 -qzerosize -o multfiles multfiles.f


      program directfiles

      implicit none

      character*80 namnin
      character*80 namnut
      character(len=5) ch
      integer,parameter :: scalar=0
      real re,xl,zl,t,pr(scalar),m1(scalar),gr(scalar)
      real dstar,rlam,spanv
      real bstart,bslope,xs
      integer i,j,n,m
      integer nx,nyp,nzc,nfzsym
      logical pou,phys,ini_bytes
      character(len=4) vars
      integer fltype,prec
      integer nproc1,nprocz1,memnx1,memny1,memnz1
      integer nproc2,nprocz2,memnx2,memny2,memnz2
      integer x,y,z,k,reclen,irec,ini_rec

      real,allocatable :: ur(:,:,:,:,:),ui(:,:,:,:,:)
      real,allocatable :: urx(:)



      call getarg(1,namnin)
      namnut = namnin(1:len(trim(namnin))-2)//'v'

      write(*,*) trim(namnin),'-->',trim(namnut)

      nproc2 = 1

      open(unit=11,file=trim(namnin)//'-00000',form='unformatted')
      rewind(11)
      
      if (scalar.ge.1) then
         read(11) re,pou,xl,zl,t,xs,(pr(i),m1(i),i=1,scalar)
      else
         read(11) re,pou,xl,zl,t,xs
      end if
      read(11) nx,nyp,nzc,nfzsym
      read(11) fltype,dstar
      if (fltype.eq.-1) read(11) rlam
      if (fltype.eq.-2) read(11) rlam,spanv
      if (fltype.ge.4) read(11) bstart,bslope,rlam,spanv
      read(11) nproc1
      close(unit=11)

      write(*,*) 'Reading ',nproc1,' files'
      write(*,*) 'fltype = ',fltype
      write(*,*) 'scalar = ',scalar
      write(*,*) 'nx,nyp,nzc = ',nx,nyp,nzc
      write(*,*) 't = ',t/dstar


      nprocz1 = (nproc1)**0.5
      memnx1  = nx/(2*nprocz1)
      memny1  = (nyp/nprocz1+1)*nprocz1
      memnz1  = nzc/nprocz1

      nprocz2 = (nproc2)**0.5
      memnx2  = nx/(2*nprocz2)
      memny2  = (nyp/nprocz2+1)*nprocz2
      memnz2  = nzc/nprocz2

      write(*,'(5i5)') nproc1,nprocz1,memnx1,memny1,memnz1
      write(*,'(5i5)') nproc2,nprocz2,memnx2,memny2,memnz2

      allocate(ur(memnx1,memny1,memnz1,3+scalar,nprocz1))
      allocate(ui(memnx1,memny1,memnz1,3+scalar,nprocz1))
      allocate(urx(nx))

c     reserve a total of 4096 bytes in the beginning.
c
c     One record is nx*prec bytes, or 
c     nx*reclen in Fortran units.
c     


      phys = .false.
      prec = kind(ur(1,1,1,1,1))
      ini_bytes = ini_rec*nx*prec
      vars = 'XYZM'

      write(*,*) 'real accuracy: ',prec
      if (prec.ne.8) stop

      inquire(iolength=reclen) ur(1,1,1,1,1)
      write(*,*) 'reclen = ',reclen,' for ',prec,' bytes'




      ini_rec = (4095/(nx*prec)+1)

      write(*,*) 'initial records for header: ',ini_rec,
     &     ' bytes: ',ini_rec*nx*prec

      open(unit=100,file=trim(namnut),form='unformatted')
      rewind(100)

c      write(100) '***BEGIN HEADER***'
      write(100) '***IN PROGRESS!***'
      write(100) 'BLA VELOCITY FIELD VERSION 2.00'
      write(100) scalar
      if (scalar.ge.1) then
         write(100) re,pou,xl,zl,t,xs,(pr(i),m1(i),i=1,scalar)
      else
         write(100) re,pou,xl,zl,t,xs
      end if
      write(100) nx,nyp,nzc,nfzsym
      write(100) fltype,dstar
      if (fltype.eq.-1) write(100) rlam
      if (fltype.eq.-2) write(100) rlam,spanv
      if (fltype.eq.4) write(100) bstart,bslope
      if (fltype.eq.5) write(100) bstart,bslope
      if (fltype.ge.6.and.fltype.le.9) 
     &     write(100) bstart,bslope,rlam,spanv
      if (abs(fltype).eq.20) write(100) (gr(i),i=1,scalar)

      write(100) prec,ini_bytes
      write(100) phys,vars
      write(100) '***END HEADER***'
      close(100)




      open(unit=100,file=namnut,access='direct',recl=reclen*nx)


c     write last record
      urx = 0.
      write(100,rec=ini_rec+(3+scalar)*nyp*nzc) urx

      do j=1,nprocz1

         do i=1,nprocz1
            n = i+(j-1)*nprocz1
            write(ch,'(i5.5)') n
            open(unit=11,file=trim(namnin)//'-'//ch,form='unformatted')
            rewind(11)
            write(*,*) 'reading : ',i,j,n,trim(namnin)//'-'//ch
            do m=1,3+scalar
               read(11) ur(:,:,:,m,i),ui(:,:,:,m,i)
            end do
            close(11)
         end do

         do m=1,3+scalar
            write(*,*) 'write m=',m
            do z=1,memnz1
               do y=1,nyp                  
                  k=1
                  do i=1,nprocz1
                     do x=1,memnx1
                        urx(k) = ur(x,y,z,m,i)
                        k=k+1
                        urx(k) = ui(x,y,z,m,i)
                        k=k+1
                     end do
                  end do
                  irec = ini_rec+
     &                 (m-1)*nyp*nzc +
     &                 (j-1)*nyp*memnz1+
     &                 (z-1)*nyp +  
     &                 (y-1)+1
                  write(100,rec=irec) urx
               end do
            end do
         end do
      end do

      close(100)


      open(unit=100,file=trim(namnin)//'.done',form='unformatted')
      close(100)


      end program directfiles
