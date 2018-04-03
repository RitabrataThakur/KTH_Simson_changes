c     To compile on BG/L
c     /opt/ibmcmp/xlf/bg/11.1/bin/xlf -qextname -qrealsize=8 -qzerosize -o multfiles multfiles.f
c     pgf90 -r8 -byteswapio -o multfiles multfiles.f
c     ifort -r8 -o multfiles multfiles.f

      program multfiles

      implicit none

      character*80 namnin
      character*80 namnut
      character pp
      character(len=5) ch
      integer scalar,scalar1
      integer  mf 
      real re,xl,zl,t,dstar,rlam,spanv
      real, allocatable :: pr(:),m1(:)
      real bstart,bslope,xs
      integer i,j,n,m
      integer nx,nyp,nzc,nfzsym
      logical pou
      integer fltype
      integer nproc1,nprocz1,memnx1,memny1,memnz1
      integer nproc2,nprocz2,memnx2,memny2,memnz2
      integer x,y,z,k

      real,allocatable :: ur(:,:,:,:,:),ui(:,:,:,:,:)
      real,allocatable :: urx(:)

c      namnin = 'end.uu'
c      namnin = 'field.01000.uu'
c      namnut = 'test.u'

c      write(*,*) 'give input name'
c      read(*,*) namnin
c      write(*,*) 'give output name'
c      read(*,*) namnut

      call getarg(1,namnin)
      call getarg(2,namnut)
c      namnut = namnin(1:len(trim(namnin))-1)

      write(*,*) 'Usage: multfiles infile.uu outfile.u [p|u]'
      write(*,*) '       with option p, the pressure fields are read'
      write(*,*) '       with option u, no scalars are included'
      write(*,*) '       without last option, all scalars are included'

      write(*,*) trim(namnin),'-->',trim(namnut)

      nproc2 = 1

      open(unit=11,file=trim(namnin)//'-00000',form='unformatted')
c
c     Determine number of scalars
c
      scalar = 0
 1111 continue
      rewind(11)
      allocate(pr(scalar),m1(scalar))
      read(11,err=1112) re,pou,xl,zl,t,xs,(pr(i),m1(i),i=1,scalar)
      scalar = scalar + 1
      deallocate(pr,m1)
      goto 1111

 1112 continue
      scalar = scalar-1
      deallocate(pr,m1)

c
c     Read file
c
      rewind(11)
      allocate(pr(scalar),m1(scalar))

      write(*,*) 'scalar = ',scalar

      mf = 3+scalar
      scalar1 = scalar

      call getarg(3,pp)
      if (pp(1:1).eq.'p') then
         write(*,*) 'Reading pressure field'
         mf = 1
      end if
      if (pp(1:1).eq.'u') then
         write(*,*) 'Reading only velocity field (no scalars)'
         mf = 3
         scalar1 = 0
      end if

      write(*,*) 'number of fields to process mf = ',mf

      read(11) re,pou,xl,zl,t,xs,(pr(i),m1(i),i=1,scalar)
      read(11) nx,nyp,nzc,nfzsym
      read(11) fltype,dstar
      if (fltype.eq.-1) read(11) rlam
      if (fltype.eq.-2) read(11) rlam,spanv
      if (fltype.ge.4) read(11) bstart,bslope,rlam,spanv
      read(11) nproc1
      close(unit=11)

      write(*,*) 'Reading ',nproc1,' files'
      write(*,*) 'fltype = ',fltype
      write(*,*) 'nx,nyp,nzc = ',nx,nyp,nzc
      write(*,*) 't = ',t

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

      allocate(ur(memnx1,memny1,memnz1,mf,nprocz1))
      allocate(ui(memnx1,memny1,memnz1,mf,nprocz1))
      allocate(urx(nx))


      do m=1,mf
         write(ch,'(i5.5)') m
         open(unit=100+m,file='dummy-'//ch,form='unformatted')
      end do


      do j=1,nprocz1

         do i=1,nprocz1
            n = i+(j-1)*nprocz1
            write(ch,'(i5.5)') n
            open(unit=11,file=trim(namnin)//'-'//ch,form='unformatted')
            rewind(11)
            write(*,*) 'reading : ',i,j,n,trim(namnin)//'-'//ch
            do m=1,mf
               read(11) ur(:,:,:,m,i),ui(:,:,:,m,i)
            end do
            close(11)
         end do

         do m=1,mf
            write(ch,'(i5.5)') m
            write(*,*) 'write : ',m,'dummy-'//ch
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
                  write(100+m) urx
               end do
            end do
         end do

      end do

      write(*,*)
      write(*,*) 'Write output to ',namnut
      open(unit=12,file=namnut,form='unformatted')
      rewind(12)
      if (scalar1.ge.1) then
         write(12) re,.false.,xl,zl,t,0.,(pr(i),m1(i),i=1,scalar1)
      else
         write(12) re,.false.,xl,zl,t,0.
      end if
      write(12) nx,nyp,nzc,nfzsym
      write(12) fltype,dstar
      if (fltype.eq.-1) write(12) rlam
      if (fltype.eq.-2) write(12) rlam,spanv
      if (fltype.ge.4) write(12) bstart,bslope,rlam,spanv

      do m=1,mf
         write(*,*) 'write field m=',m
         rewind(100+m)
         do z=1,nzc
            do y=1,nyp
               
               read(100+m) urx
               write(12) urx
               
            end do
         end do
         close(100+m)
      end do

      close(12)





      end program multfiles
