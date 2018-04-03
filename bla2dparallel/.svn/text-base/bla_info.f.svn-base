

      program bla_info

      implicit none

      character*80 namnin
      integer,parameter :: scalar=0
      real re,xl,zl,t,pr(scalar),m1(scalar),dstar,rlam,spanv
      real bstart,bslope,xs
      integer i
      integer nx,nyp,nzc,nfzsym
      logical pou
      integer fltype


      call getarg(1,namnin)

      write(*,*) 'file = >>',trim(namnin),'<<'


      open(unit=11,file=trim(namnin),form='unformatted',status='old')
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
      close(unit=11)

      write(*,*) 'fltype = ',fltype
      write(*,*) 'scalar = ',scalar
      write(*,*) 'nx,nyp,nzc = ',nx,nyp,nzc
      write(*,*) 're = ',re*dstar
      write(*,'(a,3f18.9)') ' xl,yl,zl = ',xl/dstar,2./dstar,zl/dstar
      write(*,*) 't = ',t/dstar





      end program bla_info
