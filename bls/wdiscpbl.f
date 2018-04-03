c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/bls/wdiscpbl.f $
c $LastChangedDate: 2011-04-15 11:18:19 +0200 (Fri, 15 Apr 2011) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 1660 $
c
c ***********************************************************************
      subroutine wdiscpbl(ur,re,pr,gr,xl,zl,t,xs,dstar,fltype,
     &     bstart,blength,rlam,m1,spanv,namnut,m,pln,urx,ivarvisc)
c
c     Writes m variables form ur to file namnut
c
      implicit none

      include 'par.f'

      character*80 namnut
      integer m,fltype,ivarvisc
      real re,xl,zl,t,xs,dstar,bstart,blength,rlam,spanv
      real pr(scalar),m1(scalar),gr(scalar),ri(scalar)
      complex ur(memnx,memny,memnz,3+scalar)
      complex urx(nx/2)
      complex pln(nx/2,nyp)

      integer x,y,z,i,j

      open(unit=11,file=namnut,status='unknown',form='unformatted')

c
c     Write file header
c
      if (scalar.ge.1) then
         write(11) re,.false.,xl,zl,t,xs,(pr(j),gr(j),m1(j),j=1,scalar),
     &         ivarvisc
      else
         write(11) re,.false.,xl,zl,t,xs
      end if
      write(11) nx,nyp,nzc,nfzsym
      write(11) fltype,dstar
      if (fltype.eq.-1) write(11) rlam
      if (fltype.eq.-2) write(11) rlam,spanv
      if (fltype.ge.4.and.fltype.le.9)
     &     write(11) bstart,blength,rlam,spanv


c      write(*,*) '##### In wdiscpbl.f ########'
c      write(*,*) 'gr         :', gr

      if (abs(fltype).eq.20) write(11) (gr(j),j=1,scalar)
c
c     Write velocity and scalar
c
      do i=1,m
         do z=1,nzc
            call getxyp(pln,z,i,ur)
            do y=1,nyp
               do x=1,nx/2
                  urx(x)=pln(x,y)
               end do
               write(11) urx
            end do
         end do
      end do
      close(unit=11)

      end subroutine wdiscpbl
