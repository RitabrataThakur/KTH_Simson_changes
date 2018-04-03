c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/pxyst/cont5.f $
c $LastChangedDate: 2011-01-29 15:58:58 +0100 (Sat, 29 Jan 2011) $
c $LastChangedBy: pschlatt@MECH.KTH.SE $
c $LastChangedRevision: 1646 $
c
c ***********************************************************************
      subroutine cont5(nx,ny,u1,xy,iopt,mcont,
     & slev,elev,xaxis,yaxis,heada,headb,mbox,jhard,
     & imode,inunit,utunit,lpgmr)
c
c     Calculates contours of values in an array u1 with associated
c     coordinates xy
c     Plotting is done from mxs to mxe and from mys to mye
c     make mcont contours from slev to elev if mcont positive,
c     if mcont is negative let contxx choose its own levels
c     max abs(mcont)
c     imode :
c     1 interactive
c     2 batch
c     3 input for batch
c
c     cont5 also outputs a pgmr file into plot.pgmr
c
c     originally by Dan Henningson
c     revised by Anders Lundbladh 890901
c     revised by Henningson spring 91
c     revised by Lundbladh 910515
c     cont5 920902
c
      integer imode,inunit,utunit,mcont,mbox,iopt
      integer nx,ny,naddim,nxdim
      parameter(naddim=50000,nxdim=4*naddim)
      real u1(nx,ny),xy(nx,ny,2)
      real acont(naddim),xcont(nxdim),ycont(nxdim)
      real slev,elev
      character*(*) xaxis,yaxis,heada,headb
      integer nlev(naddim),work(naddim),nad(naddim)
      integer jhard
      logical lpgmr

      real xs,xe,ys,ye,gl,gu,xmin,xmax,ymin,ymax
      character text*80
      character*1 ans
      integer ihard,is,idash,isnext,nd,ic,i,mxs,mxe,mys,mye
      integer ncopt,ncont
      common /hard/ihard

      ihard=jhard

      mxs=1
      mxe=nx
      mys=1
      mye=ny

      write(*,*)
      write(*,*) 'pxyst $Rev: 1646 $'
      write(*,*)

      if (iopt.eq.1) then
         write(*,*) 'Do you want to make a subarea selection (y/n)'
         read(inunit,9000) ans
         if (imode.ge.2) write(utunit,9000) ans
 9000    format(a1)
         if (ans.eq.'Y'.or.ans.eq.'y') then
           if (imode.eq.3) then
             write(*,*)' Give start and end for ',
     &             xaxis(1:index(xaxis,'$')-1)
           else
             write(*,*)' Give start and end for ',
     &             xaxis(1:index(xaxis,'$')-1),' (',xy(mxs,mys,1),
     &             ' - ',xy(mxe,mye,1),')'
           end if
           read(inunit,*) xs,xe
           if (imode.ge.2) write(utunit,*) xs,xe
           if (imode.ne.3) then
             do 2000 i=mxs,mxe-1
               if (xs.lt..5*(xy(i,mys,1)+xy(i+1,mys,1))) goto 1000
 2000        continue
             i=mxe
 1000        mxs=i
             do 2010 i=mxs,mxe-1
               if (xe.lt..5*(xy(i,mys,1)+xy(i+1,mys,1))) goto 1010
 2010        continue
             i=mxe
 1010        mxe=i
             write(*,*) 'nearest values', xy(mxs,mys,1),xy(mxe,mys,1)
           end if
           if (imode.eq.3) then
             write(*,*)' Give start and end for ',
     &             yaxis(1:index(yaxis,'$')-1)
           else
             write(*,*)' Give start and end for ',
     &             yaxis(1:index(yaxis,'$')-1),' (',xy(mxs,mys,2),
     &             ' - ',xy(mxe,mye,2),')'
           end if
           read(inunit,*) ys,ye
           if (imode.ge.2) write(utunit,*)ys,ye
           if (imode.ne.3) then
             do 2020 i=mys,mye-1
               if (ys.lt..5*(xy(mxs,i,2)+xy(mxs,i+1,2))) goto 1020
 2020        continue
             i=mye
 1020        mys=i
             do 2030 i=mys,mye-1
               if (ye.lt..5*(xy(mxs,i,2)+xy(mxs,i+1,2))) goto 1030
 2030        continue
             i=mye
 1030        mye=i
             write(*,*) 'nearest values', xy(mxs,mys,2),xy(mxs,mye,2)
           end if
         end if
      end if
c
c     Set levels for pgmr file
c
      if (lpgmr) then
        write(*,*) 'Give upper and lower level for greyscale'
        write(*,*) 'Equal for automatic'
        read(inunit,*) gl,gu
        if (imode.ge.2) write(utunit,*) gl,gu
      end if
      if (imode.eq.3) return

      if (mcont.gt.1) then
        do 100 i=1,mcont
          acont(i)=slev+(elev-slev)*real(i-1)/real(mcont-1)
100     continue
      end if
      if (mcont.eq.1) then
        acont(1)=slev
      end if
      ncopt=2
      if (mcont.lt.0) ncopt=1
      ncont=abs(mcont)

      call contxx(nx,ny,mxs,mxe,mys,mye,xy,u1,ncopt,ncont,acont,
     *         nxdim,xcont,ycont,naddim,nad,nlev,naddim,work)

      write(text,1222) ncont,acont(1),acont(ncont),acont(2)-acont(1)
1222  format(i3,' levels, start:',f11.8,' end:',f10.8,
     &' spacing:',f10.8,'$')
      if (abs(acont(1)).GT..1.or. abs(acont(ncont)).GT..1)
     & write(text,1223) ncont,acont(1),acont(ncont),acont(2)-acont(1)
1223  format(i3,' levels, start:',f11.4,' end:',f10.4,
     &' spacing:',f10.4,'$')
      xmin=xy(mxs,mys,1)
      xmax=xy(mxe,mye,1)
      ymin=xy(mxs,mys,2)
      ymax=xy(mxe,mye,2)
      call xterm1
      call plarea(xmin,xmax,ymin,ymax,mbox)
      call axis(xmin,0.,xmax,ymin,0.,ymax)
      call head1(heada)
      call head2(headb)
      call head3(text)
      call taxis(xaxis,yaxis)

      isnext=1
      nd=nad(1)
      do 500 ic=2,nd+1
          is=isnext
          isnext=nad(ic)
          idash=1
          if (acont(nlev(ic-1)).lt.0) idash=2
          call curve(xcont(is),ycont(is),isnext-is,idash)
 500  continue
      call plend
      if (lpgmr) call pgmr(u1,xy,mxs,mxe,mys,mye,nx,ny,mbox,gl,gu)

      return

      end subroutine cont5
