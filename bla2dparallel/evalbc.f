c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/evalbc.f $
c $LastChangedDate: 2007-10-29 14:23:13 +0100 (Mon, 29 Oct 2007) $
c $LastChangedBy: qiang@MECH.KTH.SE $
c $LastChangedRevision: 786 $
c
c ***********************************************************************
      subroutine evalbc(bcm,bcdr,bcdi,ibc,zb,
     &     dj1omr,dj1omi,dj1omh,
     &     dj1ur,dj1ui,dj1uh1,dj1uh2,dj1uh3,k,bu1jr,bu1ji,alfa,xb)
c     
c     Evaluates the boundary conditions for the homogeneous
c     and particular solutions(U,V,W: basic flow+disturbance within fringe)
c     the boundary conditions are selected as follows :
c    ibc    bc
c     0      u=v=w=0
c     1      du=dv=dw=0
c     2      d2u=d2v=d2w=0
c     3      d3u=d3v=d3w=0
c    10      du+ku=dv+kv=dw+kw=0
c    11      d2u+kdu=d2v+kdv=d2w+kdw=0
c    12      d3u+kd2u=d3v+k2dv=d3w+kd2w=0
c    20      dv+kv=0 d2v+kdv=0 omy=0
c   100      u=U, v=V w=W
c   101      du=dU dv=dV dw=dW
c   110      du+ku=dU+kU dv+kv=dV+kV dw+kw=dW+kW
c   120      dv+kv=dV+kV d2v+kdv=d2V+kdV omy=0
c   130      u=U, du=dU, w=W
c   140      u=U, dv=dV, w=W
c   150      u=U, du-vx=0, dw=0
c
c     (d is the vertical partial derivative ddy)
c
c     Layout of bcm :
c    
c     ( buh1 buh2 buh3 )
c     ( bvh1 bvh2 bvh3 )
c     ( bwh1 bwh2 bwh3 )
c     the first row corresponds to the first of the three bc
c     given for each ibc above etc.c
c     the first column corresponds to the first homogeneous 
c     solution vh1, the second to vh2, and the third to omh
c
c     dj1q(:,i,j+1) is the jth derivative of the ith velocity component of q
c     (u <=> i=1,v <=> i=2,w <=> i=3)
c     where q is one of vr,vi,vh1,vh2,omr,omi or omh
c     the other velocity components are calculated through continuity
c     (for the omh solutions v=0, for vh1, vh2 solutions omegay=0)
c
c     vh1,vh2,omh are real, thus :
c     v and y-derivatives of v are real
c     u,w and y-derivatives are imaginary
c     horizontal (first) derivatives of v are imaginary
c     horizontal (first) derivatives of u,w are real
c
      implicit none
c
      include 'par.f'
c
      integer ibc,zb
      real dj1ur(memnx*mbz,3,5),dj1ui(memnx*mbz,3,5)
      real dj1uh1(memnx*mbz,3,5),dj1uh2(memnx*mbz,3,5)
      real dj1uh3(memnx*mbz,3,4)
      real dj1omr(memnx*mbz,5),dj1omi(memnx*mbz,5)
      real dj1omh(memnx*mbz,5)
      real bcdr(memnx*mbz,3),bcdi(memnx*mbz,3),bcm(memnx*mbz,3,3)
      real bu1jr(nxp/2+1,3+scalar,3),bu1ji(nxp/2+1,3+scalar,3)
      real k(nx/2*mbz),alfa(nx/2*mbz)
      
      integer xz,nxz,i,x,j,n1,xb
      real tmp
      
      nxz=memnx*mbz

      if (ibc.ge.0.and.ibc.le.3) then
         do i=1,3
            do xz=1,nxz
               bcdr(xz,i)=dj1ur(xz,i,ibc+1)
               bcdi(xz,i)=dj1ui(xz,i,ibc+1)
               bcm(xz,i,1)=dj1uh1(xz,i,ibc+1)
               bcm(xz,i,2)=dj1uh2(xz,i,ibc+1)
               bcm(xz,i,3)=dj1uh3(xz,i,ibc+1)
            end do
         end do
      end if
      if (ibc.ge.10.and.ibc.le.12) then
         do i=1,3
            do xz=1,nxz
               x=xz+xb
               bcdr(xz,i)=dj1ur(xz,i,ibc-8)+k(x)*dj1ur(xz,i,ibc-9)
               bcdi(xz,i)=dj1ui(xz,i,ibc-8)+k(x)*dj1ui(xz,i,ibc-9)
               bcm(xz,i,1)=dj1uh1(xz,i,ibc-8)+k(x)*dj1uh1(xz,i,ibc-9)
               bcm(xz,i,2)=dj1uh2(xz,i,ibc-8)+k(x)*dj1uh2(xz,i,ibc-9)
               bcm(xz,i,3)=dj1uh3(xz,i,ibc-8)+k(x)*dj1uh3(xz,i,ibc-9)
            end do
         end do
      end if
      if (ibc.eq.20.or.ibc.eq.120) then
         do xz=1,nxz
            x=xz+xb
            bcdr(xz,1)=dj1ur(xz,2,2)+k(x)*dj1ur(xz,2,1)
            bcdi(xz,1)=dj1ui(xz,2,2)+k(x)*dj1ui(xz,2,1)
            bcdr(xz,2)=dj1ur(xz,2,3)+k(x)*dj1ur(xz,2,2)
            bcdi(xz,2)=dj1ui(xz,2,3)+k(x)*dj1ui(xz,2,2)
            bcdr(xz,3)=dj1omr(xz,1)
            bcdi(xz,3)=dj1omi(xz,1)
            bcm(xz,1,1)=dj1uh1(xz,2,2)+k(x)*dj1uh1(xz,2,1)
            bcm(xz,1,2)=dj1uh2(xz,2,2)+k(x)*dj1uh2(xz,2,1)
            bcm(xz,1,3)=0.0
            bcm(xz,2,1)=dj1uh1(xz,2,3)+k(x)*dj1uh1(xz,2,2)
            bcm(xz,2,2)=dj1uh2(xz,2,3)+k(x)*dj1uh2(xz,2,2)
            bcm(xz,2,3)=0.0
            bcm(xz,3,1)=0.0
            bcm(xz,3,2)=0.0
            bcm(xz,3,3)=dj1omh(xz,1)
         end do
      end if
      if (ibc.ge.100.and.ibc.le.101) then
         do i=1,3
            do xz=1,nxz
               bcdr(xz,i)=dj1ur(xz,i,ibc-99)
               bcdi(xz,i)=dj1ui(xz,i,ibc-99)
               bcm(xz,i,1)=dj1uh1(xz,i,ibc-99)
               bcm(xz,i,2)=dj1uh2(xz,i,ibc-99)
               bcm(xz,i,3)=dj1uh3(xz,i,ibc-99)
            end do
         end do
         if (zb.eq.1) then
            do i=1,3
               do xz=1,memnx
                  x=xz+xb
                  bcdr(xz,i)=bcdr(xz,i)-bu1jr(x,i,ibc-99)
                  bcdi(xz,i)=bcdi(xz,i)-bu1ji(x,i,ibc-99)
               end do
            end do
         end if
      end if
      if (ibc.eq.110) then 
         do i=1,3
            do xz=1,nxz
               x=xz+xb
               bcdr(xz,i)=dj1ur(xz,i,2)+k(x)*dj1ur(xz,i,1)
               bcdi(xz,i)=dj1ui(xz,i,2)+k(x)*dj1ui(xz,i,1)
               bcm(xz,i,1)=dj1uh1(xz,i,2)+k(x)*dj1uh1(xz,i,1)
               bcm(xz,i,2)=dj1uh2(xz,i,2)+k(x)*dj1uh2(xz,i,1)
               bcm(xz,i,3)=dj1uh3(xz,i,2)+k(x)*dj1uh3(xz,i,1)
            end do
         end do
         if (zb.eq.1) then
            do i=1,3
               do xz=1,memnx
                  x=xz+xb
                  bcdr(xz,i)=bcdr(xz,i)-bu1jr(x,i,2)-k(x)*bu1jr(x,i,1)
                  bcdi(xz,i)=bcdi(xz,i)-bu1ji(x,i,2)-k(x)*bu1ji(x,i,1)
               end do
            end do
         end if
      end if 
      if (ibc.eq.120.and.zb.eq.1) then 
         do xz=1,memnx
            x=xz+xb
            bcdr(xz,1)=bcdr(xz,1)-bu1jr(x,2,2)-k(x)*bu1jr(x,2,1)
            bcdi(xz,1)=bcdi(xz,1)-bu1ji(x,2,2)-k(x)*bu1ji(x,2,1)
            bcdr(xz,2)=bcdr(xz,2)-bu1jr(x,2,3)-k(x)*bu1jr(x,2,2)
            bcdi(xz,2)=bcdi(xz,2)-bu1ji(x,2,3)-k(x)*bu1ji(x,2,2)
         end do
      end if
      if (ibc.eq.130) then
         do i=1,3,2
            do xz=1,nxz
               bcdr(xz,i)=dj1ur(xz,i,1)
               bcdi(xz,i)=dj1ui(xz,i,1)
               bcm(xz,i,1)=dj1uh1(xz,i,1)
               bcm(xz,i,2)=dj1uh2(xz,i,1)
               bcm(xz,i,3)=dj1uh3(xz,i,1)
            end do
         end do
         do xz=1,nxz
            bcdr(xz,2)=dj1ur(xz,1,2)
            bcdi(xz,2)=dj1ui(xz,1,2)
            bcm(xz,2,1)=dj1uh1(xz,1,2)
            bcm(xz,2,2)=dj1uh2(xz,1,2)
            bcm(xz,2,3)=dj1uh3(xz,1,2)
         end do
         if (zb.eq.1) then
            do xz=1,memnx
               x=xz+xb
               bcdr(xz,1)=bcdr(xz,1)-bu1jr(x,1,1)
               bcdi(xz,1)=bcdi(xz,1)-bu1ji(x,1,1)
               bcdr(xz,2)=bcdr(xz,2)-bu1jr(x,1,2)
               bcdi(xz,2)=bcdi(xz,2)-bu1ji(x,1,2)
               bcdr(xz,3)=bcdr(xz,3)-bu1jr(x,3,1)
               bcdi(xz,3)=bcdi(xz,3)-bu1ji(x,3,1)
            end do
         end if
      end if
      if (ibc.eq.140) then
         do i=1,3,2
            do xz=1,nxz
               bcdr(xz,i)=dj1ur(xz,i,1)
               bcdi(xz,i)=dj1ui(xz,i,1)
               bcm(xz,i,1)=dj1uh1(xz,i,1)
               bcm(xz,i,2)=dj1uh2(xz,i,1)
               bcm(xz,i,3)=dj1uh3(xz,i,1)
            end do
         end do
         do xz=1,nxz
            bcdr(xz,2)=dj1ur(xz,2,2)
            bcdi(xz,2)=dj1ui(xz,2,2)
            bcm(xz,2,1)=dj1uh1(xz,2,2)
            bcm(xz,2,2)=dj1uh2(xz,2,2)
            bcm(xz,2,3)=dj1uh3(xz,2,2)
         end do
         if (zb.eq.1) then
            do xz=1,memnx
               x=xz+xb
               bcdr(xz,1)=bcdr(xz,1)-bu1jr(x,1,1)
               bcdi(xz,1)=bcdi(xz,1)-bu1ji(x,1,1)
               bcdr(xz,2)=bcdr(xz,2)-bu1jr(x,2,2)
               bcdi(xz,2)=bcdi(xz,2)-bu1ji(x,2,2)
               bcdr(xz,3)=bcdr(xz,3)-bu1jr(x,3,1)
               bcdi(xz,3)=bcdi(xz,3)-bu1ji(x,3,1)
            end do
         end if
      end if
      if (ibc.eq.150) then
         do xz=1,nxz
            x=x+xb
c
c     bc 1 rhs : up
c
            bcdr(xz,1)=dj1ur(xz,1,1)
            bcdi(xz,1)=dj1ui(xz,1,1)
c
c     bc 2 rhs : dup-vxp
c
            bcdr(xz,2)=dj1ur(xz,1,2)+alfa(x)*dj1ui(xz,2,1)
            bcdi(xz,2)=dj1ui(xz,1,2)-alfa(x)*dj1ur(xz,2,1)
c
c     bc 3 rhs : wyp
c
            bcdr(xz,3)=dj1ur(xz,3,2)
            bcdi(xz,3)=dj1ur(xz,3,2)
c
c     bc 1 matrix row 1 : u
c
            bcm(xz,1,1)=dj1uh1(xz,1,1)
            bcm(xz,1,2)=dj1uh2(xz,1,1)
            bcm(xz,1,3)=dj1uh3(xz,1,1)
c
c     bc 2 matrix row 2 : du-vx
c
            bcm(xz,2,1)=dj1uh1(xz,1,2)-alfa(x)*dj1uh1(xz,2,1)
            bcm(xz,2,2)=dj1uh2(xz,1,2)-alfa(x)*dj1uh2(xz,2,1)
            bcm(xz,2,3)=dj1uh3(xz,1,2)-alfa(x)*dj1uh3(xz,2,1)
c
c     bc 3 matrix row 3 : dw
c
            bcm(xz,3,1)=dj1uh1(xz,3,2)
            bcm(xz,3,2)=dj1uh2(xz,3,2)
            bcm(xz,3,3)=dj1uh3(xz,3,2)
         end do
c
c     For beta=0 subract the base flow part
c
         if (zb.eq.1) then
            do xz=1,memnx
               x=xz+xb
c     bc 1 rhs : U
               bcdr(xz,1)=bcdr(xz,1)-bu1jr(x,1,1)
               bcdi(xz,1)=bcdi(xz,1)-bu1ji(x,1,1)
            end do
         end if
      end if
c
c     Note that the first and third row in bcm are imaginary
c     for all boundary conditions except 20 and 120
c     but stored as real
c     (see comments at the beginning of this file)
c     therefore divide the first and third component of bcd
c     by i
c
      if (ibc.ne.20.and.ibc.ne.120.and.ibc.ne.130.and.ibc.ne.150) then
         n1=1
         if (zb.eq.1.and.xb.eq.0) n1=2
         do j=1,3,2
            do xz=n1,nxz
               tmp=bcdr(xz,j)
               bcdr(xz,j)=bcdi(xz,j)
               bcdi(xz,j)=-tmp
            end do
         end do
      end if
c
c     For ibc=130,150 all matrix elements are imaginary
c
      if (ibc.eq.130.or.ibc.eq.150) then
         n1=1
         if (zb.eq.1) n1=2
         do j=1,3
            do xz=n1,nxz
               tmp=bcdr(xz,j)
               bcdr(xz,j)=bcdi(xz,j)
               bcdi(xz,j)=-tmp
            end do
         end do
      end if
      
      end subroutine evalbc
