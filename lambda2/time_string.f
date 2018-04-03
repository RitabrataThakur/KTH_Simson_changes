c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/lambda2/time_string.f $
c $LastChangedDate: 2010-09-03 15:04:26 +0200 (Fri, 03 Sep 2010) $
c $LastChangedBy: mattias@MECH.KTH.SE $
c $LastChangedRevision: 1506 $
c
c ***********************************************************************
      subroutine time_string(cdt)
c
c     Construct string in the format '19-DEC-2005 22:47:06'
c
      implicit none

      integer i

      integer val(8)
      character*20 cdt
      character*3 monc

      call date_and_time(values=val)

      if (val(2).eq.1) then
         monc  = 'JAN'
      else if (val(2).eq.2) then
         monc  = 'FEB'
      else if (val(2).eq.3) then
         monc  = 'MAR'
      else if (val(2).eq.4) then
         monc  = 'APR'
      else if (val(2).eq.5) then
         monc  = 'MAY'
      else if (val(2).eq.6) then
         monc  = 'JUN'
      else if (val(2).eq.7) then
         monc  = 'JUL'
      else if (val(2).eq.8) then
         monc  = 'AUG'
      else if (val(2).eq.9) then
         monc  = 'SEP'
      else if (val(2).eq.10) then
         monc  = 'OCT'
      else if (val(2).eq.11) then
         monc  = 'NOV'
      else if (val(2).eq.12) then
         monc  = 'DEC'
      else
         monc  = 'XXX'
      end if

      write(cdt,'(i2,a1,a3,a1,i4,a1,i2,a1,i2,a1,i2)')
     &     val(3),'-',monc,'-',val(1),' ',val(5),':',val(6),':',val(7)
      do i=1,2
         if (cdt(i:i).eq.' ') then
            cdt(i:i)='0'
         end if
      end do
      do i=13,20
         if (cdt(i:i).eq.' ') then
            cdt(i:i)='0'
         end if
      end do

      end subroutine time_string
