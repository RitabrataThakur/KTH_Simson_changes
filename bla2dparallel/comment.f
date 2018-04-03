c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/comment.f $
c $LastChangedDate: 2007-10-31 10:21:52 +0100 (Wed, 31 Oct 2007) $
c $LastChangedBy: qiang@MECH.KTH.SE $
c $LastChangedRevision: 787 $
c
c ***********************************************************************
      subroutine comment(iunit)
c
c     Reads from the formatted unit iunit until the next
c     non-comment line (i.e. without starting with an '#')
c
      implicit none

      integer iunit
      character line

      do
         read(iunit,*) line
c         write(*,'(3a)') '>>',line,'<<'
c         stop
         if (line.eq.'#') then
         else
            exit
         end if
      end do
      
      backspace(iunit)
 
      end subroutine comment
