c ***********************************************************************
c
c $HeadURL: https://www.mech.kth.se/svn/simson/trunk/bla/comment.f $
c $LastChangedDate: 2007-12-13 22:41:14 +0100 (Thu, 13 Dec 2007) $
c $LastChangedBy: pschlatt@MECH.KTH.SE $
c $LastChangedRevision: 1069 $
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
         if (line.ne.'#') exit
      end do

      backspace(iunit)

      end subroutine comment
