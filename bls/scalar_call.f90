program scalar_call
implicit none

integer ith
real rt_tempflag
rt_tempflag = 1
real scalar
scalar = 1
real rt_scalar_profile
integer y
integer nyp
nyp = 11
real eta_scalar_file
real scalar_from_file
integer ivarvisc 
ivarvisc = 1
integer rt_scalar_flag
rt_scalar_flag = 1
real rt_pr
real rt_ri
character*40 rt_scalar_profile


 open(unit=10,status='old',file='fake_bls_for_running.i')


 read(10,*) ivarvisc
      do ith = 1,scalar 


         read(10,*) rt_scalar_flag(ith)

write(*,*) '  scalar flag                           :',rt_scalar_flag(ith)

         if (rt_scalar_flag(ith).ne.0) then
            read(10,*) rt_scalar_profile(ith)
         end if

         read(10,*) rt_pr(ith)
write(*,*) '  Prandtl No                          :',rt_pr(ith)
         read(10,*) rt_ri(ith)
write(*,*) '  Prandtl No                          :',rt_ri(ith)
      end do



do ith = 1,scalar
         if (rt_tempflag.ne.0)then
            open(54, status='old', file=rt_scalar_profile(ith))
            do y = 1,nyp
               read(54,*)eta_scalar_file(y,ith),scalar_from_file(y,ith)
            write(*,*) 'eta_scalar: ',eta_scalar_file(y,ith),'scalar_profile ',scalar_from_file(y,ith)
            end do
            close(54)
         end if
      end do


end program scalar_calling_testing
