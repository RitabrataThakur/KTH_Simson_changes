# ***********************************************************************
#
# $HeadURL: https://www.mech.kth.se/svn/simson/trunk/xys_add/Makefile $
# $LastChangedDate: 2012-10-10 18:31:16 +0200 (Wed, 10 Oct 2012) $
# $LastChangedBy: taras@MECH.KTH.SE $
# $LastChangedRevision: 1807 $
#
# ***********************************************************************

# compiler flags
F90   = ifort
FLAGS = -r8 -i4 -O3 -warn truncated_source -warn unused

# main program, FFTs and plot routines
main  = xys_add
image = $(main)

# classified subroutines
sub = read_xys.o wxys.o
dep = $(main).o $(sub) 
obj = $(dep) 

# compilation rules
%.o: %.f
	 $(F90) $(FLAGS) -c $*.f 

$(image): $(obj) Makefile
	@echo " " ; echo " Link object files into executable image '$(image)'."
	@echo " "
	$(F90) $(FLAGS) -o $@ $(obj)
	@echo " " ; echo " Done!" ; echo " "

$(dep):	par.f

clean:
	@echo " " ; echo " Clean up all object files and '$(image)'."
	@-rm -f *.o *.F *~ $(image) ; echo " "
