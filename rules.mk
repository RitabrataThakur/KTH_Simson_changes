# ***********************************************************************
#
# $HeadURL: https://www.mech.kth.se/svn/simson/trunk/rules.mk $
# $LastChangedDate: 2007-11-12 11:26:42 +0100 (Mon, 12 Nov 2007) $
# $LastChangedBy: mattias@MECH.KTH.SE $
# $LastChangedRevision: 846 $
#
# ***********************************************************************

# Mark standard targets as phony (not real file targets).
# This helps to avoid problems if there is also a file with the same name
.PHONY: all clean dist doc install release test

# Disable implicit rules (some of them may otherwise cause problems)
.SUFFIXES:

# Rule for building fixed form fortran
%.o: %.f
	@echo " Compiling $*.f with $(F90) $(FLAGS) $(ENDIAN) -c" ;\
	$(F90) $(FLAGS) $(ENDIAN) -c $*.f ;\

# Rules for recursion into sub directories for standard targets
%.all:
	@cd $* && make -f Makefile.config all

%.clean:
	@cd $* && make -f Makefile.config clean

%.dist:
	@cd $* && make -f Makefile.config dist

%.install:
	@cd $* && make -f Makefile.config install

%.test:
	@cd $* && make -f Makefile.config test
