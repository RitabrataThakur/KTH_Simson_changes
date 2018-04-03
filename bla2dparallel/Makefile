# ***********************************************************************
# 
# $HeadURL: https://www.mech.kth.se/svn/simson/branches/bla2dparallel/Makefile $
# $LastChangedDate: 2011-01-15 17:41:30 +0100 (Sat, 15 Jan 2011) $
# $LastChangedBy: ilak@MECH.KTH.SE $
# $LastChangedRevision: 1634 $
#
# ***********************************************************************
#  **********************************************************************
#  ** Makefile for BLA                                                 **
#  ** Philipp Schlatter, November 2005                                 **
#  **********************************************************************
#
#  usage: make [omp=yes][mpi=yes][arch=Linux|Linux_PGI|AIX|BG][wrapper=yes]
#
#  tested compilers:
#  - Intel ifort v8, v9
#  - PGI Fortran
#  - IBM XLF


#  **********************************************************************
#  **  OBJECT  FILES                                                   **
#  **********************************************************************

# main program
ifeq "$(wrapper)" "yes"
	PROG = wrapper
else
	PROG = bla
endif

# files dependent on par.f
ppar     = ppar.o
start    = preprbl.o prhs.o getdt.o wall_time.o rkcoeff.o stopnow.o rescale.o
input    = rparambl.o comment.o
bld      = bflow.o step.o dstep.o xlim.o blfou.o blshift.o cbflow.o \
           splfre.o rbflow.o bfvort.o rbla.o bfgrad.o
read     = rdiscbl.o rwavebl.o lenbox.o rstreakbl.o rwavesbl.o
write    = wdiscbl.o 
fftsub   = fft2db.o fft2df.o
fft      = cvecfft_acc.o
stat     = boxamp.o boxcfl.o wampbl.o wcflbl.o wplbl.o wplbl_serial.o \
           boxext.o plnmin.o plnmax.o wextbl.o boxcflbf.o
xys      = wxys.o boxxys.o
get      = getxy.o getxz.o 
put      = putxy.o putxz.o 
bc       = evald.o evalbc.o solvbc.o cwallbc.o
high     = linearbl.o nonlinbl.o 
forc     = fring.o locf.o gtrip.o trip.o tzran.o ffun.o
galil    = xysh.o xzsh.o fshift.o
les      = init_filter.o bandsolver.o hpf_y.o lpf_y.o rparamles.o les.o \
           sgs_stress.o filter_cutoff.o mhd_poi.o
pressure = fringp.o linearp.o nonlinp.o wdiscp.o
mpi_dep  = getpxz_z.o getpxz_x.o putpxz_z.o putpxz_x.o getpxz.o putpxz.o
wallr    = #rparamwallrough.o auxwallrough.o steadywallrough.o \
           #  updatewallrough.o
jet      = jetinit.o
3dbase   = bfvort3.o pre_nlin_3d.o
wrapperdep  = bla.o

# files not dependent on par.f
diffint  = dcheb.o icheb.o d2cheb.o rdcheb.o
solv     = setmat.o trid.o setmacrhs.o setmac.o tridc.o crhsc.o 
bl       = cubip.o 
mdep     = ctim.o cflush.o time_string.o
corr     = corrd.o corrp.o corrd2.o
misc     = ran2.o cnxtmp.o
spl      = spline.o spliet.o
mpi_idep = mpi_findout.o

# dependencies
dep_  = $(PROG).o $(start) $(read) $(write) $(stat) $(get) $(put) \
        $(high) $(fftsub) $(bld) $(bc) $(forc) $(galil) $(xys) $(input) \
        $(les) $(fft) $(pressure) $(wallr) $(jet) $(3dbase)
idep_ = $(diffint) $(solv) $(bl) $(mdep) $(corr) $(misc) $(spl)

ifeq "$(wrapper)" "yes"
	dep_ += $(wrapperdep)
endif

# compilation mode
ifeq "$(omp)" "yes"
	ifeq "$(mpi)" "yes"
		dep  = $(dep_) $(mpi_dep)
		idep = $(idep_) $(mpi_idep)
		MODE = (MPI & OpenMP)
	else	
		dep  = $(dep_)
		idep = $(idep_) 
		MODE = (OpenMP)
	endif
else
	ifeq "$(mpi)" "yes"
		dep  = $(dep_) $(mpi_dep)
		idep = $(idep_) $(mpi_idep)
		MODE = (MPI)
	else	
		dep  = $(dep_)
		idep = $(idep_) 
		MODE = (serial)
	endif
endif

obj   = $(dep) $(idep) $(ppar)



#  **********************************************************************
#  **  TARGET  MACHINES                                                **
#  **********************************************************************

# determining target architecture
arch=$(shell uname)

# Special machines: 
# Mekanik with 32bit AMD processors (e.g. gampen)
ifeq "$(shell uname -r)" "2.6.32-5-amd64"
  arch=Linux_PGI
endif
# Mekanik with 64bit AMD processors (annapurna)
ifeq "$(shell uname -r)" "2.6.18-6-amd64"
  arch=Linux_PGI
endif
# Mekanik with 64bit AMD processors (damavand)
ifeq "$(shell uname -r)" "2.6.26-2-amd64"
  arch=Linux_PGI
endif
# ekman
ifeq "$(shell uname -r)" "2.6.18-128.7.1.el5.pdc5noib"
  arch=Linux_PGI
endif
# IBM BlueGene
ifeq "$(shell uname -r)" "2.6.5-7.287.3-pseries64"
  arch=BG
endif
# Lindgren
ifeq "$(shell uname -r)" "2.6.27.48-0.1.1_1.0301.5616-cray_gem_s"
  arch=Cray_PGI
endif



#  **********************************************************************
#  **  COMPILER  OPTIONS                                               **
#  **********************************************************************

ifeq "$(arch)" "Linux"
# cvecfft has to be compiled with -openmp if multiple threads are used.
	CFLAGS = -r8 -i4 -fpp -O3 -warn truncated_source -warn unused -DLINUX_INTEL #-convert big_endian
	ifeq "$(mpi)" "yes"
		FCC = mpif90
		ifeq "$(shell uname -r)" "2.6.9+sc+papi-lenngren1"
			FCC = mpif90 -i_dynamic -mcmodel=large
		endif
		ifeq "$(shell uname -r)" "2.6.18-53.1.14.el5"
			FCC = ifort -Nmpi
		endif
		ifeq "$(shell uname -r)" "2.6.32-7-amd64"
			FCC = ifort -Nmpi
		endif
		ifeq "$(omp)" "yes"
			FC = $(FCC) $(CFLAGS) -openmp -DOPENMP -DMPI
		else	
			FC = $(FCC) $(CFLAGS) -i_dynamic -mcmodel=large -DMPI
		endif
	else	
		FCC = ifort
		ifeq "$(shell uname -r)" "2.6.9+sc+papi-lenngren1"
			FCC = ifort -i_dynamic -mcmodel=large
		endif
		ifeq "$(omp)" "yes"
			FC = $(FCC) $(CFLAGS) -openmp -DOPENMP
		else
			FC = $(FCC) $(CFLAGS) 
		endif
	endif
	LIBS =
endif

#  **********************************************************************

ifeq "$(arch)" "Linux_PGI"	
	CFLAGS = -r8 -i4 -Mpreprocess -O3 -pc 64 -Kieee -fast -fastsse -Mcache_align -DLINUX_PGI -mcmodel=medium # -Mipa=fast,inline:10,reshape -Munroll -Minline=levels:10,reshape -Mvect -Kieee
        #CFLAGS = -r8 -i4 -O2 -byteswapio -cpp 
	ifeq "$(mpi)" "yes"
		FCC = mpif90 #-ccl pgf90
		ifeq "$(omp)" "yes"
			FC = $(FCC) $(CFLAGS) -openmp -DOPENMP -DMPI
		else	
			FC = $(FCC) $(CFLAGS) -DMPI
		endif
	else	
		FCC = pgf90
		ifeq "$(omp)" "yes"
			FC = $(FCC) $(CFLAGS) -mp -DOPENMP
		else
			FC = $(FCC) $(CFLAGS) 
		endif
	endif
	LIBS = 
endif

#  **********************************************************************

ifeq "$(arch)" "Cray_PGI"
        CFLAGS = -mcmodel=medium -r8 -i4 -Mpreprocess -O3 -fast -fastsse -Mcache_align -Mvect -DLINUX_PGI -byteswapio -pc 64 # -Mipa=fast,inline:10,reshape -Munroll -Minline=levels:10,reshape -Mvect -Kieee 
        ifeq "$(mpi)" "yes"
                FCC = ftn
                ifeq "$(omp)" "yes"
                        FC = $(FCC) $(CFLAGS) -openmp -DOPENMP -DMPI
                else
                        FC = $(FCC) $(CFLAGS) -DMPI
                endif
        else
                FCC = pgf90
                ifeq "$(omp)" "yes"
                        FC = $(FCC) $(CFLAGS) -mp -DOPENMP
                else
                        FC = $(FCC) $(CFLAGS)
                endif
        endif
        LIBS =
endif

#  **********************************************************************

ifeq "$(arch)" "AIX"
# typical command: poe bla -procs 4 -hostfile hostfile
	CFLAGS = -q64 -qrealsize=8 -qnosave -O3 -qstrict -WF,-DIBM
	ifeq "$(mpi)" "yes"
		ifeq "$(omp)" "yes"
			FCC = mpxlf_r
			FC = $(FCC) $(CFLAGS) -qsmp=omp -WF,-DOPENMP -WF,-DMPI
		else	
			FCC = mpxlf_r
			FC = $(FCC) $(CFLAGS) -WF,-DMPI
		endif
	else	
		ifeq "$(omp)" "yes"
			FCC = xlf_r
			FC = $(FCC) $(CFLAGS) -qsmp=omp -WF,-DOPENMP
		else
			FCC = xlf
			FC = $(FCC) $(CFLAGS) 
		endif
	endif
	LIBS =
endif

#  **********************************************************************

ifeq "$(arch)" "BG"
	CFLAGS = -qfixed -qrealsize=8 -qzerosize -qtune=440 -qarch=440d -qnosave -qnostrict -qfloat=norngchk -qreport -qhot -O3 -WF,-DIBM -qsuffix=cpp=f
	ifeq "$(mpi)" "yes"
		ifeq "$(omp)" "yes"
			FCC = XXX_NEEDS_MPI_XXX
			FC = $(FCC) $(CFLAGS) -qsmp=omp -WF,-DOPENMP -WF,-DMPI
		else	
			FCC = mpixlf90
			FC = $(FCC) $(CFLAGS) -WF,-DMPI
		endif
	else	
		ifeq "$(omp)" "yes"
			FCC = XXX_NEEDS_MPI_XXX
			FC = $(FCC) $(CFLAGS) -qsmp=omp -WF,-DOPENMP
		else
			FCC = XXX_NEEDS_MPI_XXX
			FC = $(FCC) $(CFLAGS) 
		endif
	endif
	LIBS =
endif

ifeq "$(wrapper)" "yes"
	FC += -DWRAPPER
endif

#  **********************************************************************
#  **  COMPILATION  RULES                                              **
#  **********************************************************************

prog: FIRST_ACTION $(PROG)

$(PROG): $(obj) Makefile 
	@echo '***************************************'
	@echo "** linking $(MODE) $@ on $(arch) with $(FCC)"
	@echo '***************************************'
	@rm -f $@
	@$(FC) -o $@ $(obj) $(LIBS) $(ARPACK)

%.o : %.f
	@echo '***************************************'
	@echo "** compiling $(MODE) $*.f on $(arch) with $(FCC)"
	@echo '***************************************'
	@$(FC) -c $*.f

$(dep):	par.f
$(ppar): par.f

latest_revision: $(dep) $(idep)
	@echo "      write(*,*) 'Latest Revision in Repository: " `grep LastChangedRevision *.f | cut -d " " -f 3 | sort -n | tail -n 1` "'" > 'latest_revision'

latest_revision: ppar.f

ppar.o : latest_revision

.PHONY : clean new FIRST_ACTION 

clean:
	@echo '***************************************'
	@echo "** clean up directory                **"
	@echo '***************************************'
	@rm -f *.o *.F bla wrapper
	@rm -f *~ ifc* latest_revision

new: clean $(PROG)

.SUFFIXES:

FIRST_ACTION:
	@tcsh -f -c 'echo \\033[1\;34m\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\\033[m\\017'
	@tcsh -f -c 'echo \\033[1\;34m\*\* Updating $(PROG) \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \*\*\\033[m\\017'
	@tcsh -f -c 'echo \\033[1\;34m\*\* usage: make \[mpi=yes\]\[omp=yes\]\[arch=Linux\|Linux_PGI\|AIX\|BG\]]\[wrapper=yes\]\ \*\*\\033[m\\017'
	@tcsh -f -c 'echo \\033[1\;34m\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\\033[m\\017'
#	@echo '*************************************************'
#	@echo "** Updating BLA                                **"
#	@echo "** usage: make [mpi=yes][omp=yes][wrapper=yes] **"
#	@echo '*************************************************'
