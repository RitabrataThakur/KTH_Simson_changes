#
# Compilers and flags
#
FC                 = 
F90                = ifort
FLAGS              = -r8 -i4 -fpp -O3 -fp-model precise -DLINUX_INTEL
LDFLAGS            = -r8 -i4 -fpp -O3 -fp-model precise -DLINUX_INTEL
#
# MPI compiler and flags
#
MPI_F90            = 
MPI_FLAGS          = 
MPI_LDFLAGS        = 
#
# OpenMP compiler and flags
#
OMP_F90            = ifort
OMP_FLAGS          = 
OMP_LDFLAGS        = 
#
# Additional libraries
#
LIBS               = 
#
# Big and little endian options
#
LEND_FLAG          = -convert little_endian
BEND_FLAG          = -convert big_endian
#
# Flag showing whether OpenMP should be used or not
# Options: [yes,no] (empty is interpreted as no)
# Note that this option can be overriden by Makefile option as
# make [mpi=yes][omp=yes]
#
omp                =
#
# Flag showing whether MPI should be used or not
# Options: [yes,no] (empty is interpreted as no)
# Note that this option can be overriden by Makefile option as
# make [mpi=yes][omp=yes]
#
mpi                =
#
# Flag showing if little endian should be used or not
# make [lendian=yes][lendian=no] empty means that big endian will be used
#
lendian            = yes
#
# Global Makefile settings.
#
# prefix is the base directory where Simson will be installed.
# (Binaries will go into /home/jose/SIMSONsvn/x86_64/bin/)
#
prefix             = /home/jose/SIMSONsvn/x86_64
#
# EXESUFFIX is an extension added to each executable e.g. _dbl for
# double precision.
#
EXESUFFIX          = 
#
# MACHINE identifies the host machine type
# Used as default directory for installation
#
MACHINE            = x86_64
#
# Pdflatex and bibtex (used to generate documentation)
#
PDFLATEX           = /usr/bin/pdflatex
BIBTEX             = /usr/bin/bibtex
#
# Install program to be used
#
INSTALL            = install -C
#
# Compile time parameters
#
PARAM              = par.f
#
# FFT package to be used
#
FFT                = ../fft/
FFTPACK            = cvecfft_acc
#
# SHELL environment variable
#
SHELL              = /bin/bash
