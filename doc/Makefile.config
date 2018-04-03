# ***********************************************************************
#
# $HeadURL: https://www.mech.kth.se/svn/simson/trunk/doc/Makefile $
# $LastChangedDate: 2010-10-07 23:31:36 +0200 (Thu, 07 Oct 2010) $
# $LastChangedBy: mattias@MECH.KTH.SE $
# $LastChangedRevision: 1562 $
#
# ***********************************************************************

# Makefile for documentation

# Gobal configuration file
include ../config.mk

# Gobal Makefile rules
include ../rules.mk

TEXFILES=\
simson-user-guide.tex \
simson-svn-guide.tex \
simson-paraview-guide.tex

PDFFILES=\
simson-user-guide.pdf \
simson-svn-guide.pdf \
simson-paraview-guide.pdf

# Printed version

#TEXFILESv4=\
#simson-user-guide-v4.0.tex

#PDFFILESv4=\
#simson-user-guide-v4.0.pdf



###########################################################################

#all: $(PDFFILES) $(PDFFILESv4)

all: $(PDFFILES)

.SUFFIXES: .pdf .tex

%.pdf: %.tex simson.cls simson.bst fluids.bib Makefile
	pdflatex $<
	bibtex $*
	pdflatex $<
	pdflatex $<
	pdflatex $<

clean:
	@echo " " ; echo " Cleaning up all auxiliary files"
	@-rm -f *.dvi *.log *.aux *.pdf *.out *.toc *.lof *.bbl *.blg ; echo " "
