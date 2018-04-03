# ***********************************************************************
#
# $HeadURL: https://www.mech.kth.se/svn/simson/trunk/test/Makefile $
# $LastChangedDate: 2010-09-03 14:59:12 +0200 (Fri, 03 Sep 2010) $
# $LastChangedBy: mattias@MECH.KTH.SE $
# $LastChangedRevision: 1505 $
#
# ***********************************************************************

include ../config.mk
include ../rules.mk

SUBDIRS= \
fltype_-2 \
fltype_3 \

###########################################################################

all:
	@echo " " ; echo " Nothing to be done for all" ; echo " "

test: $(SUBDIRS:%=%.test)

clean: $(SUBDIRS:%=%.clean)
