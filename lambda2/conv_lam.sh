#!/bin/tcsh
# ***********************************************************************
# 
# $HeadURL: https://www.mech.kth.se/svn/simson/trunk/lambda2/conv_lam.sh $
# $LastChangedDate: 2007-04-01 16:50:18 +0200 (Sun, 01 Apr 2007) $
# $LastChangedBy: pschlatt@MECH.KTH.SE $
# $LastChangedRevision: 510 $
#
# ***********************************************************************

foreach f (field.02??.u)
lambda2 $f $f:r
end
