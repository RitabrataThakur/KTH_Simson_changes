#!/bin/tcsh
# ***********************************************************************
# 
# $HeadURL: https://www.mech.kth.se/svn/simson/trunk/lambda2/conv_jpg.sh $
# $LastChangedDate: 2007-04-01 16:50:18 +0200 (Sun, 01 Apr 2007) $
# $LastChangedBy: pschlatt@MECH.KTH.SE $
# $LastChangedRevision: 510 $
#
# ***********************************************************************

foreach f (*.tiff)
echo doing $f
convert -quality 100 $f $f:r.jpg
rm -f $f
end
