% ***********************************************************************
%
% $HeadURL: https://www.mech.kth.se/svn/simson/trunk/matlab/textline.m $
% $LastChangedDate: 2006-09-22 11:01:19 +0200 (Fri, 22 Sep 2006) $
% $LastChangedBy: mattias@MECH.KTH.SE $
% $LastChangedRevision: 147 $
%
% ***********************************************************************
function [str2] = textline(text)
%        1234567890123456789012345678901234567890
  str = '                                        ';
  str2 = [str str];

  for i = 1:min(80,length(text))
    str2(i)=text(i);
  end