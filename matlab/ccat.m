% ***********************************************************************
%
% $HeadURL: https://www.mech.kth.se/svn/simson/trunk/matlab/ccat.m $
% $LastChangedDate: 2006-09-22 11:01:19 +0200 (Fri, 22 Sep 2006) $
% $LastChangedBy: mattias@MECH.KTH.SE $
% $LastChangedRevision: 147 $
%
% ***********************************************************************
function res=ccat(dim,varargin);
res=[];
for ind=1:length(varargin);
  res=cat(dim,res,varargin{ind});
end
