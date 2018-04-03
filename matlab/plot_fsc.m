% ***********************************************************************
%
% $HeadURL: https://www.mech.kth.se/svn/simson/trunk/matlab/plot_fsc.m $
% $LastChangedDate: 2006-12-04 20:05:01 +0100 (Mon, 04 Dec 2006) $
% $LastChangedBy: pschlatt@MECH.KTH.SE $
% $LastChangedRevision: 392 $
%
% ***********************************************************************
clear all
close all

a=load('fsc.txt');

figure
hold on
plot(a(:,1),a(:,1)-a(:,2),'b-' );
plot(a(:,1),a(:,3)       ,'b--');
plot(a(:,1),a(:,4)       ,'b:' );
plot(a(:,1),a(:,5)       ,'b-.');

plot(a(:,1),a(:,6),'r-' );
plot(a(:,1),a(:,7),'r--');
plot(a(:,1),a(:,8),'r:' );

if (size(a,2)==11) 
  plot(a(:,1),1-a(:,9) ,'g-' );
  plot(a(:,1), -a(:,10),'g--');
  plot(a(:,1), -a(:,11),'g:' );
  legend('\eta-f','f\prime','f\prime\prime','f\prime\prime\prime','g','g\prime','g\prime\prime','1-\theta','-\theta\prime','-\theta\prime\prime')
else
  legend('\eta-f','f\prime','f\prime\prime','f\prime\prime\prime','g','g\prime','g\prime\prime')
end 
xlabel('\eta')

