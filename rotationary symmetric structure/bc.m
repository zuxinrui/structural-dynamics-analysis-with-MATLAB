function [kk,mm,ff]=femApplybc1(kk,mm,ff,bcdof,bcval)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
 ni=length(bcdof);

 
 for ii=1:ni
   if bcdof(ii)==1
     for j=1:ni                      % attention!!!this must change if dofs are changed!!!
       kk(ii,j)=0;
       kk(j,ii)=0;
       mm(ii,j)=0;
       mm(j,ii)=0;
       ff(j)=ff(j)-bcval(ii)*kk(j,ii);
     end
     kk(ii,ii)=1;
     ff(ii)=bcval(ii);
   end
 end

end

