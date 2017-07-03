function [kk]=femAssemble1(kk,k,index)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
eldof = length(index);

for i=1:eldof
   ii=index(i);
   for j=1:eldof
      jj=index(j);
      kk(ii,jj)=kk(ii,jj)+k(i,j);
   end
end

end

