function [index]=femEldof(nd,No_nel,No_dof)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
   k=0;
   for i=1:No_nel
     start = (nd(i)-1)*No_dof;
       for j=1:No_dof
         k=k+1;
         index(k)=start+j;
       end
   end

end

