

E=2.1e11;
A=1e-4;
density=7.3e3; 
node_number=26;
elment_number=61;
nc=zeros(26,2);
nc(1,1)=0;nc(1,2)=0;
nc(2,1)=0;nc(2,2)=1;
for ig=1:12
    for jg=1:2
        nc(ig*2+jg,1)=nc(jg,1)+ig;          
        nc(ig*2+jg,2)=nc(jg,2);
    end
end
en=zeros(61,2);
en(1,1)=1;    en(1,2)=2;                   
en(2,1)=1;    en(2,2)=4;
en(3,1)=2;    en(3,2)=3;
en(4,1)=1;    en(4,2)=3;
en(5,1)=2;    en(5,2)=4;
for i=1:11
   for j=1:5
      en(i*5+j,1)=en(j,1)+2*i;
      en(i*5+j,2)=en(j,2)+2*i;
   end
end
en(61,1)=25;    en(61,2)=26;
ed(1:node_number,1:2)=1;                  
constraint=[1,1;1,2;25,2];
for loopi=1:length(constraint);
    ed(constraint(loopi,1),constraint(loopi,2))=0;
end
dof=0;
for loopi=1:node_number
     for loopj=1:2
        if ed(loopi,loopj)~=0
            dof=dof+1;
            ed(loopi,loopj)=dof;
        end
    end
end
ek=E*A*[1 0 -1 0;0 0 0 0;-1 0 1 0;0 0 0 0];  
em=(density*A)/6*[2 0 1 0; ...            
                   0 2 0 1; ...
                   1 0 2 0; ...
                  0 1 0 2];
k(1:dof,1:dof)=0;                        
m=k;                                  
theta(1:61)=0;
el(1:61)=0;
e2s(1:4)=0;                    
for loopi=1:elment_number
    for zi=1:2
        e2s((zi-1)*2+1)=ed(en(loopi,zi),1);
        e2s((zi-1)*2+2)=ed(en(loopi,zi),2);
    end
    el(loopi)=sqrt((nc(en(loopi,1),1)-nc(en(loopi,2),1))^2+(nc(en(loopi,1),2)-nc(en(loopi,2),2))^2);
    theta(loopi)=asin((nc(en(loopi,1),2)-nc(en(loopi,2),2))/el(loopi));
    lmd=[cos(theta(loopi)) sin(theta(loopi)); -sin(theta(loopi)) cos(theta(loopi))]; 
    t=[lmd zeros(2); zeros(2) lmd];
    dk=t'*ek*t/el(loopi);
    dm=t'*em*t*el(loopi);
    for jx=1:4
        for jy=1:4
            if(e2s(jx)*e2s(jy)~=0)
                k(e2s(jx),e2s(jy))=k(e2s(jx),e2s(jy))+dk(jx,jy);
                m(e2s(jx),e2s(jy))=m(e2s(jx),e2s(jy))+dm(jx,jy);
            end
        end
    end
end
load('MR1.mat', 'MR');
load('KR1.mat', 'KR');

k = KR;
m = MR;

load('MR2.mat');
load('KR2.mat');
k2 = KR2;
m2 = MR2;

load('MR3.mat');
load('KR3.mat');
k3 = KR;
m3 = MR;

load('K.mat');
load('M.mat');
k0 = kk1;
m0 = mm1;




p=140;
epsilon=1e-5;
[v1,d1]=Ssiter(k,m,p,epsilon);
frequency1=sqrt(d1)/(2*pi);

[v2,d2]=Ssiter(k2,m2,p,epsilon);
frequency2=sqrt(d2)/(2*pi);

[v3,d3]=Ssiter(k3,m3,p,epsilon);
frequency3=sqrt(d3)/(2*pi);

[v0,d0]=Ssiter(k0,m0,p,epsilon);
frequency0=sqrt(d0)/(2*pi);

[frequency0,indexf]=sort(frequency0);
frequency0=real(frequency0);
plot(frequency0);hold on;
plot(frequency1);hold on;
plot(frequency2);hold on;
plot(frequency3);hold on;