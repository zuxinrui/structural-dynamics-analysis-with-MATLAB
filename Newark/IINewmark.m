E=2.0e11;
A=2.5e-3;
density=7860; 
node_number=24;
elment_number=46;
nc=[0 0;2 0;4 0;6 0;8 0;10 0;0 2;2 2;4 2;6 2;8 2;10 2;2 4;4 4;6 4;8 4;0 8;2 8;4 8;6 8;8 8;10 8;4 12;6 12];
en=[1 7;1 8;1 2;2 8;2 3;3 9;3 4;4 10;4 5;5 11;6 11;5 6;6 12;11 12;11 16;10 11;10 16;10 15;9 10;9 14;9 13;8 9;8 13;7 8;13 17;13 18;13 14;14 18;14 19;14 15;15 20;15 21;15 16;16 21;16 22;21 22;20 21;21 24;20 24;20 23;19 20;19 23;18 19;17 18;18 23;23 24];
ed(1:node_number,1:2)=1;
dof=0;
for loopi=1:node_number
     for loopj=1:2
            dof=dof+1;
            ed(loopi,loopj)=dof;
    end
end
ek=E*A*[1 0 -1 0;0 0 0 0;-1 0 1 0;0 0 0 0];
em=(density*A)/6*[2 0 1 0;            
                   0 2 0 1;
                  1 0 2 0;
                  0 1 0 2];   
k(1:dof,1:dof)=0; 
m=k; 
theta(1:9)=0;
el(1:9)=0;
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
                k(e2s(jx),e2s(jy))=k(e2s(jx),e2s(jy))+dk(jx,jy);
                m(e2s(jx),e2s(jy))=m(e2s(jx),e2s(jy))+dm(jx,jy);
        end
    end
end
c=0*k+0*m;
nt=10000;dt=0.0001;
time=0:dt:nt*dt;

q0=zeros(dof,1);
dq0=zeros(dof,1);
bcdof=[1;1;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
fd=zeros(dof,nt);
for i=1:nt
    fd(17,i)=-200*cos(100*i);
end
[acc,vel,dsp]=Newmark(k,c,m,fd,bcdof,nt,dt,q0,dq0);
plot(time, dsp(35,:))  
 xlabel('Time(seconds)')
 ylabel('Horizontal displ.(m)')
 title('alpha=1;beta=0.5;nt=10000;dt=0.0001')