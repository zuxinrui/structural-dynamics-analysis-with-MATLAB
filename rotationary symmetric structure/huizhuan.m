clear; clc;
asd = sqrt(3);

Prob_title='static analysis of a 2-D frame structure';

No_el=54;                                             
No_nel=2;                                      
No_dof=2;                                        
No_node=30;                                   
Sys_dof=No_node*No_dof;                               
prop(1)=2.1e11;                                         
prop(2)=0.3;                                                
prop(3)=7860;                              
prop(4)=0.05;                       
prop(5)=0.05;                            
prop(6)=0.05*0.05;                          
prop(7)=0.05*0.05^3/12;              
prop(8)=0.05*0.05^3/12;              
 
prop(9)=1;                                            
                   
prop(10)=1;                                                
                                   
Opt_section=1;                         

gcoord=[
        1    0.0      2.0      0.0;
        2    asd      1.0      0.0;
        3    asd      3.0      0.0;
        4    2*asd      2.0      0.0;
        5    2*asd      4.0      0.0;
        6    3*asd      3.0      0.0;
        7    asd      -1.0      0.0;
        8    2*asd      0.0      0.0;
       9    2*asd      -2.0      0.0;
       10    3*asd      -1.0      0.0;
       11    3*asd      -3.0      0.0;
       12    0.0      -2.0      0.0;
       13    asd      -3.0      0.0;
       14    0.0      -4.0      0.0;
       15    asd      -5.0      0.0;
       16    0.0      -6.0      0.0;
       17    -asd      -1.0      0.0;
       18    -asd      -3.0      0.0;
       19    -2*asd      -2.0      0.0;
       20    -2*asd      -4.0      0.0;
       21    -3*asd      -3.0      0.0;
       22   -asd   1.0    0.0;
       23   -2*asd    0.0    0.0;
       24   -2*asd    2.0    0.0;
       25   -3*asd    1.0    0.0;
       26   -3*asd    3.0    0.0;
       27   -asd    3.0    0.0;
       28   0.0    4.0    0.0;
       29   -asd    5.0    0.0;
       30   0.0     6.0    0.0]
gcoord=[gcoord(:,2:4)];

nodes=[ 1          1          2;
        2          1          3;
        3          2          3;
        4          2          4;
        5          3          4;
        6          3          5;
        7          4          5;
        8          4          6;
        9          5         6;
       10         2         7;
       11         2         8;
       12   7   8;
       13   7   9;
       14   8   9;
       15   8   10;
       16   9   10;
       17   9   11;
       18   10   11;
       19   7   12;
       20   7   13;
       21   12   13;
       22   12   14;
       23   13   14;
       24   13   15;
       25   14   15;
       26   14  16;
       27   15  16;
       28   12  17;
       29   12   18;
       30   17   18;
       31   17   19;
       32   18   19;
       33   18   20;
       34   19   20;
       35   19  21;
       36   20  21;
       37   17  22;
       38   17  23;
       39   22  23;
       40   22  24;
       41   23  24;
       42   23  25;
       43   24  25;
       44   24  26;
       45   25  26;
       46   1  22;
       47   22  27;
       48   1  27;
       49   1  28;
       50   27  28;
       51   27  29;
       52   28  29;
       53   28  30;
       54   29  30
       ];
nodes=[nodes(:,2:3)];
P=[-9,1,1];
ConNode=[ 6,  1,  1,  1;
          11,  1,  1,  1;
          16,  1,  1,  1;
          21,  1,  1,  1;
          26,  1,  1,  1;
          30,  1,  1,  1
          ];
ConVal =[ 6,  0,  0,  0;
          11,  0,  0,  0;
          16,  0,  0,  0;
          21,  0,  0,  0;
          26,  0,  0,  0;
          30,  0,  0,  0
          ];                           

Opt_beam=1;                                        
                                                      
Opt_mass=2;                                            
Opt_section=1;                                   
Opt_graphics1=1;                       
Opt_graphics2=1;                       

k=zeros(No_nel*No_dof,No_nel*No_dof);                   
m=zeros(No_nel*No_dof,No_nel*No_dof);                    

kk=zeros(Sys_dof,Sys_dof);                  
mm=zeros(Sys_dof,Sys_dof);                  
ff=zeros(Sys_dof,1);                           

bcdof=zeros(Sys_dof,1);                               
bcval=zeros(Sys_dof,1);                               

index=zeros(No_nel*No_dof,1);                        
[n1,n2]=size(ConNode);
                                                 
for ni=1:n1
  ki=ConNode(ni,1);
  bcdof((ki-1)*No_dof+1:ki*No_dof)=ConNode(ni,2:No_dof+1);
                                                 
  bcval((ki-1)*No_dof+1:ki*No_dof)=ConVal(ni,2:No_dof+1);
                                                
end

[n1,n2]=size(P);

for ni=1:n1
  ff(No_dof*(P(ni,2)-1)+P(ni,3))=P(ni,1);
end

for iel=1:No_el                              
  nd(1)=nodes(iel,1);                     
  nd(2)=nodes(iel,2);                   
  x(1)=gcoord(nd(1),1); y(1)=gcoord(nd(1),2); z(1)=gcoord(nd(1),3);
                                                      
  x(2)=gcoord(nd(2),1); y(2)=gcoord(nd(2),2); z(2)=gcoord(nd(2),3);
                                                     
  leng=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2+(z(2)-z(1))^2);
                                                       
  if (x(2)-x(1))==0;
    beta=pi/2; 
  elseif ((y(2)-y(1))==0)&((x(2)-x(1))<0)
    beta=pi;
  else
    beta=atan((y(2)-y(1))/(x(2)-x(1)));
  end 

  if Opt_beam==1
    [k,m]=F21(prop,leng,beta,Opt_section,Opt_mass);
                                
  elseif Opt_beam==2
    [k,m]=F22(prop,leng,beta,Opt_section,Opt_mass);
                                  
  end

  index=dof(nd,No_nel,No_dof);
                                     
  kk=semble(kk,k,index);                   % assemble system stiffness matrix
  mm=semble(mm,m,index);                   % assemble system mass matrix
end
[kk1,mm1,ff1]=bc(kk,mm,ff,bcdof,bcval);
                            
displmt=kk1\ff1;                 
sdof1=60;
No_dof1=No_dof/(Sys_dof/sdof1);
for ii=1:No_node
  for ij=1:No_dof1
    displmtnode(ii,ij)=displmt((ii-1)*No_dof1+ij,1);
  end
end

if Opt_graphics1==1
  for iel=1:No_el                            
    nd(1)=nodes(iel,1);                  
    nd(2)=nodes(iel,2);                 
    x(1)=gcoord(nd(1),1); y(1)=gcoord(nd(1),2); z(1)=gcoord(nd(1),3);
                                                     
    x(2)=gcoord(nd(2),1); y(2)=gcoord(nd(2),2); z(2)=gcoord(nd(2),3);
                                                     
    figure(1)
      plot(x,y), xlabel('x'), ylabel('y'), hold on;
      axis([-12 12 -12 12]);
  end

  if Opt_graphics2~=1
    title('nodal connectivity of elements'), hold off;
  end
end

if Opt_graphics2==1
  Ampl=10;
  gcoordA=gcoord+[Ampl*displmtnode(:,1:2),zeros(No_node,1)];
  for iel=1:No_el                            
    nd(1)=nodes(iel,1);                  
    nd(2)=nodes(iel,2);                
    x(1)=gcoordA(nd(1),1); y(1)=gcoordA(nd(1),2); z(1)=gcoordA(nd(1),3);
                                                       
    x(2)=gcoordA(nd(2),1); y(2)=gcoordA(nd(2),2); z(2)=gcoordA(nd(2),3);
                                                      
    figure(1)
      plot(x,y,'r'), title('structure'), hold on;
  end
  hold off
end

disp('The calculation is use of')

if Opt_beam==1
  disp('Euler-Bernoulli beam element')
elseif Opt_beam==2
  disp('Timoshenko beam element')
end
disp(' ')
disp('                  numrical solution')
disp('     node        x        y          ')
num=1:1:No_node;
displmts=[num' displmtnode]                            
