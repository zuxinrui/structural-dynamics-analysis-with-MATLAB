function [acc,vel,dsp]=Newmark(kk,cc,mm,fd,bcdof,nt,dt,q0,dq0)
[sdof,n2]=size(kk);
 
dsp=zeros(sdof,nt);
vel=zeros(sdof,nt);
acc=zeros(sdof,nt);
dsp(:,1)=q0; 
vel(:,1)=dq0;

alpha=1; beta=0.5;
acc(:,1)=inv(mm)*(fd(:,1)-kk*dsp(:,1)-cc*vel(:,1));
ekk=kk+mm/(alpha*dt^2)+cc*beta/(alpha*dt);
for i=1:sdof
  if bcdof(i)==1
    dsp(i,1)=0;
    vel(i,1)=0;
    acc(i,1)=0;
  end
end

for it=1:nt
  cfm=dsp(:,it)/(alpha*dt^2)+vel(:,it)/(alpha*dt)+acc(:,it)*(0.5/alpha-1);
  cfc=dsp(:,it)*beta/(alpha*dt)+vel(:,it)*(beta/alpha-1)...
     +acc(:,it)*(0.5*beta/alpha-1)*dt;
  efd=fd(:,it)+mm*cfm+cc*cfc;

  dsp(:,it+1)=inv(ekk)*efd;
  acc(:,it+1)=(dsp(:,it+1)-dsp(:,it))/(alpha*dt^2)-vel(:,it)/(alpha*dt)...
              -acc(:,it)*(0.5/alpha-1);
  vel(:,it+1)=vel(:,it)+acc(:,it)*(1-beta)*dt+acc(:,it+1)*beta*dt;
  for i=1:sdof
    if bcdof(i)==1
      dsp(i,it+1)=0;
      vel(i,it+1)=0;
      acc(i,it+1)=0;
    end
  end

end

if cc(1,1)==0
  disp('The transient response results of undamping system')
else
  disp('The transient response results of damping system')
end
 
disp('The method is Newmark integration')
