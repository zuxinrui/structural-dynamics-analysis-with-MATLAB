q=min(2*p,p+8);        
[n1,n2]=size(K);
v0=zeros(n1,q);
MK=zeros(n1,1);
d0=zeros(p,1);
for i=1:n1
    MK(i)=M(i,i)/K(i,i);
end
[Y, I]=sort(MK,1,'descend');
v0(:,1)=1;
for i=2:q
    v0(I(i-1),i)=1/Y(i-1);
end

iter=0;
flag=1;
while flag
    iter=iter+1;
    yiter=M*v0;
    xiter1=K\yiter;
    Km=(xiter1)'*K*xiter1;
    Mm=(xiter1)'*M*xiter1;
    [vm,dm]=eig(Km,Mm);   
    [lambda,k]=sort(diag(dm));
    vm=vm(:,k);
    Factor=diag(vm'*Mm*vm);
    Vnorm=vm*inv(sqrt(diag(Factor)));               
    v0=xiter1*Vnorm;
    diter=lambda(1:p);
    if abs(norm(diter)-norm(d0))/norm(diter)<epsilon      
        flag=0;
    else
        flag=1;
    end
    d0=diter;
end
v=v0(:,1:p);
d=d0;
