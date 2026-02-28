function u=BvpFE(node,elem,D,q,f,alpha,beta)
n=length(node);
x=zeros(n,1);
for i=1:n-1
x(i) = node(elem(i,1));
x(i+1) = node(elem(i,2));
end
A=zeros(n,n);
for i=2:n-1
    A(i,i)=(1/(x(i)-x(i-1)))*(D(i)+D(i-1))/2 ...
    +(1/(x(i+1)-x(i)))*(D(i+1)+D(i))/2;
end
A(1,1)=(1/(x(2)-x(1)))*(D(2)+D(1))/2;
A(n,n)=(1/(x(n)-x(n-1)))*(D(n)+D(n-1))/2;
for j=1:n-1
    A(j+1,j)=-(1/(x(j+1)-x(j)))*(D(j+1)+D(j))/2;
end

for j=2:n
    A(j-1,j)=-(1/(x(j)-x(j-1)))*(D(j)+D(j-1))/2;
end
RHS=zeros(n,1);
for i=2:n-1
    RHS(i)=f(i)*(x(i+1)-x(i-1))/2;
end
RHS(1)=f(1)*(x(2)-x(1))/2;
RHS(n)=f(n)*(x(n)-x(n-1))/2;
P=zeros(n,n);
for i=2:n-1
    P(i,i)=q(i)*(x(i+1)-x(i-1))/2;
end
P(1,1)=q(1)*(x(2)-x(1))/2;
P(n,n)=q(n)*(x(n)-x(n-1))/2;
RHS(1)=RHS(1)-alpha*D(1);
RHS(n)=beta;
LHS=A+P;
LHS(end,:)=0;
LHS(end,end)=1;
u=LHS\RHS;














