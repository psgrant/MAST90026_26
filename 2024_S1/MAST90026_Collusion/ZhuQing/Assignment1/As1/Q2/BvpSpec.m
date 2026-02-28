function [x,u]=BvpSpec(D,pt,qt,ft,alpha,beta,N)
[~,x]=cheb(N-1);
W=D*D;
P=diag(ones(N,1));
Q=diag(ones(N,1));
for i=1:N
    P(i,i)=pt(i);
    Q(i,i)=qt(i);
end
A=W+P*D+Q;
for i=1:N
    A(N,i)=D(N,i);
end
rhs=ones(N,1);
s=zeros(N,1);
s(1)=beta;
for i=1:N
    rhs(i)=ft(i);
end
rhs(N)=alpha;
rhs=rhs-A*s;
s(2:N)=A(2:N,2:N)\rhs(2:N);
u=s;
u= flip(u);













