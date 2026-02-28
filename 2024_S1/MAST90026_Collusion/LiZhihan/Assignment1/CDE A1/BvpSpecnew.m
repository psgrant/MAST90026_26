%solves the BVP using spectral collocation at Chebyshev points.
function [x,u]=BvpSpecnew(D,pt,qt,ft,alpha,beta,N)
[~,x]=cheb(N-1);

%matrices related to the second derivative：
W=D*D;
P=diag(ones(N,1));%diagonal matrices
Q=diag(ones(N,1));
for i=1:N
    P(i,i)=pt(i);
    Q(i,i)=qt(i);
end
A=W+P*D+Q;%constructed using these matrices and adding them up.
for i=1:N
    A(N,i)=D(N,i); %overwrites the last row of A with the last row of the differentiation matrix D
end
rhs=ones(N,1);%initialized and filled with values from ft.
s=zeros(N,1);
s(1)=beta;%boundary condition
for i=1:N
    rhs(i)=ft(i);
end
rhs(N)=alpha;%boundary condition
rhs=rhs-A*s;
s(2:N)=A(2:N,2:N)\rhs(2:N);
u=s;
u= flip(u);
plot(u)
return












