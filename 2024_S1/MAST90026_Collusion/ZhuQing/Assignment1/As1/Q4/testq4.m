a=-1;
b=1;
N=101;
p=1;
x=linspace(a,b,N);
%[~,x]=cheb(N-1);
%x= flip(x);
Dx = @(x) 1;
qx = @(x) 32*x^(6)/(p+x^4)^2;
fx = @(x) 12*x^(2)/(p+x^4)^2;
D=zeros(1,N);
q=zeros(1,N);
f=zeros(1,N);
for i=1:length(x)
    D(i) = Dx(x(i));
    q(i) = qx(x(i));
    f(i) = fx(x(i));
end
alpha1=1;
alpha2=-1;
beta=0.5;
u=BvpFE(x,[linspace(1,N-1,N-1);linspace(2,N,N-1)]',D,q,f,alpha1,beta);
plot(x,u)