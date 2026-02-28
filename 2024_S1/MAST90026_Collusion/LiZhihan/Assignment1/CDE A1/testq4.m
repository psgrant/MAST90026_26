%test the function BvpFEnew:
a=-1;
b=1;
N=201;
p=1;
x=linspace(a,b,N);
%defines the functions Dx, qx, and fx, which are the values of D(x), q(x), and f(x), respectively.
Dx = @(x) 1;
qx = @(x) 32*x^(6)/(p+x^4)^2;
fx = @(x) 12*x^(2)/(p+x^4)^2;
%initializes the vectors D, q, and f and fills them with the evaluated function values at the points x.
D=zeros(1,N);
q=zeros(1,N);
f=zeros(1,N);
for i=1:length(x)
    D(i) = Dx(x(i));
    q(i) = qx(x(i));
    f(i) = fx(x(i));
end
%boundary conditions
alpha=1;
beta=0.5;
u2=BvpFEnew(linspace(-1,1,N),[linspace(1,N-1,N-1);linspace(2,N,N-1)]',D,q,f,alpha,beta);