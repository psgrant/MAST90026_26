%After defining the functions, test them:
N=101;
[D,x]=cheb(N-1);
p=0.01;

%defines the functions px, qx, fx for the given BVP.
px = @(x) 0;
qx = @(x) -32*x^(6)/(p+x^4)^2;
fx = @(x) -12*x^(2)/(p+x^4)^2;

%initializes vectors pt, qt, ft and fills them with the evaluated function values at each Chebyshev point.
pt=zeros(1,length(x));
qt=zeros(1,length(x));
ft=zeros(1,length(x));
for i=1:length(x)
    pt(i)=  px(x(i));
    qt(i) = qx(x(i));
    ft(i) = fx(x(i));
end
alpha=1;%boundary conditions
beta=0.5;
N=length(x);
[x,u]=BvpSpecnew(D,pt,qt,ft,alpha,beta,N);%compute the solution u to the BVP at the Chebyshev points x.
plot(x,u)