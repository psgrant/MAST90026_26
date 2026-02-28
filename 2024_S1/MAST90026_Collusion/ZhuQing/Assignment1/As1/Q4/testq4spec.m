N=100;
[D,x]=cheb(N-1);
p=1;
px = @(x) 0;
qx = @(x) -32*x^(6)/(p+x^4)^2;
fx = @(x) -12*x^(2)/(p+x^4)^2;
pt=zeros(1,length(x));
qt=zeros(1,length(x));
ft=zeros(1,length(x));
for i=1:length(x)
    pt(i)=  px(x(i));
    qt(i) = qx(x(i));
    ft(i) = fx(x(i));
end
alpha=4/((p+1)^2);
beta=1/(1+p);
N=length(x);
[x,u]=BvpSpec(D,pt,qt,ft,alpha,beta,N);
plot(x,u)