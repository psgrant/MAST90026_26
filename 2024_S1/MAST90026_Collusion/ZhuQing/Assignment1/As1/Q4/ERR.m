function e=ERR(N,p)
s=zeros(N-2,1);
for i=3:N
    [D,x]=cheb(i-1);
px = @(x) 0;
qx = @(x) -32*x^(6)/(p+x^4)^2;
fx = @(x) -12*x^(2)/(p+x^4)^2;
pt=zeros(1,length(x));
qt=zeros(1,length(x));
ft=zeros(1,length(x));
    for j=1:length(x)
    pt(j)=  px(x(j));
    qt(j) = qx(x(j));
    ft(j) = fx(x(j));
    end
alpha=4/((p+1)^2);
beta=1/(1+p);
    [~,s1]=BvpSpec(D,pt,qt,ft,alpha,beta,length(x));
    s2=TRUE2(i,p);
    s(i-2)=norm(s2-s1,'inf');
end
e=s;
x=zeros(N-2,1);
y=zeros(N-2,1);
for i=3:N
    x(i-2)=log(i);
    y(i-2)=log(e(i-2));
end
plot(x,y)