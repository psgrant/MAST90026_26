function [order,e]=loglogg2(N,p)
s=zeros(N-2,1);
for j=3:N
    [~,x]=cheb(j-1);
    x=flip(x);
    Dx = @(x) 1;
    qx = @(x) 32*x^(6)/(p+x^4)^2;
    fx = @(x) 12*x^(2)/(p+x^4)^2;
    D=zeros(1,j);
    q=zeros(1,j);
    f=zeros(1,j);
    for i=1:length(x)
    D(i) = Dx(x(i));
    q(i) = qx(x(i));
    f(i) = fx(x(i));
    end
    alpha=4/((p+1)^2);
    beta=1/(1+p);
    s1=BvpFE(x,[linspace(1,j-1,j-1);linspace(2,j,j-1)]',D,q,f,alpha,beta);
    s2=TRUE(x,[linspace(1,j-1,j-1);linspace(2,j,j-1)]',p);
    s(j-2)=norm(s1-s2,'inf');
end
e=s;
w=zeros(N-2,1);
for i=3:N
    w(i-2)=log(i);
end
z=zeros(N-2,1);
for i=3:N
    z(i-2)=log(s(i-2));
end
k = polyfit(w, z, 1);
order=-k(1);
plot(w,z)