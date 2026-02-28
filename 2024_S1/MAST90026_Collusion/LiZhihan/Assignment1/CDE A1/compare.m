N=500;
p=0.01;
s=zeros(N-2,1);
for j=3:N
    x=linspace(-1,1,j);
    %[~,x]=cheb(N-1);
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
w=zeros(N-2,1);
for i=3:N
    w(i-2)=log(i);
end
z=zeros(N-2,1);
for i=3:N
    z(i-2)=log(s(i-2));
end
snew=zeros(N-2,1);
for j=3:N
    [~,xnew]=cheb(j-1);
    xnew=flip(xnew);
    Dx = @(x) 1;
    qx = @(x) 32*x^(6)/(p+x^4)^2;
    fx = @(x) 12*x^(2)/(p+x^4)^2;
    D=zeros(1,j);
    q=zeros(1,j);
    f=zeros(1,j);
    for i=1:length(xnew)
    D(i) = Dx(xnew(i));
    q(i) = qx(xnew(i));
    f(i) = fx(xnew(i));
    end
    alpha=4/((p+1)^2);
    beta=1/(1+p);
    s1=BvpFE(xnew,[linspace(1,j-1,j-1);linspace(2,j,j-1)]',D,q,f,alpha,beta);
    s2=TRUE(xnew,[linspace(1,j-1,j-1);linspace(2,j,j-1)]',p);
    snew(j-2)=norm(s1-s2,'inf');
end
wnew=zeros(N-2,1);
for i=3:N
    wnew(i-2)=log(i);
end
znew=zeros(N-2,1);
for i=3:N
    znew(i-2)=log(snew(i-2));
end

% Plotting both lines in one figure
plot(w, z, 'b-', wnew, znew, 'r--')
xlabel('log(N)')
ylabel('log(error)')
title('Convergence Behavior')
legend('linespace', 'Cheb', 'Location', 'best')


