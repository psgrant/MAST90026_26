function [uu, x, y] = FDM2DMOD(f, a, b, c, d, n, m)

% Define the grid and step sizes
hx = (b-a)/(n+1); % step size in x direction or \Detla x in ppt
hy = (d-c)/(m+1); % step size in y direcito or \Delta y in ppt
xx = a:hx:b;
yy = c:hy:d;


% Sparse Matrix Setup for FDM
e = ones(n,1);
D = spdiags(-1/hy^2*e, 0, n, n);
C = spdiags([-1/hx^2*e 2*(1/hx^2+1/hy^2)*e -1/hx^2*e], [-1 0 1], n, n);
r = ones(m,1);
I = spdiags(r, 0, m, m);
E = spdiags([r r], [-1, 1], m, m);
A = kron(I, C) + kron(E, D);
Q=zeros(n*m,n*m); 

% form the right hand side
[x, y] = meshgrid(xx, yy);
rhs = zeros(n*m,1);
for i=1:m
    for j=(i-1)*n+1:(i-1)*n+n
        Q(j,j)=(x(i+1,j-(i-1)*n+1)*x(i+1,j-(i-1)*n+1)+y(i+1,j-(i-1)*n+1)*y(i+1,j-(i-1)*n+1));
    end
end
A=A+Q;
for j = 1:m
    for i=1:n
    rhs((j-1)*n+i) = f(x(j+1,i+1), y(j+1,i+1));
    end
    rhs((j-1)*n+1)=rhs((j-1)*n+1)+exp(x(j+1,1)^2+y(j+1,1)^2)/(hx)^2;
    rhs((j-1)*n+n)=rhs((j-1)*n+n)+exp(x(j+1,n+2)^2+y(j+1,n+2)^2)/(hx)^2;
end
b1=zeros(n,1);
c1=zeros(n,1);
for i=1:n
    b1(i)=rhs(i);
    c1(i)=exp(x(1,i+1)^2+y(1,i+1)^2);
end
b1=b1-D*c1;
for i=1:n
    rhs(i)=b1(i);
end

b2=zeros(n,1);
c2=zeros(n,1);
for i=n*m-n+1:n*m
    b2(i-(m-1)*n)=rhs(i);
end
for i=1:n
    c2(i)=exp(x(m+2,i+1)^2+y(m+2,i+1)^2);
end
b2=b2-D*c2;
for i=n*m-n+1:n*m
    rhs(i)=b2(i-(n*m-n));
end

% solve the linear system 
u = A\rhs;


% Assembling the Solution with Boundary Values
uu = zeros(n+2, m+2);
uu(2:end-1, 2:end-1) = reshape(u, n, m);
uu = uu';
for i=1:n+2
    uu(1,i)=exp(x(1,i)^2+y(1,i)^2);
    uu(m+2,i)=exp(x(m+2,i)^2+y(m+2,i)^2);
end
for i=2:m+1
    uu(i,1)=exp(x(i,1)^2+y(i,1)^2);
    uu(i,n+2)=exp(x(i,n+2)^2+y(i,n+2)^2);
end

% Plot the solution
surf(x, y, uu)