function [uu, x, y] = FDM2D(f, a, b, c, d, n, m)

hx = (b-a)/(n+1); % step size in x direction or \Detla x in ppt
hy = (d-c)/(m+1); % step size in y direcito or \Delta y in ppt

xx = a:hx:b;
yy = c:hy:d;

% form the matrix
e = ones(n,1);
D = spdiags(-1/hy^2*e, 0, n, n);
C = spdiags([-1/hx^2*e 2*(1/hx^2+1/hy^2)*e -1/hx^2*e], [-1 0 1], n, n);
r = ones(m,1);
I = spdiags(r, 0, m, m);
E = spdiags([r r], [-1, 1], m, m);
A = kron(I, C) + kron(E, D);

% form the right hand side
[x, y] = meshgrid(xx, yy);
rhs = zeros(n*m,1);
for j = 1:m
    rhs((j-1)*n+1:(j-1)*n+n) = f(x(j+1,2:n+1), y(j+1,2:n+1));
end
%fij = f(x(2:end-1, 2:end-1), y(2:end-1, 2:end-1));
%rhs = fij(:);

% solve the linear system 
u = A\rhs;

uu = zeros(n+2, m+2);
uu(2:end-1, 2:end-1) = reshape(u, n, m);
uu = uu';

surf(x, y, uu)