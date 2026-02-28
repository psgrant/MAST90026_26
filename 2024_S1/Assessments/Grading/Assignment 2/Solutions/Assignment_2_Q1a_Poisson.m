n = 10;
m = 20;

xs = 0;
xe = 1;
ys = 0;
ye = 2;

x = linspace(xs, xe, n+2)';
y = linspace(ys, ye, m+2)';
[X, Y] = meshgrid(x, y);
    
Xin = X(2:m+1, 2:n+1)';
Yin = Y(2:m+1, 2:n+1)';

uExact = @(xa, ya)exp(xa + ya);
gFunc = uExact;

hx = (xe - xs)/(n + 1);
hy = (ye - ys)/(m + 1);

em = ones(m, 1);
en = ones(n, 1);
Cmat = spdiags([-en/hx^2, 2*en*(1/hx^2 + 1/hy^2), -en/hx^2], -1:1, n, n);
Dmat = -(1/hy^2)*speye(n);
Emat = spdiags([em, em], [-1,1], m, m);

A = kron(speye(m), Cmat) + kron(Emat, Dmat);
f = -2*exp(Xin(:) + Yin(:));

f(1:n) = f(1:n) - Dmat*gFunc(x(2:n+1), ys);
for i = 1:m
    f((i - 1)*n + 1) = f((i - 1)*n + 1) + gFunc(xs, y(i + 1))/(hx^2);
    f(i*n) = f(i*n) + gFunc(xe, y(i + 1))/(hx^2);
end
f((m-1)*n + 1:m*n) = f((m-1)*n + 1:m*n) - Dmat*gFunc(x(2:n+1), ye);

u = A\f;
u = reshape(u, [n, m])';

uTot = zeros(m+2, n+2);
uTot(1, :) = gFunc(x, ys);
uTot(end, :) = gFunc(x, ye);
uTot(:, 1) = gFunc(xs, y);
uTot(:, end) = gFunc(xe, y);
uTot(2:m+1, 2:n+1) = u;
surf(X, Y, uTot)
xlabel('x')
ylabel('y')
zlabel('u')

% Calculate the L2 norm
err = sqrt(trapz(y, trapz(x, (uTot - uExact(X, Y)).^2, 2)));