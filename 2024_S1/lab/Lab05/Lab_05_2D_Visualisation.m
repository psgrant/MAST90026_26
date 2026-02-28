% Number of points and domain
N = 5;
xs = -1;
xe = 1;
ys = -1;
ye = 1;

% Vectors of x and y values
x = linspace(xs, xe, N+2);
y = linspace(ys, ye, N+2);

% Create matrices of x and y values on a grid
[X, Y] = meshgrid(x, y);

% Create function and evaluate it on the 2D Mesh
u_exact_fun = @(x, y)sin(pi*x).*cos(pi*y/2);
u_exact = u_exact_fun(X, Y);

% Create plots using a variety of ways to repersent 2D data
subplot(2, 2, 1)
surf(X, Y, u_exact)
xlabel('x')
ylabel('y')
title('Surface plot of exact solution', 'Grid displayed')
subplot(2, 2, 2)
surf(X, Y, u_exact, 'MeshStyle', 'none')
xlabel('x')
ylabel('y')
title('Surface plot of exact solution', 'Without grid')
subplot(2, 2, 3)
contourf(X, Y, u_exact)
colorbar
pbaspect([1 1 1])
xlabel('x')
ylabel('y')
title('Contour plot of exact solution', 'With colorbar')
subplot(2, 2, 4)
contour(X, Y, u_exact, 'ShowText', 'on')
pbaspect([1 1 1])
xlabel('x')
ylabel('y')
title('Contour plot of exact solution', 'With values on contours')