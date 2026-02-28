clear
close all
clc

% Boundary of rectangular domain [a, b]x[c, d]
a = -1;
b = 1;
c = 0;
d = 1;
% Number of points in the x and y directions
nx = 10;
ny = 20;
% Define which quadrature rule to use
quad_order = 2;
% List of x and y points
x = linspace(a, b, nx);
y = linspace(c, d, ny);
% Create grid of x and y values
[X, Y] = meshgrid(x, y);
% Turn grid into list we can use to represent node
X = X(:);
Y = Y(:);
node = [X, Y];
% Generate the triangulation from the nodes
elem = delaunay(node);

% Define the coefficient functions, boundary conditions and exact solution
DFunc = @(X)exp(X(1) + X(2));
qFunc = @(X)2*exp(X(1) + X(2));
fFunc = @(X)-2*exp(2*(X(1) + X(2)));
uExact = @(X)exp(X(:, 1) + X(:, 2));
gFunc = uExact;

% Solve the PDE using the finite element method
u = FEM_Elliptic_2D_Dirichlet(node, elem, DFunc, qFunc, fFunc, gFunc, quad_order);

% Plot the student result and analytical solution
subplot(1, 2, 1)
%trisurf(elem, node(:, 1), node(:, 2), u)
trimesh(elem, node(:, 1), node(:, 2), u)
title('Student FEM Solution')
xlabel('x')
ylabel('y')
pbaspect([1 1 1])
subplot(1, 2, 2)
%trisurf(elem, node(:, 1), node(:, 2), uExact(node))
trimesh(elem, node(:, 1), node(:, 2), uExact(node))
title('Exact Solution')
xlabel('x')
ylabel('y')
pbaspect([1 1 1])
