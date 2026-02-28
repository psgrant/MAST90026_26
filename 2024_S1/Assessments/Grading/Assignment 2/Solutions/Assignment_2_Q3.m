% Boundary of rectangular domain [a, b]x[c, d]
a = 0;
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

% Initial guess
uk = ones(nx*ny, 1);
% Maximum iterations before terminating
MAX_ITER = 100;
% Tolerance on the inf-norm before ending iteration
TOL = 1e-6;

% Define functions which don't change on each iteration
DFunc = @(X)1;
uExact = @(X)exp(X(:, 1) + 2*X(:, 2));
gFunc = uExact;
for j = 1:MAX_ITER
    % Define the coefficient functions which change on each iteration
    ukFunc = scatteredInterpolant(node(:, 1), node(:, 2), uk);
    qFunc = @(X)2*ukFunc(X(1), X(2));
    fFunc = @(X)exp(X(1) + 2*X(2))*(exp(X(1) + 2*X(2)) - 5) + ukFunc(X(1), X(2))^2;
    
    uk1 = FEM_Elliptic_2D_Dirichlet(node, elem, DFunc, qFunc, fFunc, gFunc, quad_order);

    err = max(abs(uk1 - uk));

    uk = uk1;
    if err < TOL
        disp(strcat('Returned solution is converged: ||E|| =  ', num2str(err)))
        break
    end
end

if j == MAX_ITER
    disp(strcat('Returned solution is NOT converged: ||E|| = ', num2str(err)))
end

% Plot the student result and analytical solution
subplot(1, 2, 1)
%trisurf(elem, node(:, 1), node(:, 2), u)
trimesh(elem, node(:, 1), node(:, 2), uk)
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
