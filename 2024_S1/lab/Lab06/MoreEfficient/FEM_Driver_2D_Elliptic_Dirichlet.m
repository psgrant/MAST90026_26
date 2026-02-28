% Boundary of rectangular domain [a, b]x[c, d]
xS = -1;
xE = 1;
yS = -1;
yE = 1;
% Number of points in the x and y directions
nx = 20;
ny = 20;

[node, elem, bndNode, freeNode] = createTriMeshSquare(xS, xE, yS, yE, nx, ny);

% Define the coefficient functions, boundary conditions and exact solution
DFunc = @(X)1;
qFunc = @(X)1;
fFunc = @(X)-exp(X(1) + X(2));
uExact = @(X)exp(X(:, 1) + X(:, 2));

bndVal = uExact(node(bndNode, :));

% Solve the PDE using the finite element method
u = FEM_Elliptic_2D_Dirichlet(node, elem, DFunc, qFunc, fFunc, bndNode, bndVal, freeNode);

% Calculate the error using triangle quadrature on the elements
err = calcErr(node, elem, u, uExact);

% Plot the mesh and the results
subplot(1, 3, 1)
triplot(elem, node(:, 1), node(:, 2))
pbaspect([1 1 1])
title('Mesh')
xlabel('x')
ylabel('y')
subplot(1, 3, 2)
%trisurf(elem, node(:, 1), node(:, 2), u)
trimesh(elem, node(:, 1), node(:, 2), u)
title('FEM Solution')
xlabel('x')
ylabel('y')
pbaspect([1 1 1])
subplot(1, 3, 3)
%trisurf(elem, node(:, 1), node(:, 2), u)
trimesh(elem, node(:, 1), node(:, 2), uExact(node))
title('Exact Solution')
xlabel('x')
ylabel('y')
pbaspect([1 1 1])