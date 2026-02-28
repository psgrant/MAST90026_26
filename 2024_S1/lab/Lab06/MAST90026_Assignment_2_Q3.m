% Boundary of rectangular domain [a, b]x[c, d]
a = 0;
b = 1;
c = 0;
d = 1;
% Number of points in the x and y directions
nx = 20;
ny = 20;
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

DFunc = @(X)1;
fFunc = @(X)exp(X(1) + 2*X(2))*(exp(X(1) + 2*X(2)) - 5);
uExact = @(X)exp(X(:, 1) + 2*X(:, 2));
gFunc = uExact;

max_iter = 100;
tol = 1e-10;
uk = 1*ones(size(node, 1), 1);
for i = 1:max_iter
    ukInterpolate = scatteredInterpolant(node, uk);
    ukFunc = @(X)ukInterpolate([X(1) X(2)]);
    fkFunc = @(X)(fFunc(X) + ukFunc(X).^2);
    qkFunc = @(X)2*ukFunc(X);
    
    uk1 = FEM_Elliptic_2D_Dirichlet(node, elem, DFunc, qkFunc, fkFunc, gFunc);

    errk = max(abs(uk1 - uk));
    if errk < tol
        break
    else
        uk = uk1;
    end
end

err = max(abs(uk1 - uExact(node)));
trimesh(elem, node(:, 1), node(:, 2), uk1)