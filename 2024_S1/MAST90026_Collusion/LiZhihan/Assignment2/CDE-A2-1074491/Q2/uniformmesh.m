% Define the number of points in each dimension
n = 50;
m = 100;

% Create a uniform grid of points inside the unit square
x = linspace(0, 1, n);
y = linspace(0, 1, n);
x1 = linspace(0, 1, m);
y1 = linspace(0, 1, m);
[X, Y] = meshgrid(x, y);
[X1, Y1] = meshgrid(x1, y1);

% Reshape into column vectors
nodes = [x(:), y(:)];
nodes1 = [x1(:), y1(:)];

% Create the Delaunay triangulation
elem = delaunay(node);
elem1 = delaunay(node1);

quad_order=2; % The quad_order

appro=FEM_Elliptic_2D_Dirichlet(node, elem, @(x,y) 1, @(x,y) 0, @(x,y) 17*sin(x).*cos(4*y), @(x,y) sin(x).*cos(4*y), quad_order);
appro1=FEM_Elliptic_2D_Dirichlet(node1, elem1, @(x,y) 1, @(x,y) 0, @(x,y) 17*sin(x).*cos(4*y), @(x,y) sin(x).*cos(4*y), quad_order);
exact=zeros(n*n,1);
exact1=zeros(m*m,1);
for i=1:n^2 %% Exact solutions
    exact(i)=sin(node(i,1))*cos(4*node(i,2));
end
for i=1:m^2
    exact1(i)=sin(node1(i,1))*cos(4*node1(i,2));
end
max_err= norm(appro-exact,inf);
disp(max_err); % % Calculate maximum grid errors
max_err1= norm(appro1-exact1,inf);
p=log(max_err/max_err1)/log(2); %compute the order
disp(p)

