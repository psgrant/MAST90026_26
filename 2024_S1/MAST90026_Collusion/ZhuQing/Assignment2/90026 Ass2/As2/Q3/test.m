% Define the number of points in each dimension
n = 100;

% Create a uniform grid of points inside the unit square
x = linspace(0, 1, n);
y = linspace(0, 1, n);
[X, Y] = meshgrid(x, y);

% Reshape into column vectors
x = X(:);
y = Y(:);
node = [x, y];
% Create the Delaunay triangulation
elem = delaunay(node);
N=length(node(:,1));
exact=zeros(N,1);
for i=1:N
    exact(i)=exp(node(i,1)+2*node(i,2));
end
[it,appro]=solvenew(@(x,y) -5*exp(x+2*y)+exp(2*x+4*y),@(x,y) exp(x+2*y),@(x,y) 1,elem,node);
error = norm(appro-exact,inf);
disp(error)
disp(it)
