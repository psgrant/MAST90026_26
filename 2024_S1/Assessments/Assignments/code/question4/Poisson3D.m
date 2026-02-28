% Poisson3D
%   sets up the matrix for a 3D Poisson problem
%   on the unit cube with homogeneous Dirichlet BC
%   n X n X n internal meshpoints

% adapted from mit18086_Poisson.m by Benjamin Seibold

n=80;
K1D = spdiags(ones(n,1)*[-1 2 -1],-1:1,n,n);  % 1d Poisson matrix
I1D = speye(size(K1D));                       % 1d identity matrix
K2D = kron(K1D,I1D)+kron(I1D,K1D);            % 2d Poisson matrix
I2D = speye(size(K2D));                       % 2d identity matrix
K3D = kron(K2D,I1D)+kron(I2D,K1D);            % 3d Poisson matrix

xtrue = rand(n^3,1); 
b = K3D*xtrue; 
x = K3D\b;