function [h, err, x, u] =  Homework01_2bii(N, eps)

% Start and end coordinates
xs = 0;
xe = 1;

% Boundary conditions
us = 1;
ue = 3;

% create vector of x values
x = linspace(xs, xe, N+2)';

h = (xe-xs)/(N+1);

% Lower diagonal coefficient
lC = eps/(h^2) + 1/h;
% Diagonal coefficient
dC = -(1/h + 2*eps/(h^2));
% Upper diagonal coefficient
uC = eps/(h^2);

% Create the FD matrix
e = ones(N, 1);
A = spdiags([lC*e, dC*e, uC*e], -1:1, N, N);

% Create RHS vector
b = -ones(N, 1);
b(1) = b(1) - lC*us;
b(N) = b(N) - uC*ue;

% Solve the system and append boundary conditions
u = A\b;
u = [us; u; ue];

% Calculate the grid 2-norm
err = sqrt(h*sum((u - Homework01_ExactSolution(x, eps)).^2));