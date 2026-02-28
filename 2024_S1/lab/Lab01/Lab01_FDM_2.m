function [h, err, x, u] =  Lab01_FDM_2(N)

% Start and end coordinates
xs = -1;
xe = 1;

% Boundary conditions
us = 1;
ue = 2;

% create vector of x values
x = linspace(xs, xe, N+2)';

h = (xe-xs)/(N+1);



% Create indices and values of the matrix entries for use by sparse.
% R is row, C is column, diag is the diagonal entries, lowDiag is the lower
% diagonal entries and upDiag are the upper diagonal entries.
diagR = linspace(1, N, N)';
diagC = linspace(1, N, N)';
diag = -(abs(x(2:N+1)) + (2/(h^2))*ones(N, 1));

lowDiagR = linspace(2, N, N-1)';
lowDiagC = linspace(1, N-1, N-1)';
lowDiag = (1/(h^2) + 1/(2*h))*ones(N-1, 1);

upDiagR = linspace(1, N-1, N-1)';
upDiagC = linspace(2, N, N-1)';
upDiag = (1/(h^2) - 1/(2*h))*ones(N-1, 1);

% Collect all matrix entries into single vector
elemR = [diagR; lowDiagR; upDiagR];
elemC = [diagC; lowDiagC; upDiagC];
elem = [diag; lowDiag; upDiag];

% Create the FD matrix
A = sparse(elemR, elemC, elem);

% Create RHS vector
b = 30*x(2:N+1);
b(1) = b(1) - us/(h^2) - us/(2*h);
b(N) = b(N) - ue/(h^2) + ue/(2*h);

% Solve the system and append boundary conditions
u = A\b;
u = [us; u; ue];

% Calculate the grid 2-norm
err = sqrt(h*sum((u - Lab01_ExactSolution(x)).^2));