N = 10;

% Start and end coordinates
xs = 0;
xe = 1;

% Boundary conditions
us = 1;
ue = 0;

% create vector of x values
x = linspace(xs, xe, N+2)';
x = x(2:N+1);

h = (xe-xs)/(N+1);

% Create indices and values of the matrix entries for use by sparse.
% R is row, C is column, diag is the diagonal entries, lowDiag is the lower
% diagonal entries and upDiag are the upper diagonal entries.
diagR = linspace(1, N, N)';
diagC = linspace(1, N, N)';
diag = -(2/(h^2))*ones(N, 1);

lowDiagR = linspace(2, N, N-1)';
lowDiagC = linspace(1, N-1, N-1)';
lowDiag = (1/(h^2))*ones(N-1, 1);

upDiagR = linspace(1, N-1, N-1)';
upDiagC = linspace(2, N, N-1)';
upDiag = (1/(h^2))*ones(N-1, 1);

% Collect all matrix entries into single vector
elemR = [diagR; lowDiagR; upDiagR];
elemC = [diagC; lowDiagC; upDiagC];
elem = [diag; lowDiag; upDiag];

% Create the FD matrix
A = sparse(elemR, elemC, elem);

% Use a linear function between the two boundary conditions as the initial guess
un = linspace(us, ue, N+2)';
un = un(2:N+1);

% Maximum number of iterations before declaring 'Doesn't Converge'
MAX_ITER = 100;
% Tolerance on the L2 norm before ending iteration
TOL = 1e-6;
for j = 1:MAX_ITER
    F = A*un + exp(un) - x;
    F(1) = F(1) + us/(h^2);
    F(N) = F(N) + ue/(h^2);

    J = A + spdiags(exp(un), 0, N, N);
    deltn = -J\F;
    un1 = un + deltn;
    err = sqrt(h*sum(deltn.^2));
    un = un1;
    if err < TOL
        disp(strcat('Returned solution is converged: ||E|| =  ', num2str(err)))
        break
    end
end

if j == MAX_ITER
    disp('Returned solution is NOT converged: ||E|| = ')
end

x = [xs; x; xe];
u = [us; un1; ue];

plot(x, u)