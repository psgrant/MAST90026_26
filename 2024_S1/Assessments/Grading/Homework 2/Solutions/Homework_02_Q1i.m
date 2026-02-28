function [h, err, x, u] = Homework_02_Q1i(N, p, q, r, a, b, alpha, beta)

x = linspace(a, b, N+2)';

h = (b-a)/(N+1);

% Lower diagonal coefficient
lC = 1/(h^2) - p/(2*h);
% Diagonal coefficient
dC = -2/(h^2) + q;
% Upper diagonal coefficient
uC = 1/(h^2) + p/(2*h);

% Create the FD matrix
e = ones(N, 1);
A = spdiags([[lC*e; 0; 0], [1 + 1/h; dC*e; 1], [0; -1/h; uC*e]], -1:1, N+2, N+2);

% Create RHS vector
RHS = [alpha; r*ones(N, 1); beta];

% Solve the system and append boundary conditions
u = A\RHS;

% Calculate the grid 2-norm
err = sqrt(h*sum((u - Homework_02_Exact_Solution(x, p, q, r, a, b, alpha, beta)).^2));

end

