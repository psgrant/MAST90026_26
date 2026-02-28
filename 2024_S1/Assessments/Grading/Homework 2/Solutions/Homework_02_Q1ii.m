function [h, err, x, u] = Homework_02_Q1ii(N, p, q, r, a, b, alpha, beta)

x = linspace(a, b, N+2)';

h = (b-a)/(N+1);

% Lower diagonal coefficient
lC = 1/(h^2) - p/(2*h);
% Diagonal coefficient
dC = -2/(h^2) + q;
% Upper diagonal coefficient
uC = p/(2*h) + 1/(h^2);

% Create the FD matrix
e = ones(N, 1);
A = spdiags([[lC*e; 0], [-2*(1+h)/(h^2) + p + q; dC*e], [0; 2/h^2; uC*e(1:N-1)]], -1:1, N+1, N+1);

% Create RHS vector
RHS = [r + alpha*(p - 2/h); r*ones(N - 1, 1); r - beta*(1/(h^2) + p/(2*h))];

% Solve the system and append boundary conditions
u = A\RHS;

u = [u; beta];

% Calculate the grid 2-norm
err = sqrt(h*sum((u - Homework_02_Exact_Solution(x, p, q, r, a, b, alpha, beta)).^2));

end

