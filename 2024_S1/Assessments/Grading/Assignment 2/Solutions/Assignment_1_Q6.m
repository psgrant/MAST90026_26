% Number of internal points
N = 100;

% Value of the parameter
L = 50;

% Start and end coordinates
xs = 0;
xe = 1;

% Boundary conditions
us = 0;
ue = 0;

% create vector of x values
x = linspace(xs, xe, N+2)';
x = x(2:N+1);

h = (xe-xs)/(N+1);

% Create the FD matrix
e = ones(N, 1);
A = spdiags([e, -2*e, e]/h^2, -1:1, N, N);

% Use a constant initial guess
un = -2*ones(N+2, 1);
un = un(2:N+1);

% Maximum number of iterations before declaring 'Doesn't Converge'
MAX_ITER = 100;
% Tolerance on the L2 norm before ending iteration
TOL = 1e-6;
for j = 1:MAX_ITER
    F = -A*un + L^2*(un.*(un - ones(size(un))));
    F(1) = F(1) + us/(h^2);
    F(N) = F(N) + ue/(h^2);

    J = -A + L^2*spdiags(2*un - ones(size(un)), 0, N, N);
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
    disp(strcat('Returned solution is NOT converged: ||E|| = ', num2str(err)))
end

x = [xs; x; xe];
u = [us; un1; ue];

plot(x, u)
xlabel('$$\tilde{x}$$', 'Interpreter', 'LaTeX')
ylabel('u')