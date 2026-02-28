N = 10;

% Start and end coordinates
xs = 0;
xe = 1;

% Boundary conditions
us = 0;
ue = 1;

% create vector of x values
x = linspace(xs, xe, N+2)';
x = x(2:N+1);

h = (xe-xs)/(N+1);

% Create the FD matrix
e = ones(N, 1);
A = spdiags([e, -2*e, e]/h^2, -1:1, N, N);

% Use a linear function between the two boundary conditions as the initial guess
un = linspace(us, ue, N+2)';
un = un(2:N+1);

% Maximum number of iterations before declaring 'Doesn't Converge'
MAX_ITER = 100;
% Tolerance on the L2 norm before ending iteration
TOL = 1e-6;
for j = 1:MAX_ITER
    % Create the non-linear equation to be satisfied
    F = A*un - sin(un) + ones(size(un));
    F(1) = F(1) + us/(h^2);
    F(N) = F(N) + ue/(h^2);

    % Create the Jacobian matrix
    J = A - spdiags(cos(un), 0, N, N);
    % Solve the linearised system and perform the Newton update
    deltn = -J\F;
    un1 = un + deltn;
    % Calculate the error on the residual
    err = sqrt(h*sum(deltn.^2));
    un = un1;
    if err < TOL
        % If solution has converged end the for loop
        disp(strcat('Returned solution is converged: ||E|| =  ', num2str(err)))
        break
    end
end

% If for loop ran to max iteration give warning solution is not converged
if j == MAX_ITER
    disp('Returned solution is NOT converged: ||E|| = ')
end

% Create final solution vector including BCs and then plot
x = [xs; x; xe];
u = [us; un1; ue];

plot(x, u)
xlabel('x')
ylabel('u')
title('')