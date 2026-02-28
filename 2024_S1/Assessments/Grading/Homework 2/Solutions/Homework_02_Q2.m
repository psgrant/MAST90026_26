% Create a list of the number of Chebyshev Gauss Lobatto points to use
NSim = 50;
N = linspace(2, NSim+1, NSim);

% Define the ODE parameters
p = -3;
q = -1;
r = 5;
alpha = 1;
beta = -30;


err = zeros(1, NSim);

for i = 1:NSim
    % Set the number of collocation points
    NN = N(i);
    % Calculate the location of the collocation points and the differentiation matrix
    [D, x] = cheb(NN);

    % Construct the collocation matrix. Rows corresponding to boundary
    % conditions are the corresponding row of the identity matrix
    A = D^2 + p*D + q*eye(NN + 1);
    A(1, 1) = 1;
    A(1, 2:NN+1) = 0;
    A(NN+1, NN+1) = 1;
    A(NN+1, 1:NN) = 0;

    % Construct the RHS of the linear system
    b = r*ones(NN+1, 1);
    b(1) = beta;
    b(NN+1) = alpha;
    
    % Solve the collocation linear system
    u = A\b;

    % Approxmate the error using Chebyshev Gauss Lobatto quadrature
    err(i) = sqrt((pi/NN)*sum(sqrt(1 - x(2:end-1).^2).*(u(2:end-1) - Homework_02_Q2_Exact_Solution(x(2:end-1), p, q, r, alpha, beta)).^2));
end

% Plot the error vs number of collocation points
loglog(N, err)
title('Error in the collocation scheme vs number of collocation points')
xlabel('N')
ylabel('||E||_2')


