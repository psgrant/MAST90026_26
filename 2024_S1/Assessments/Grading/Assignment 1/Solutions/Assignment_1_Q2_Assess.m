clear
close all
clc

% Create a list of the number of Chebyshev Gauss Lobatto points to use
NSim = 100;
N = linspace(2, NSim+1, NSim);

% Define the ODE parameters
alpha = (3/8)*(cos(1) - 5*cos(3));
beta = -(1/4)*(1 + 5*cos(2))*sin(1);

ptfun = @(xfun)6*tan(xfun);
qtfun = @(xfun)-9*ones(size(xfun));
ftfun = @(xfun)21*sin(xfun);

err = zeros(1, NSim);

for i = 1:NSim
    
    % Create function inputs
    [D,x] = cheb(N(i));
    pt = ptfun(x);
    qt = qtfun(x);
    ft = ftfun(x);

    % Run the students collocation code for increasing N
    [x, u] = BvpSpec(D,pt,qt,ft,alpha,beta,N(i));

    % Approxmate the 2-norm error using Chebyshev Gauss Lobatto quadrature
    err(i) = sqrt((pi/N(i))*sum(sqrt(1 - x.^2).*(u - Assignment_1_Exact(x)).^2));
end

% Plot the error vs number of collocation points
loglog(N, err)
title('Error in the collocation scheme vs number of collocation points')
xlabel('N')
ylabel('||E||_2')


