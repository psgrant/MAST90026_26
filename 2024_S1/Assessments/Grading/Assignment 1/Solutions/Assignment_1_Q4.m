clear
close all
clc

% Number of simulations to run. h is halved each time
NStart = 1;
NSim = 10;
N = 2.^linspace(NStart, NStart + NSim, NSim + 1);

% Value of p to analyse
p = 1;

% Coefficient functions
ptfun = @(xfun)zeros(size(xfun));
qtfun = @(xfun)-32*xfun.^6./(p+xfun.^4).^2;
ftfun = @(xfun)-12*xfun.^2./(p+xfun.^4).^2;

% Analytic solution
uAn = @(xfun)1./(p + xfun.^4);

% BC locations
a = -1;
b = 1;

% BCs
alpha = 4/(1+p)^2;
beta = 1/(1+p);

% Run FEM and calculate 2-norm of error for evenly spaced points
h = zeros(1, NSim);
err1 = zeros(1, NSim);
for i = 1:NSim
    node = linspace(a, b, N(i)+1)';
    h(i) = max(node(2:end) - node(1:end-1));
    elem = [(1:N(i))', (2:N(i)+1)'];
    % Convert to self adjoint form
    [D, q, f] = toSelfAdjointForm(node, ptfun(node), qtfun(node), ftfun(node));
    u = BvpFE(node, elem, D, q, f, alpha, beta);
    err1(i) = sqrt(trapz(node, (u - uAn(node)).^2));
end

% Run FEM and calculate 2-norm of error for CGL points
hmax = zeros(1, NSim);
err2 = zeros(1, NSim);
for i = 1:NSim
    [~, node] = cheb(N(i));
    node = node(end:-1:1);
    node = ((b - a)*node + a + b)/2;
    hmax(i) = max(node(2:end) - node(1:end-1));
    elem = [(1:N(i))', (2:N(i)+1)'];
    % Convert to self adjoint form
    [D, q, f] = toSelfAdjointForm(node, ptfun(node), qtfun(node), ftfun(node));
    u = BvpFE(node, elem, D, q, f, alpha, beta);
    err2(i) = sqrt(trapz(node, (u - uAn(node)).^2));
end

figure(1)
loglog(h, err1, hmax, err2)
legend('FEM evenly spaced', 'FEM CGL points')
xlabel('h_{max}')
ylabel('||E||_2')


% Create a list of the number of Chebyshev Gauss Lobatto points to use
NSim = 100;
N = linspace(2, NSim+1, NSim);

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
    err(i) = sqrt((pi/N(i))*sum(sqrt(1 - x.^2).*(u - uAn(x)).^2));
end

% Plot the error vs number of collocation points
figure(2)
loglog(N, err)
title('Error in the collocation scheme vs number of collocation points')
xlabel('N')
ylabel('||E||_2')