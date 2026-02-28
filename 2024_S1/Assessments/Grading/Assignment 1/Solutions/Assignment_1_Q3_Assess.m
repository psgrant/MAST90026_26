clear
close all
clc

% Number of simulations to run. h is halved each time
NStart = 3;
NSim = 10;
N = 2.^linspace(NStart, NStart + NSim, NSim + 1);

Dfun = @(xfun)sec(xfun).^6;
qfun = @(xfun)9*sec(xfun).^6;
ffun = @(xfun)-21*sec(xfun).^5.*tan(xfun);

a = -3/4;
b = 1/2;

% Fix BCs
alpha = (3/8)*(cos(3/4) - 5*cos(9/4));
beta = -(1/4)*(1 + 5*cos(1))*sin(1/2);

hmax = zeros(1, NSim);
err = zeros(1, NSim);
t = zeros(1, NSim);
% run the simulations and only save the error
for i = 1:NSim
    tic
    [~, node] = cheb(N(i));
    node = node(end:-1:1);
    node = ((b - a)*node + a + b)/2;
    hmax(i) = max(node(2:end) - node(1:end-1));
    elem = [(1:N(i))', (2:N(i)+1)'];
    u = BvpFE(node, elem, Dfun(node), qfun(node), ffun(node), alpha, beta);
    err(i) = sqrt(trapz(node, (u - Assignment_1_Exact(node)).^2));
    t(i) = toc;
end

% loglog plot of the error.
subplot(1, 2, 1)
loglog(hmax, err)
title('Error vs grid size')
xlabel('h')
ylabel('||E||_2')
subplot(1, 2, 2)
loglog(t, err)
title('Error vs execution time')
xlabel('t')
ylabel('||E||_2')