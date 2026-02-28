clear
close all
clc

% Test functions
ptfun = @(xfun)6*tan(xfun);
qtfun = @(xfun)-9*ones(size(xfun));
ftfun = @(xfun)21*sin(xfun);

% Create a Chebyshev-Gauss-Lobatto grid to evaluate student function
Ntest = 50;
[~, x] = cheb(Ntest);

% Run students function
[D, q, f] = toSelfAdjointForm(x, ptfun(x), qtfun(x), ftfun(x));

% This calculates the mean value of (Numerical D/sec(x)^6)
studentconst = (1/2)*(pi/Ntest)*sum(sqrt(1 - x.^2).*(D./sec(x).^6));
% Analytical D accounting for students choice of constant
Dan = studentconst*sec(x).^6;
% Analytical q
qan = 9*Dan;
% Analyical f
fan = -21*sin(x).*Dan;

% Plot analytical and numerical results to test for correctly running code
subplot(1, 3, 1)
plot(x, D./Dan);
title('Numerical D/Analytical D')
subtitle('This plot should be approximately 1 everywhere')
xlabel('x')
subplot(1, 3, 2)
plot(x, q./qan)
title('Numerical q/Analytical q')
subtitle('This plot should be approximately 1 everywhere')
xlabel('x')
subplot(1, 3, 3)
plot(x, f./fan)
title('Numerical f/Analytical f')
subtitle('This plot should be approximately 1 everywhere')
xlabel('x')