function u = Homework01_ExactSolution(x, eps)
% The exact solution to u''(x) - u'(x) - 6u(x) = 30x, u(0) = 1, u(1) = 2
    u = 1 + x + (exp(x/eps) - 1)/(exp(1/eps) - 1);