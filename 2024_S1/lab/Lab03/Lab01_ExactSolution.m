function u = Lab01_ExactSolution(x)
% The exact solution to u''(x) - u'(x) - 6u(x) = 30x, u(0) = 1, u(1) = 2
    u = exp(-2*x).*((29+37*exp(4))*exp(3+5*x)-exp(2)*(37+29*exp(6))-5*(exp(10)-1)*exp(2*x).*(6*x-1))/(6*(exp(10)-1));