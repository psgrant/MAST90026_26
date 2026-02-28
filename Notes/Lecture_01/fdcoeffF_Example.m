% Range of x values
xS = 0;
xE = 1;

% Number of points
N = 5;

% Create vector of x evenly spaced values, h = (xE - xS)/(N - 1)
x = linspace(0, 1, N);

% 2nd order 
fdcoeffF(2, x(4), x(3:5))