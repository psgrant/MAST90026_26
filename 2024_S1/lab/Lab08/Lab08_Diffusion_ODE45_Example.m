% Set the domain of the PDE
xs = -2;
xe = 1;

% Number of internal points
N = 40;

% List of x values including end points
x = linspace(xs, xe, N + 2);
% Grid size
h = (xe - xs)/(N + 1);

% Set the initial condition and use this to define the boundary conditions
u0 = x.^2;

us = u0(1);
ue = u0(N + 2);

% Extract the internal points
xIn = x(2:N + 1);
u0In = u0(2:N + 1);

% range of time values to solve for
tspan = [0 5];

% Create the RHS of the discretised system
e = ones(N, 1);
A = spdiags([e/(h^2), -2*e/(h^2), e/(h^2)], -1:1, N, N);
b = zeros(N, 1);
b(1) = us/(h^2);
b(N) = ue/(h^2);

% Solve the discretised PDE using ODE45
options = odeset('RelTol', 1e-5);
[t, u] = ode45(@(t, u)A*u + b, tspan, u0In, options);

% Reformat results and plot
[X, T] = meshgrid(x, t);
u = [us*ones(size(t)) u ue*ones(size(t))];
figure(1)
mesh(X, T, u)