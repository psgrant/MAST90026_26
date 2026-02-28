% Set the domain of the PDE
xs = 0;
xe = 10;

% Number of internal points
N = 100;

% List of x values including end points
x = linspace(xs, xe, N + 2);
% Grid size
h = (xe - xs)/(N + 1);

% Set the initial condition and use this to define the boundary conditions
% Try all these initial conditions to see the kinds of behaviour we can expect
u0 = ones(size(x)) - (x > 1);
%u0 = ones(size(x)) + sin(x);
%u0 = -2*ones(size(x));

% The dirichlet boundary condition
us = 1;

% Extract the internal points
xIn = x(2:N + 2);
u0In = u0(2:N + 2);

% range of time values to solve for
tspan = [0 10];

% Create the RHS of the discretised system
e = ones(N + 1, 1);
A = spdiags([e/(h^2), -2*e/(h^2), e/(h^2)], -1:1, N+1, N+1);
A(N+1, N) = 2*A(N+1, N);

b = zeros(N+1, 1);
b(1) = us/(h^2);
RHS = @(u)A*u + u.*(1 - u) + b;

% Solve the discretised PDE using ODE45
options = odeset('RelTol', 1e-5);
[t, u] = ode45(@(t, u)RHS(u), tspan, u0In, options);

% Reformat results and plot
[X, T] = meshgrid(x, t);
u = [us*ones(size(t)) u];
figure(1)
mesh(X, T, u)
xlabel('x')
ylabel('t')
zlabel('u')