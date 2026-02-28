% Spatial discretisation
xs = -1;
xe = 2;
Nx = 100;

x = linspace(xs, xe, Nx + 1);
xIn = x(2:Nx + 1);
h = (xe - xs)/Nx;

% Time discretisation
ts = 0;
te = 10;
tspan = [ts te];

% Wave number
c = 1;

% Function for initial condition

IC = @(xA)heaviside(xA + 0.5) - heaviside(xA - 0.5);

e = ones(Nx, 1);
A = spdiags([c*e/(2*h), -c*e/(2*h)], [-1, 1], Nx, Nx);
A(1, Nx) = c/(2*h);
A(Nx, 1) = -c/(2*h);

% Solve the discretised PDE using ODE45
options = odeset('RelTol', 1e-5);
[t, u] = ode45(@(t, u)A*u, tspan, IC(xIn), options);

u = [u(:, Nx) u];

Nt = length(t);

plot(x, IC(x))
ylim([-0.1, 1.1])
pbaspect([1, 1, 1])
xlabel('x')
ylabel('u')
plotNum = ceil((Nt/(90*te*c)));
for i = 1:Nt
    if mod(i, plotNum) == 0
        plot(x, u(i, :))
        ylim([-0.1, 1.1])
        pbaspect([1, 1, 1])
        xlabel('x')
        ylabel('u')
        pause(3*te*c/Nt)
    end
end