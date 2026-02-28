% Spatial discretisation
xs = -1;
xe = 2;
Nx = 2000;

x = linspace(xs, xe, Nx + 1)';
xIn = x(2:Nx + 1);
h = (xe - xs)/Nx;

% Time discretisation
ts = 0;
te = 10;
Nt = 10000;
k = (te - ts)/Nt;

% Wave/Courrant numbers
c = 1;
nu = c*k/h;

% Function for initial condition
IC = @(xA)heaviside(xA + 0.5) - heaviside(xA - 0.5);

% Create A matrix
e = ones(Nx, 1);
A = spdiags([-e, e], [-1, 1], Nx, Nx);
A(1, Nx) = -1;
A(Nx, 1) = 1;
A = (nu/2)*A;

% Plot initial condition
u = IC(xIn);
plot(x, IC(x))
ylim([-0.1, 1.1])
pbaspect([1, 1, 1])
xlabel('x')
ylabel('u')
% Time step through the algorithm and only display plots at certain times
plotNum = ceil((Nt/(90*te*c)));
for i = 1:Nt
    if mod(i, plotNum) == 0
        plot(x, [u(Nx); u])
        ylim([-0.1, 1.1])
        xlim([xs, xe])
        pbaspect([1, 1, 1])
        xlabel('x')
        ylabel('u')
        pause(3*te*c/Nt)
    end
    % The actual algorithm
    u = u-A*u;
end