% Spatial discretisation
xs = -1;
xe = 2;
Nx = 200;

x = linspace(xs, xe, Nx + 1)';
xIn = x(2:Nx + 1);
h = (xe - xs)/Nx;

% Time discretisation
ts = 0;
te = 5;
Nt = 1000;
k = (te - ts)/Nt;

% Wave/Courrant numbers
c = 1;
nu = c*k/h;
disp(nu)

% Function for initial condition

IC = @(xA)heaviside(xA + 0.5) - heaviside(xA - 0.5);

e = ones(Nx, 1);
if c > 0
    Askew = spdiags([-e, e], [-1, 0], Nx, Nx);
    Askew(1, Nx) = -1;
    Askew = nu*Askew;
elseif c < 0
    Askew = spdiags([-e, e], [0, 1], Nx, Nx);
    Askew(Nx, 1) = 1;
    Askew = nu*Askew;
end

u = IC(xIn);
plot(x, IC(x))
ylim([-0.1, 1.1])
pbaspect([1, 1, 1])
xlabel('x')
ylabel('u')
plotNum = ceil((Nt/(90*te*c)));
for i = 1:Nt
    if mod(i, plotNum) == 0
        plot(x, [u(Nx); u])
        ylim([-0.1, 1.1])
        pbaspect([1, 1, 1])
        pause(3*te*c/Nt)
    end

    u = u - Askew*u;
end