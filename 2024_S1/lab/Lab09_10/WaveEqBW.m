% Spatial discretisation
xs = -1;
xe = 2;
Nx = 10000;

x = linspace(xs, xe, Nx + 1)';
xIn = x(2:Nx + 1);
h = (xe - xs)/Nx;

% Time discretisation
ts = 0;
te = 10;
Nt = 100000;
k = (te - ts)/Nt;

% Wave/Courrant numbers
c = 1;
nu = c*k/h;
disp(nu)

% Function for initial condition

IC = @(xA)heaviside(xA + 0.5) - heaviside(xA - 0.5);

e = ones(Nx, 1);
if c > 0
    Aadv = spdiags([e, -4*e, 3*e], -2:0, Nx, Nx);
    Aadv(1, Nx) = -4;
    Aadv(1, Nx - 1) = 1;
    Aadv(2, Nx) = 1;
    Aadv = (nu/2)*Aadv;
    Adiff = spdiags([e, -2*e, e], -2:0, Nx, Nx);
    Adiff(1, Nx) = -2;
    Adiff(1, Nx - 1) = 1;
    Adiff(2, Nx) = 1;
    Adiff = (nu^2/2)*Adiff;
elseif c < 0
    Aadv = spdiags([-3*e, 4*e, -e], 0:2, Nx, Nx);
    Aadv(Nx, 1) = 4;
    Aadv(Nx, 2) = -1;
    Aadv(Nx - 1, 1) = -1;
    Aadv = nu*Aadv;
    Adiff = spdiags([e, -2*e, e], 0:2, Nx, Nx);
    Adiff(Nx, 1) = -2;
    Adiff(Nx, 2) = 1;
    Adiff(Nx - 1, 1) = 1;
    Adiff = (nu^2/2)*Adiff;
end
B = speye(Nx) - Aadv + Adiff;

u = IC(xIn);
plot(x, IC(x))
ylim([-0.1, 1.1])
pbaspect([1, 1, 1])
xlabel('x')
ylabel('u')
plotNum = ceil((1/30)/(3*te*abs(c)/(Nt)));
for i = 1:Nt
    if mod(i, plotNum) == 0
        plot(x, [u(Nx); u])
        ylim([-0.1, 1.1])
        xlim([xs, xe])
        pbaspect([1, 1, 1])
        xlabel('x')
        ylabel('u')
        pause(3*te*abs(c)/Nt)
    end
    u = B*u;
end