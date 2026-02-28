% Spatial discretisation
xs = -1;
xe = 2;
Nx = 100;

x = linspace(xs, xe, Nx + 1);
h = (xe - xs)/Nx;

% Time discretisation
ts = 0;
te = 10;

Nt = 3000;
t = linspace(ts, te, Nt + 1);
k = (te - ts)/Nt;

% Wave number
c = 1;

% Function for initial condition

IC = @(xA)heaviside(xA + 0.5) - heaviside(xA - 0.5);

plot(x, IC(x))
ylim([-0.1, 1.1])
pbaspect([1, 1, 1])
xlabel('x')
ylabel('u')
plotNum = ceil((Nt/(90*te*c)));
for i = 1:Nt
    if mod(i, plotNum) == 0
        plot(x, IC(mod(x - c*t(i) - xs, xe - xs) + xs))
        ylim([-0.1, 1.1])
        pbaspect([1, 1, 1])
        xlabel('x')
        ylabel('u')
        pause(3*te*c/Nt)
    end
end