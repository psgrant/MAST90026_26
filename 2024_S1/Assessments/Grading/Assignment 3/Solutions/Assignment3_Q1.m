close all
clear
clc

% Set the domain of the PDE
xs = 0;
xe = 2;
ts = 0;
te = 5;

% Number of spacial points to solve over
N = 1000;

% List of x values including end points
x = linspace(xs, xe, N + 1)';
% Grid size
h = (xe - xs)/N;

% Number of time steps (satisfies CFL condition)
T1 = 5*N;
k1 = (te - ts)/T1;
nu1 = k1/h;

% Number of time steps (doesn't satisfy CFL condition)
T2 = 5*N/4;
k2 = (te - ts)/T2;
nu2 = k2/h;

% Set the two initial conditions
inita = (x >= 0.2) - (x >= 0.4);
initb = exp(-10*(4*x - 1).^2);

% u1 is initial condition a solved with Lax-Wendroff with CFL condition satisfied
uaLxWSt = inita(2:N+1);
% u2 is initial condition b solved with Lax-Wendroff with CFL condition satisfied
ubLxWSt = initb(2:N+1);
% u1 is initial condition a solved with Upwind with CFL condition satisfied
uaUWSt = inita(2:N+1);
% u1 is initial condition b solved with Upwind with CFL condition satisfied
ubUWSt = initb(2:N+1);
% u1 is initial condition a solved with Lax-Wendroff with CFL condition NOT satisfied
uaLxWUSt = inita(2:N+1);
% u1 is initial condition b solved with Lax-Wendroff with CFL condition NOT satisfied
ubLxWUSt = initb(2:N+1);
% u1 is initial condition a solved with Upwind with CFL condition NOT satisfied
uaUWUSt = inita(2:N+1);
% u1 is initial condition b solved with Upwind with CFL condition NOT satisfied
ubUWUSt = initb(2:N+1);

% Matrix advection term Lax-Wendroff
BLxW = spdiags([-1, 1], [-1, 1], N, N);
BLxW(1, N) = -1;
BLxW(N, 1) = 1;
% Matrix diffusion term Lax-Wendroff
ALxW = spdiags([1, -2, 1], -1:1, N, N);
ALxW(1, N) = 1;
ALxW(N, 1) = 1;
% Matrix advection term Upwind
BUW = spdiags([-1, 1], -1:0, N, N);
BUW(1, N) = -1;

% Solve for both initial conditions when the CFL condition is satisfied
for i = 1:T1
    uaLxWSt = uaLxWSt - nu1*BLxW*uaLxWSt/2 + nu1^2*ALxW*uaLxWSt/2;
    ubLxWSt = ubLxWSt - nu1*BLxW*ubLxWSt/2 + nu1^2*ALxW*ubLxWSt/2;
    uaUWSt = uaUWSt - nu1*BUW*uaUWSt;
    ubUWSt = ubUWSt - nu1*BUW*ubUWSt;
end

% Solve for both initial conditions when the CFL condition is NOT satisfied
for i = 1:T2
    uaLxWUSt = uaLxWUSt - nu2*BLxW*uaLxWUSt/2 + nu2^2*ALxW*uaLxWUSt/2;
    ubLxWUSt = ubLxWUSt - nu2*BLxW*ubLxWUSt/2 + nu2^2*ALxW*ubLxWUSt/2;
    uaUWUSt = uaUWUSt - nu2*BUW*uaUWUSt;
    ubUWUSt = ubUWUSt - nu2*BUW*ubUWUSt;
end

uaLxWSt = [uaLxWSt(end); uaLxWSt];
ubLxWSt = [ubLxWSt(end); ubLxWSt];
uaUWSt = [uaUWSt(end); uaUWSt];
ubUWSt = [ubUWSt(end); ubUWSt];
uaLxWUSt = [uaLxWUSt(end); uaLxWUSt];
ubLxWUSt = [ubLxWUSt(end); ubLxWUSt];
uaUWUSt = [uaUWUSt(end); uaUWUSt];
ubUWUSt = [ubUWUSt(end); ubUWUSt];

% Plot the results at t = 5. This is what your results should looks like (depending on how you've discretised)
figure(1)
subplot(2, 2, 1)
plot(x, uaLxWSt)
subplot(2, 2, 2)
plot(x, ubLxWSt)
subplot(2, 2, 3)
plot(x, uaUWSt)
subplot(2, 2, 4)
plot(x, ubUWSt)

figure(2)
subplot(2, 2, 1)
plot(x, uaLxWUSt)
subplot(2, 2, 2)
plot(x, ubLxWUSt)
subplot(2, 2, 3)
plot(x, uaUWUSt)
subplot(2, 2, 4)
plot(x, ubUWUSt)