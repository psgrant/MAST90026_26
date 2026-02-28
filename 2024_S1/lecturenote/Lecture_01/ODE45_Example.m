% Parameters
PI = 0;
beta = 0.0095;
delt = 0.0001;
zeta = 0.1;
alph = 0.005;

% Initial Conditions
S0 = 1000;
Z0 = 0;
R0 = 0;
X0 = [S0; Z0; R0];

% ODE
dXdt = @(t, X) [PI - beta*X(1)*X(2) - delt*X(1); ...
    beta*X(1)*X(2) + zeta*X(3) - alph*X(1)*X(2); ...
    delt*X(1) + alph*X(1)*X(2) - zeta*X(3)];

% Length of simulation
tEnd = 50;

% Solve ODE
[t, X] = ode45(dXdt, [0, tEnd], X0);

% Plot results
plot(t, X)
