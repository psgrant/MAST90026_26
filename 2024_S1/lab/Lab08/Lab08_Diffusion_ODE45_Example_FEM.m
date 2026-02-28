% Set the domain of the PDE
xs = -2;
xe = 1;

% Number of elements
N = 50;

% range of time values to solve for
tspan = [0 5];

% List of x values including end points
x = linspace(xs, xe, N+1)';

% Set the initial condition and use this to define the boundary conditions
U0 = x.^2;

us = U0(1);
ue = U0(N + 1);

% Extract the internal points
xIn = x(2:N);
U0In = U0(2:N);

D = ones(N+1, 1)';
q = ones(N+1, 1)';
f = zeros(N+1, 1)';
h = x(2:N+1) - x(1:N);
elem = [linspace(1, N, N)', linspace(2, N+1, N)'];

diagRA = linspace(1, N+1, N+1)';
diagCA = linspace(1, N+1, N+1)';
diagA = zeros(N+1, 1);

upDiagRA = linspace(1, N, N)';
upDiagCA = linspace(2, N+1, N)';
upDiagA = zeros(N, 1);

lowDiagRA = linspace(2, N+1, N)';
lowDiagCA = linspace(1, N, N)';
lowDiagA = zeros(N, 1);

diagRM = linspace(1, N+1, N+1)';
diagCM = linspace(1, N+1, N+1)';
diagM = zeros(N+1, 1);

%A = spalloc(N+1, N+1, 3*N+1);
b = zeros(N+1, 1);
bOld = zeros(N+1, 1);
for l = 1:N
    hl = h(l);
    Dl = D(l);
    Dl1 = D(l+1);
    ql = q(l);
    ql1 = q(l+1);
    fl = f(l);
    fl1 = f(l+1);
    U0l = U0(l);
    U0l1 = U0(l+1);
    KEl = elem_stiff(hl, Dl, Dl1);
    MEl = elem_mass(hl, ql, ql1);
    FEl = elem_load(hl, fl, fl1);
    FICEl = elem_IC(hl, Dl, Dl1, U0l, U0l1);

    b(elem(l, :), 1) = b(elem(l, :), 1) + FEl + FICEl;

    diagA(l) = diagA(l) + KEl(1, 1);
    diagA(l+1) = diagA(l+1) + KEl(2, 2);
    upDiagA(l) = upDiagA(l) + KEl(1, 2);
    lowDiagA(l) = lowDiagA(l) + KEl(2, 1);

    diagM(l) = diagM(l) + MEl(1, 1);
    diagM(l+1) = diagM(l+1) + MEl(2, 2);
end

% Collect all matrix entries into single vector
elemRowA = [diagRA; lowDiagRA; upDiagRA];
elemColA = [diagCA; lowDiagCA; upDiagCA];
elemValA = [diagA; lowDiagA; upDiagA];

% Create the FD matrix
A = sparse(elemRowA, elemColA, elemValA);

% Create the FD matrix
M = sparse(diagRM, diagCM, diagM);

A = A(2:N, 2:N);
M = M(2:N, 2:N);
b = b(2:N);

% Solve the discretised PDE using ODE45
options = odeset('RelTol', 1e-5);
[t, u] = ode45(@(t, u)M\(-A*u + b), tspan, zeros(N-1, 1), options);

% Reformat results and plot
[X, T] = meshgrid(x, t);
u = [us*ones(size(t)), u + (U0In*ones(1, length(t)))', ue*ones(size(t))];
figure(2)
mesh(X, T, u)