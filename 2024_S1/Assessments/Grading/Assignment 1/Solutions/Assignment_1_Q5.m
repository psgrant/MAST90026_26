% Number of points
N = 40;

% BC locations
a = -1;
b = 1;

% BCs
alpha = 2;
beta = 0;

node = linspace(a, b, N+1)';
elem = [(1:N)', (2:N+1)'];
% Convert to self adjoint form
[pt, qt, ft] = Assignment_1_Q5_BVP_Functions(node);
[D, q, f] = toSelfAdjointForm(node, pt, qt, ft);
u = BvpFE(node, elem, D, q, f, alpha, beta);

plot(node, u)