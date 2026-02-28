N = 20;

xs = -1;
xe = 1;

us = 1;
ue = 2;

x = linspace(xs, xe, N+1)';
D = exp(-x);
q = -6*exp(-x);
f = -30*x.*exp(-x);
h = x(2:N+1) - x(1:N);
elem = [linspace(1, N, N)', linspace(2, N+1, N)'];

diagR = linspace(1, N+1, N+1);
diagC = linspace(1, N+1, N+1);
diag = zeros(1, N+1);

upDiagR = linspace(1, N, N);
upDiagC = linspace(1, N, N);
upDiag = zeros(N, 1);

lowDiagR = linspace(1, N, N);
lowDiagC = linspace(1, N, N);
lowDiag = zeros(N, 1);

A = spalloc(N+1, N+1, 3*N+1);
b = zeros(N+1, 1);
for l = 1:N
    hl = h(l);
    Dl = D(l);
    Dl1 = D(l+1);
    ql = q(l);
    ql1 = q(l+1);
    fl = f(l);
    fl1 = f(l+1);
    KEl = elem_stiff(hl, Dl, Dl1);
    MEl = elem_mass(hl, ql, ql1);
    FEl = elem_load(hl, fl, fl1);
    A(elem(l, :), elem(l, :)) = A(elem(l, :), elem(l, :)) + KEl;
    A(elem(l, :), elem(l, :)) = A(elem(l, :), elem(l, :)) + MEl;
    b(elem(l, :), 1) = b(elem(l, :), 1) + FEl;
end

A(1, 1) = 1;
A(1, 2) = 0;
A(N+1, N) = 0;
A(N+1, N+1) = 1;

b(1) = us;
b(N+1) = ue;

u = A\b;
err = sqrt(sum(h.*((u(1:N) - Lab01_ExactSolution(x(1:N))).^2 + (u(2:N+1) - Lab01_ExactSolution(x(2:N+1))).^2)/2));

xEx = linspace(xs, xe, 1000);
plot(x, u, 'o', xEx, Lab01_ExactSolution(xEx))