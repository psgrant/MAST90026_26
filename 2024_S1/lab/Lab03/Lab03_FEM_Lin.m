N = 10;

xs = -1;
xe = 1;

us = 1;
ue = 2;

x = linspace(xs, xe, N+1)';
D = exp(-x);
q = 6*exp(-x);
f = -30*x.*exp(-x);
h = x(2:N+1) - x(1:N);
elem = [linspace(1, N, N)', linspace(2, N+1, N)'];

diagR = linspace(1, N+1, N+1)';
diagC = linspace(1, N+1, N+1)';
diag = zeros(N+1, 1);

upDiagR = linspace(1, N, N)';
upDiagC = linspace(2, N+1, N)';
upDiag = zeros(N, 1);

lowDiagR = linspace(2, N+1, N)';
lowDiagC = linspace(1, N, N)';
lowDiag = zeros(N, 1);

%A = spalloc(N+1, N+1, 3*N+1);
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

    b(elem(l, :), 1) = b(elem(l, :), 1) + FEl;

    diag(l) = diag(l) + KEl(1, 1) + MEl(1, 1);
    diag(l+1) = diag(l+1) + KEl(2, 2) + MEl(2, 2);
    upDiag(l) = upDiag(l) + KEl(1, 2);
    lowDiag(l) = lowDiag(l) + KEl(2, 1);
end

% Collect all matrix entries into single vector
elemRow = [diagR; lowDiagR; upDiagR];
elemCol = [diagC; lowDiagC; upDiagC];
elemVal = [diag; lowDiag; upDiag];

% Create the FD matrix
A = sparse(elemRow, elemCol, elemVal);

A = A(2:N, 2:N);
b = b(2:N);

KE1 = elem_stiff(h(1), D(1), D(2));
b(1) = b(1) - us*KE1(1, 2);
KEN = elem_stiff(h(N), D(N), D(N+1));
b(N-1) = b(N-1) - ue*KEN(1, 2);

u = A\b;
u = [us; u; ue];
err = sqrt(sum(h.*((u(1:N) - Lab01_ExactSolution(x(1:N))).^2 + (u(2:N+1) - Lab01_ExactSolution(x(2:N+1))).^2)/2));


xEx = linspace(xs, xe, 1000);
plot(x, u, 'o', xEx, Lab01_ExactSolution(xEx))