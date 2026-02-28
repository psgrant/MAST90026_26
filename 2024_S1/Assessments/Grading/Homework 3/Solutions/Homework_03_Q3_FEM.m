function [h, err, x, u] =  Homework_03_Q3_FEM(N, q, r, a, b, alpha, beta)

x = linspace(a, b, N+1)';
Dx = ones(N+1, 1);
qx = q*ones(N+1, 1);
fx = r*ones(N+1, 1);
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
RHS = zeros(N+1, 1);
for l = 1:N
    hl = h(l);
    Dl = Dx(l);
    Dl1 = Dx(l+1);
    ql = qx(l);
    ql1 = qx(l+1);
    fl = fx(l);
    fl1 = fx(l+1);
    KEl = elem_stiff(hl, Dl, Dl1);
    MEl = elem_mass(hl, ql, ql1);
    FEl = elem_load(hl, fl, fl1);

    RHS(elem(l, :), 1) = RHS(elem(l, :), 1) + FEl;
    if l == 1
        RHS(2) = RHS(2) - alpha*(KEl(1, 2) + MEl(1, 2));
    elseif l == N
        RHS(N) = RHS(N) - beta*(KEl(2, 1) + MEl(2, 1));
    end

    diag(l) = diag(l) + KEl(1, 1) + MEl(1, 1);
    diag(l+1) = diag(l+1) + KEl(2, 2) + MEl(2, 2);
    upDiag(l) = upDiag(l) + KEl(1, 2) + MEl(1, 2);
    lowDiag(l) = lowDiag(l) + KEl(2, 1) + MEl(2, 1);
end

% Collect all matrix entries into single vector
elemRow = [diagR; lowDiagR; upDiagR];
elemCol = [diagC; lowDiagC; upDiagC];
elemVal = [diag; lowDiag; upDiag];

% Create the FD matrix
A = sparse(elemRow, elemCol, elemVal);

A = A(2:N, 2:N);
RHS = RHS(2:N);

u = A\RHS;
u = [alpha; u; beta];
err = sqrt(trapz(x, (u - Homework_03_Q3_Exact_Solution(x, q, r, a, b, alpha, beta)).^2));
h = max(h);

end