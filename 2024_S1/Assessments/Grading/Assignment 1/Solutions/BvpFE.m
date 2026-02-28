function u = BvpFE(node, elem, D, q, f, alpha, beta)
    % Calculate the number of elements from the length of node
    N = length(node) - 1;
    h = node(2:end) - node(1:end-1);

    diagR = linspace(1, N+1, N+1)';
    diagC = linspace(1, N+1, N+1)';
    diag = zeros(N+1, 1);
    
    upDiagR = linspace(1, N, N)';
    upDiagC = linspace(2, N+1, N)';
    upDiag = zeros(N, 1);
    
    lowDiagR = linspace(2, N+1, N)';
    lowDiagC = linspace(1, N, N)';
    lowDiag = zeros(N, 1);
    
    RHS = zeros(N+1, 1);

    for l = 1:N
        KEl = elem_stiff(h(l), D(l), D(l+1));
        MEl = elem_mass(h(l), q(l), q(l+1));
        FEl = elem_load(h(l), f(l), f(l+1));
    
        RHS(elem(l, :), 1) = RHS(elem(l, :), 1) + FEl;

        if l == N
            RHS(N) = RHS(N) - beta*(KEl(2, 1) + MEl(2, 1));
        end
    
        diag(elem(l, 1)) = diag(elem(l, 1)) + KEl(1, 1) + MEl(1, 1);
        diag(elem(l, 2)) = diag(elem(l, 2)) + KEl(2, 2) + MEl(2, 2);
        upDiag(elem(l, 1)) = upDiag(elem(l, 1)) + KEl(1, 2) + MEl(1, 2);
        lowDiag(elem(l, 1)) = lowDiag(elem(l, 1)) + KEl(2, 1) + MEl(2, 1);
    end

    % Collect all matrix entries into single vector
    elemRow = [diagR; lowDiagR; upDiagR];
    elemCol = [diagC; lowDiagC; upDiagC];
    elemVal = [diag; lowDiag; upDiag];

    % Create the sparse FD matrix out of the first N rows and columns
    A = sparse(elemRow, elemCol, elemVal);
    A = A(1:N, 1:N);
    % Create the RHS from the first N entries and incorporate the Neumann BC
    RHS = RHS(1:N);
    RHS(1) = RHS(1) - D(1)*alpha;

    % Solve for u
    u = A\RHS;
    u = [u; beta];
end