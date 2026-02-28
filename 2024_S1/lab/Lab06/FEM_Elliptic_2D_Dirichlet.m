function u = FEM_Elliptic_2D_Dirichlet(node, elem, DFunc, qFunc, fFunc, gFunc)
    N = size(node, 1);
    NT = size(elem, 1);
    A = sparse(N, N);
    RHS = zeros(N, 1);
    for j = 1:NT
        P1 = node(elem(j, 1), :);
        P2 = node(elem(j, 2), :);
        P3 = node(elem(j, 3), :);
        Kej = elem_stiff_1(P1, P2, P3, DFunc);
        Mej = elem_mass_1(P1, P2, P3, qFunc);
        Fej = elem_load_1(P1, P2, P3, fFunc);
        A(elem(j, :), elem(j, :)) = A(elem(j, :), elem(j, :)) + Kej + Mej;
        RHS(elem(j, :), 1) = RHS(elem(j, :), 1) + Fej;
    end
    
    a = min(node(:, 1));
    b = max(node(:, 1));
    c = min(node(:, 2));
    d = max(node(:, 2));
    eps = 1e-10;
    isLeftBnd = abs(node(:, 1) - a) < eps;
    isRightBnd = abs(node(:, 1) - b) < eps;
    isBottomBnd = abs(node(:, 2) - c) < eps;
    isTopBnd = abs(node(:, 2) - d) < eps;
    isBndNode = false(N, 1);
    isBndNode(isLeftBnd) = true;
    isBndNode(isRightBnd) = true;
    isBndNode(isBottomBnd) = true;
    isBndNode(isTopBnd) = true;
    bndNode = find(isBndNode);
    freeNode = find(~isBndNode);
    u = zeros(N, 1);
    u(bndNode) = gFunc(node(bndNode, :));

    RHS = RHS - A*u;
    u(freeNode) = A(freeNode, freeNode)\RHS(freeNode);
end