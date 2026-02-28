function u = FEM_Elliptic_2D_Dirichlet(node, elem, DFunc, qFunc, fFunc, bndNode, bndVal, freeNode)
    N = size(node, 1);
    NT = size(elem, 1);
    A = sparse(N, N);
    RHS = zeros(N, 1);
    for j = 1:NT
        P1 = node(elem(j, 1), :);
        P2 = node(elem(j, 2), :);
        P3 = node(elem(j, 3), :);
        Kej = elem_stiff(P1, P2, P3, DFunc);
        Mej = elem_mass(P1, P2, P3, qFunc);
        Fej = elem_load(P1, P2, P3, fFunc);
        A(elem(j, :), elem(j, :)) = A(elem(j, :), elem(j, :)) + Kej + Mej;
        RHS(elem(j, :), 1) = RHS(elem(j, :), 1) + Fej;
    end
    
    u = zeros(N, 1);
    u(bndNode) = bndVal;

    RHS = RHS - A*u;
    u(freeNode) = A(freeNode, freeNode)\RHS(freeNode);
end