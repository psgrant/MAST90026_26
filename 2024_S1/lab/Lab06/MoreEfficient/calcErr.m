function err = calcErr(node, elem, u, uExact)

    NT = size(elem, 1);
    err = 0;
    for j = 1:NT
        P1 = node(elem(j, 1), :);
        P2 = node(elem(j, 2), :);
        P3 = node(elem(j, 3), :);
        JT = [P2(1) - P1(1), P3(1) - P1(1); P2(2) - P1(2), P3(2) - P1(2)];
        dJT = abs(det(JT));
        err = err + ((u(elem(j, 1)) - uExact(P1)).^2 + (u(elem(j, 2)) - uExact(P2)).^2 + (u(elem(j, 3)) - uExact(P3)).^2)*dJT/6;
    end
    err = sqrt(err);
end