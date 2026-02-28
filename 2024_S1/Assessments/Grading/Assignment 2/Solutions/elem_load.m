function F = elem_load(P1, P2, P3, f, quad_order)
    JT = [P2(1) - P1(1), P3(1) - P1(1); P2(2) - P1(2), P3(2) - P1(2)];
    dJT = abs(det(JT));
    N1 = @(XI) 1 - XI(1) - XI(2);
    N2 = @(XI)XI(1);
    N3 = @(XI)XI(2);
    Ft = @(XI)JT*[XI(1); XI(2)] + [P1(1); P1(2)];
    int = @(XI)f(Ft(XI))*[N1(XI); N2(XI); N3(XI)]*dJT;
    if quad_order == 1
        F = (int([0, 0]) + int([1, 0]) + int([0, 1]))/6;
    elseif quad_order == 2
        F = (int([0.5, 0]) + int([0.5, 0.5]) + int([0, 0.5]))/6;
    end
end