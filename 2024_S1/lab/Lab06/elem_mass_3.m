function M = elem_mass_3(P1, P2, P3, q)
    JT = [P2(1) - P1(1), P3(1) - P1(1); P2(2) - P1(2), P3(2) - P1(2)];
    dJT = abs(det(JT));
    N1 = @(XI) 1 - XI(1) - XI(2);
    N2 = @(XI)XI(1);
    N3 = @(XI)XI(2);
    Ft = @(XI)JT*[XI(1); XI(2)] + [P1(1); P1(2)];
    int = @(XI)q(Ft(XI))*[N1(XI); N2(XI); N3(XI)]*[N1(XI), N2(XI), N3(XI)]*dJT;
    M = (int([0, 0]) + int([1, 0]) + int([0, 1]))/6;
end