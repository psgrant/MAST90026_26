function K = elem_stiff_3(P1, P2, P3, D)
    JT = [P2(1) - P1(1), P3(1) - P1(1); P2(2) - P1(2), P3(2) - P1(2)];
    Ft = @(XI)JT*[XI(1); XI(2)] + [P1(1); P1(2)];
    dJT = det(JT);
    iJTT = [JT(2, 2), -JT(2, 1); -JT(1, 2), JT(1, 1)]/dJT;
    dJT = abs(dJT);
    gradN1 = [-1; -1];
    gradN2 = [1; 0];
    gradN3 = [0; 1];
    int11 = (iJTT*gradN1)'*(iJTT*gradN1);
    int12 = (iJTT*gradN1)'*(iJTT*gradN2);
    int13 = (iJTT*gradN1)'*(iJTT*gradN3);
    int22 = (iJTT*gradN2)'*(iJTT*gradN2);
    int23 = (iJTT*gradN2)'*(iJTT*gradN3);
    int33 = (iJTT*gradN3)'*(iJTT*gradN3);
    int = @(XI)D(Ft(XI))*[int11, int12, int13; int12, int22, int23; int13, int23, int33]*dJT;
    K = (int([0, 0]) + int([1, 0]) + int([0, 1]))/6;
end