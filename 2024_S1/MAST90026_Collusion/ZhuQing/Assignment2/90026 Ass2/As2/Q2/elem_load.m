function F = elem_load(P1, P2, P3, f, quad_order)
    j11 = P2(1) - P1(1);
    j12 = P3(1) - P1(1);
    j21 = P2(2) - P1(2);
    j22 = P3(2) - P1(2);
    K = abs(j11 * j22 - j12 * j21);
    Load = zeros(3, 1);
    if (quad_order == 1)
        Load(1) = f(P1(1),P1(2))*K/6;
        Load(2) = f(P2(1),P2(2))*K/6;
        Load(3) = f(P3(1),P3(2))*K/6;
        F = Load;
    else 
        Load(1) = (f((P1(1) + P2(1))/2,(P1(2) + P2(2))/2)*K/2+f((P1(1) + P3(1))/2,(P1(2) + P3(2))/2)*K/2)/6;   
        Load(2) = (f((P1(1) + P2(1))/2,(P1(2) + P2(2))/2)*K/2+f((P2(1) + P3(1))/2,(P2(2) + P3(2))/2)*K/2)/6;
        Load(3) = (f((P1(1) + P3(1))/2,(P1(2) + P3(2))/2)*K/2+f((P2(1) + P3(1))/2,(P2(2) + P3(2))/2)*K/2)/6;
        F = Load;
    end
end
