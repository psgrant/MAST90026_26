function [JR, bR] = rectangle_map(P1, P2, P3, P4)
    % Initialize the Jacobian matrix JR and the translation vector bR
    JR = [P2(1) - P1(1), P4(1) - P1(1);
          P2(2) - P1(2), P4(2) - P1(2)];
    bR = P1(:); % Direct assignment using column vector format of P1
end
