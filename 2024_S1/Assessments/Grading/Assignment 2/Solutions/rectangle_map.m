function [JR, bR] = rectangle_map(P1, P2, ~, P4)
    JR = [P2(1) - P1(1), P4(1) - P1(1); P2(2) - P1(2), P4(2) - P1(2)];
    bR = [P1(1); P1(2)];
end