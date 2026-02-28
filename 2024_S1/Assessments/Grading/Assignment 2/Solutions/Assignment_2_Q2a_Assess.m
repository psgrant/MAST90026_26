clear
close all
clc

% Test element
x1 = 1;
y1 = 2.5;
x2 = 2;
y2 = 2;
x3 = 2;
y3 = 3;

P1 = [x1, y1];
P2 = [x2, y2];
P3 = [x3, y3];

% Test functions
D = @(X)sin(X(1).*X(2));
q = @(X)exp(sqrt(X(1) + X(2)));
f = @(X)X(1) + X(2);

% Calculate using linearly accurate quadrature
Fej1 = elem_load(P1, P2, P3, f, 1);
Mej1 = elem_mass(P1, P2, P3, q, 1);
Kej1 = elem_stiff(P1, P2, P3, D, 1);

% Calculate using quadratically accurate quadrature
Fej2 = elem_load(P1, P2, P3, f, 2);
Mej2 = elem_mass(P1, P2, P3, q, 2);
Kej2 = elem_stiff(P1, P2, P3, D, 2);

% Compare to correct solution
disp('Linearly Accurate Results')
disp(' ')
disp(['Student:  Fej = [', num2str(Fej1(1), '%.4f'), '; ', num2str(Fej1(2), '%.4f'), '; ', num2str(Fej1(3), '%.4f'), ']'])
disp(['Solution: Fej = [', num2str(0.5833), '; ', num2str(0.6667), '; ', num2str(0.8333), ']'])
disp(' ')
disp(['                ⌈', num2str(Mej1(1, 1), '%.4f'), ' ', num2str(Mej1(1, 2), '%.4f'), ' ', num2str(Mej1(1, 3), '%.4f'), '⌉'])
disp(['Student:  Mej = |', num2str(Mej1(2, 1), '%.4f'), ' ', num2str(Mej1(2, 2), '%.4f'), ' ', num2str(Mej1(2, 3), '%.4f'), '|'])
disp(['                ⌊', num2str(Mej1(3, 1), '%.4f'), ' ', num2str(Mej1(3, 2), '%.4f'), ' ', num2str(Mej1(3, 3), '%.4f'), '⌋'])
disp(['                ⌈', num2str(1.0823), ' ', num2str(0, '%.4f'), ' ', num2str(0, '%.4f'), '⌉'])
disp(['Solution: Mej = |', num2str(0, '%.4f'), ' ', num2str(1.2315), ' ', num2str(0, '%.4f'), '|'])
disp(['                ⌊', num2str(0, '%.4f'), ' ', num2str(0, '%.4f'), ' ', num2str(1.5594), '⌋'])
disp(' ')
disp(['                ⌈', num2str(Kej1(1, 1), '%.4f'), '  ', num2str(Kej1(1, 2), '%.4f'), '  ', num2str(Kej1(1, 3), '%.4f'), '⌉'])
disp(['Student:  Kej = | ', num2str(Kej1(2, 1), '%.4f'), ' ', num2str(Kej1(2, 2), '%.4f'), '  ', num2str(Kej1(2, 3), '%.4f'), '|'])
disp(['                ⌊ ', num2str(Kej1(3, 1), '%.4f'), '  ', num2str(Kej1(3, 2), '%.4f'), ' ', num2str(Kej1(3, 3), '%.4f'), '⌋'])
disp(['                ⌈', num2str(-0.0730, '%.4f'), '  ', num2str(0.0365), '  ', num2str(0.0365), '⌉'])
disp(['Solution: Kej = | ', num2str(0.0365), ' ', num2str(-0.0912), '  ', num2str(0.0547), '|'])
disp(['                ⌊ ', num2str(0.0365), '  ', num2str(0.0547), ' ', num2str(-0.0912), '⌋'])
disp(' ')
disp(' ')
disp('Quadratically Accurate Results')
disp(' ')
disp(['Student:  Fej = [', num2str(Fej2(1), '%.4f'), '; ', num2str(Fej2(2), '%.4f'), '; ', num2str(Fej2(3), '%.4f'), ']'])
disp(['Solution: Fej = [', num2str(0.6667), '; ', num2str(0.6875), '; ', num2str(0.7292), ']'])
disp(' ')
disp(['                ⌈', num2str(Mej2(1, 1), '%.4f'), ' ', num2str(Mej2(1, 2), '%.4f'), ' ', num2str(Mej2(1, 3), '%.4f'), '⌉'])
disp(['Student:  Mej = |', num2str(Mej2(2, 1), '%.4f'), ' ', num2str(Mej2(2, 2), '%.4f'), ' ', num2str(Mej2(2, 3), '%.4f'), '|'])
disp(['                ⌊', num2str(Mej2(3, 1), '%.4f'), ' ', num2str(Mej2(3, 2), '%.4f'), ' ', num2str(Mej2(3, 3), '%.4f'), '⌋'])
disp(['                ⌈', num2str(0.6164), ' ', num2str(0.2889), ' ', num2str(0.3274), '⌉'])
disp(['Solution: Mej = |', num2str(0.2889), ' ', num2str(0.6365), ' ', num2str(0.3476), '|'])
disp(['                ⌊', num2str(0.3274), ' ', num2str(0.3476), ' ', num2str(0.6750, '%.4f'), '⌋'])
disp(' ')
disp(['                ⌈', num2str(Kej2(1, 1), '%.4f'), '  ', num2str(Kej2(1, 2), '%.4f'), '  ', num2str(Kej2(1, 3), '%.4f'), '⌉'])
disp(['Student:  Kej = | ', num2str(Kej2(2, 1), '%.4f'), ' ', num2str(Kej2(2, 2), '%.4f'), '  ', num2str(Kej2(2, 3), '%.4f'), '|'])
disp(['                ⌊ ', num2str(Kej2(3, 1), '%.4f'), '  ', num2str(Kej2(3, 2), '%.4f'), ' ', num2str(Kej2(3, 3), '%.4f'), '⌋'])
disp(['                ⌈', num2str(-0.3371), '  ', num2str(0.1686), '  ', num2str(0.1686), '⌉'])
disp(['Solution: Kej = | ', num2str(0.1686), ' ', num2str(-0.4214), '  ', num2str(0.2528), '|'])
disp(['                ⌊ ', num2str(0.1686), '  ', num2str(0.2528), ' ', num2str(-0.4214), '⌋'])
