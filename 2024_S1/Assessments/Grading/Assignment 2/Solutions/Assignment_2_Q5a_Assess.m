clear
close all
clc

% Test element
x1 = 1;
y1 = 1;
x2 = 2;
y2 = 3;
x3 = 1;
y3 = 7/2;
x4 = 0;
y4 = 3/2;

P1 = [x1, y1];
P2 = [x2, y2];
P3 = [x3, y3];
P4 = [x4, y4];

% Call students function
[JR, bR] = rectangle_map(P1, P2, P3, P4);

% Display student and correct result
disp(['Student:  bR = [', num2str(bR(1), '%.2f'), '; ', num2str(bR(2), '%.2f'), ']'])
disp(['Solution: bR = [', num2str(1, '%.2f'), '; ', num2str(1, '%.2f'), ']'])
disp(' ')
disp(['Student:  JR = ⌈', num2str(JR(1, 1), '%.2f'), ' ', num2str(JR(1, 2), '%.2f'), '⌉'])
disp(['               ⌊', num2str(JR(2, 1), '%.2f'), '  ', num2str(JR(2, 2), '%.2f'), '⌋'])
disp(['Solution: JR = ⌈', num2str(1, '%.2f'), ' ', num2str(-1, '%.2f'), '⌉'])
disp(['               ⌊', num2str(2, '%.2f'), '  ', num2str(0.5, '%.2f'), '⌋'])