% Test if student function is working correctly
f = @(x)15*sin(x.^3);
xj = 0.76;
xjp1 = 0.892;
FEj = elem_load(f, xj, xjp1);
disp(['Student solution FEj = [', num2str(FEj(1), 4), '; ', num2str(FEj(2), 4), ']']);
disp('Correct solution FEj = [0.4927; 0.5677]')