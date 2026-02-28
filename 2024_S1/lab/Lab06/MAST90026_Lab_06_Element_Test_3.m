x1 = 0.5;
y1 = 0.5;
x2 = 1;
y2 = 1;
x3 = 0.5;
y3 = 1;

P1 = [x1, y1];
P2 = [x2, y2];
P3 = [x3, y3];

D = @(X)5;
q = @(X)1;
f = @(X)3;

%%% Problem 2: Try these functions after you have it working with constant D, q and f
D = @(X)1 + X(1).*X(2);
q = @(X)exp(sqrt(X(1).^2 + X(2).^2));
f = @(X)sin(X(1)) + cos(X(2));

Kej3 = elem_stiff_3(P1, P2, P3, D);
Mej3 = elem_mass_3(P1, P2, P3, q);
Fej3 = elem_load_3(P1, P2, P3, f);