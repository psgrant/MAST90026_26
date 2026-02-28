function u = FEMNO(node, elem, Dfunc, fFunc, gFunc,Q)
N = size(node,1); NT=size(elem,1);
A = zeros(N, N); b = zeros(N,1);
for j = 1:NT
P1 = node(elem(j, 1), :);
P2 = node(elem(j, 2), :);
P3 = node(elem(j, 3), :);
Q1 = elem(j, 1);
Q2 = elem(j, 2);
Q3 = elem(j, 3);
KE = elem_stiffnew(P1, P2, P3, Dfunc); % compute local stiffness matrix
FE = elem_loadnew(P1, P2, P3, fFunc,Q1,Q2,Q3, Q); % compute local load vector
ME = elem_massnew(P1, P2, P3,Q1,Q2,Q3, Q); % compute local mass matrix
A(elem(j, :), elem(j, :)) = A(elem(j, :), elem(j, :)) + KE+ME ;
b(elem(j, :), 1) = b(elem(j, :), 1) + FE ;
end

isLeftBnd = abs(node(:, 1)-0) < 10^(-10);
isRightBnd = abs(node(:, 1)-1) < 10^(-10);
isButtomBnd = abs(node(:, 2)-0) < 10^(-10);
isTopBnd = abs(node(:, 2)-1) < 10^(-10);
isBndNode = false(N, 1);
isBndNode(isLeftBnd) = true;
isBndNode(isRightBnd) = true;
isBndNode(isButtomBnd) = true;
isBndNode(isTopBnd) = true;
bndNode = find(isBndNode);
freeNode = find(~isBndNode);
v = zeros(N, 1);
v(bndNode) = gFunc(node(bndNode, 1),node(bndNode,2));
b = b-A*v;
v(freeNode) = A(freeNode, freeNode)\b(freeNode);
u=v;






