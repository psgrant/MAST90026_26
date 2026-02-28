function u = FEM_Elliptic_2D_Dirichlet(node, elem, Dfunc, qFunc, fFunc, gFunc, quad_order)
N = size(node,1); NT=size(elem,1);
A = zeros(N, N); b = zeros(N,1);
for j = 1:NT
P1 = node(elem(j, 1), :);
P2 = node(elem(j, 2), :);
P3 = node(elem(j, 3), :);
KE = elem_stiff(P1, P2, P3, Dfunc, quad_order); % compute local stiffness matrix
FE = elem_load(P1, P2, P3, fFunc, quad_order); % compute local load vector
ME = elem_mass(P1, P2, P3, qFunc, quad_order); % compute local mass matrix
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






