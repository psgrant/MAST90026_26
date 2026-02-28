quad_order=2;
appro=FEM_Elliptic_2D_Dirichlet(node, elem, @(x,y) 1, @(x,y) 0, @(x,y) 17*sin(x).*cos(4*y), @(x,y) sin(x).*cos(4*y), quad_order);
exact=zeros(length(node(:,1)),1);
for i=1:length(node(:,1))
    exact(i)=sin(node(i,1))*cos(4*node(i,2));
end
max_error=norm(appro-exact,inf);
disp(max_error);
[node1,elem1] = uniformrefine(node,elem);
appro1=FEM_Elliptic_2D_Dirichlet(node1, elem1, @(x,y) 1, @(x,y) 0, @(x,y) 17*sin(x).*cos(4*y), @(x,y) sin(x).*cos(4*y), quad_order);
exact1=zeros(length(node1(:,1)),1);
for i=1:length(node1(:,1))
    exact1(i)=sin(node1(i,1))*cos(4*node1(i,2));
end
max_error1= norm(appro1-exact1,inf);
p=log(max_error/max_error1)/log(2);%%compute the order
disp(p)
