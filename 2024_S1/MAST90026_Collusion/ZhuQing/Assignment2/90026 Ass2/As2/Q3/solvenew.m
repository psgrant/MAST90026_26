function [it,sol]=solvenew(f,g,D,elem,node)
Q=zeros(length(node(:,1))*length(node(:,1)),1);
change=1;
count=1;
u=FEMNO(node, elem, D, f, g,Q);
while change > 1e-10 
    Q=u;
    unew=FEMNO(node, elem, D, f, g,Q);
    change=norm(unew-u,'inf');
    u=unew;
    count=count+1;
end
sol=u;
it=count;


