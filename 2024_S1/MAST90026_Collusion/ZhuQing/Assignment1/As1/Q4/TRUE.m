function exa=TRUE(node,elem,p)
n=length(node);
x=zeros(n,1);
for i=1:n-1
x(i) = node(elem(i,1));
x(i+1) = node(elem(i,2));
end
exact=zeros(n,1);
for i=1:n
    exact(i)=1/(p+(x(i))^4);
end
exa=exact;

