function exa=TRUE2(N,p)
[~,x]=cheb(N-1);
x=flip(x);
exact=zeros(N,1);
for i=1:N
    exact(i)=1/(p+(x(i))^4);
end
exa=exact;