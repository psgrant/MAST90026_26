function [errr,p]=Errr(N)
f=@(x,y) -2*exp(x+y);
Error=zeros(N-9,1);
for k=10:N
hx = 1/(k+1); 
hy = 2/(k+1); 
xx = 0:hx:1;
yy = 0:hy:2;
[x, y] = meshgrid(xx, yy);
Exact=zeros(k+2,k+2);
for i=1:k+2
    for j=1:k+2
        Exact(i,j)=exp(x(i,j)+y(i,j));
    end
end
numsol=FDM2D(f, 0, 1, 0, 2, k, k);
error=Exact-numsol;
error=reshape(error,(k+2)*(k+2),1);
Error(k-9)=norm(error,'inf');
end
errr=Error;
s=10:1:N;
loglog(s,errr,'*')
p=log(Error(N-59)/Error(N-9))/(log(2));





