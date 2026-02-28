function u=che(p)
L = chebop(-1, 1);%define the chebop domain 
L.op = @(x,u) diff(u,2)*(p+x^2)^2 -8*x^2*u;%define the exact form of chebop 
L.lbc = 1/(p+1); L.rbc = 1/(p+1);%define the operator boundary conditions
u = L\(-2); %solve the BVP problem and get the polynomial Chebfun with number of cheb points 
plot(u), grid on %plot the solution
residue=norm(L(u)+2); %calculate the residue in ||Lu-f||
disp(residue)
exact=zeros(201,1);
numsol=zeros(201,1);
xx=(-1:0.01:1);%give uniform mesh of spacing 0.01 in [-1,1]
for i=1:201
    exact(i)=1/(p+xx(i)*xx(i));
    numsol(i)=u(xx(i));
end %define the exact solution on [-1,1]
difference=numsol-exact;
error=norm(difference,'inf');%calculate the global error in inf norm
disp(error)





