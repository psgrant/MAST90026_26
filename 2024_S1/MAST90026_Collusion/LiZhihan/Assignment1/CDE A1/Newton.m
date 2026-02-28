%using finite difference(FD) discretization and Newton's method for solving nonlinear equations.
%sets the number of discretization points N and the domain length L.
N = 100;
L=40;
h=1/N;%step size

%scales the second derivative matrix by (hL)^2 to account for the domain size.
D2=-2*diag(ones(N-1,1))+diag(ones(N-2,1),-1)+diag(ones(N-2,1),1);
D2=-D2/((h*L)^2);
u =ones(N-1,1);%initializes the solution vector u
change = 1; it = 0;% initializes the change (used for convergence) to 1 and the iteration counter it to 0.

%executes the Newton-Raphson iteration 
while change > 1e-15                  
 J=D2+2*diag(u)-diag(ones(N-1,1));%calculate jacobi matrix
 M=u;
 for i=1: N-1
       M(i)=u(i)*u(i);
end
 F=D2*u+M-u;%calculates the residual vector F of the nonlinear equation system.
 deltau=-J\F;%solves the linear system to update the solution.
 unew=u+deltau; %obtain unew
 change = norm(unew-u,inf);
 u = unew; it = it+1;%increments the iteration counter.
end
plot(u)

  