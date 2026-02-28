function u=Newton(L)
  N = 100;
  h=1/N;
  D2=-2*diag(ones(N-1,1))+diag(ones(N-2,1),-1)+diag(ones(N-2,1),1);
  D2=-D2/((h*L)^2);
  u =ones(N-1,1);
  change = 1; it = 0;
  while change > 1e-15                  
    J=D2+2*diag(u)-diag(ones(N-1,1));%calculate jacobi matrix
    M=u;
    for i=1: N-1
       M(i)=u(i)*u(i);
    end
    F=D2*u+M-u;
    deltau=-J\F;
    unew=u+deltau; %obtain unew
    change = norm(unew-u,inf);
    u = unew; it = it+1;
  end
  u=[0;u;0];



  