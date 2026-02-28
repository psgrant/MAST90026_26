function [h, err, X, Y, uTot] = Lab05_FDM_Poisson_2(N)
    x = linspace(-1, 1, N+2);
    y = linspace(-1, 1, N+2);
    [X, Y] = meshgrid(x, y);
    
    Xin = X(2:N+1, 2:N+1);
    Yin = Y(2:N+1, 2:N+1);
    
    h = 2/(N+1);
    IN = speye(N);
    
    e = ones(N, 1);
    Cmat = (1/h^2)*spdiags([-e, 4*e, -e], -1:1, N, N);
    Dmat = -(1/h^2)*speye(N);
    Emat = spdiags([e, e], [-1,1], N, N);
    
    A = kron(IN, Cmat) + kron(Emat, Dmat);
    f = (5*pi^2/4)*sin(pi*Xin(:)).*cos(pi*Yin(:)/2);
    
    u = A\f;
    u = reshape(u, [N, N]);
    
    uTot = zeros(N+2); 
    uTot(2:N+1, 2:N+1) = u;
    
    figure(1)
    surf(X, Y, uTot)
    xlabel('x')
    ylabel('y')

    uExact = sin(pi*X).*cos(pi*Y/2);
    err = sqrt(trapz(y, trapz(x, (uTot - uExact).^2, 2)));
end