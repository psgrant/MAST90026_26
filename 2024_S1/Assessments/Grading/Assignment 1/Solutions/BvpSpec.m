function [x, u] = BvpSpec(D,pt,qt,ft,alpha,beta,N)
    
    % Calculate the chebyshev point for the output. Note this is somewhat
    % redundant since they would have already been calculated to get D.
    [~, x] = cheb(N);

    % Construct the collocation matrix. Rows corresponding to boundary
    % conditions are the corresponding row of the identity matrix
    A = D^2 + diag(pt)*D + diag(qt);
    A(1, 1) = 1;
    A(1, 2:N+1) = 0;
    A(N+1, :) = D(N+1, :);

    % Construct the RHS of the linear system
    b = ft;
    b(1) = beta;
    b(N+1) = alpha;
    
    % Solve the collocation linear system
    u = A\b;
end