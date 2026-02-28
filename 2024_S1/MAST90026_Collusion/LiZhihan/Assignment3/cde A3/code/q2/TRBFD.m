function [solution, error] = TRBFD(t, N, K)
    h = 1 / N;
    k = t / K;
    
    % Construct matrix A
    A = diag(-2 * ones(N + 1, 1)) + diag(ones(N, 1), -1) + diag(ones(N, 1), 1);
    A(1, 2) = 2;
    A(N + 1, N) = 2;
    A = A / h^2;
    
    % Construct matrices F, B, and G
    F = diag(ones(N + 1, 1)) - k * A / 4;
    B = diag(ones(N + 1, 1)) + k * A / 4;
    G = diag(ones(N + 1, 1)) - k * A / 3;
    
    % Initialize x0 and u
    x0 = (0:N)' * h;
    u = cos(pi * x0);
    
    % Initialize solution matrix uu
    uu = zeros(N + 1, K + 1);
    uu(:, 1) = u;
    
    % Time-stepping loop
    for i = 2:K + 1
        v = B * u;
        w = F \ v;
        u = G \ ((4 * w - u) / 3);
        uu(:, i) = u;
    end
    
    solution = uu';
    
    % Compute the exact solution
    exact = zeros(N + 1, K + 1);
    for i = 1:N + 1
        for j = 1:K + 1
            exact(i, j) = cos(pi * x0(i)) * exp(-pi^2 * t * (j - 1) / K);
        end
    end
    
    % Calculate error at the midpoint and final time step
    mid_point = N/4+1;
    error = abs(solution(K + 1, mid_point) - exact(mid_point, K + 1));
end




