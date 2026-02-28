function solution = UP(N, K)
    h = 2 / N;
    k = 5 / K;
    v = k / h;
    t = 0:k:5;
    x = (1:N)' * h;  % Column vector for spatial grid points

    % Initial condition
    u0 = 0.5 * (sign(x - 0.2) - sign(x - 0.4));

    % Coefficient matrix
    A1 = diag(-ones(N, 1)) + diag(ones(N-1, 1), -1);
    A1(1, N) = 1;
    A = diag(ones(N, 1)) + v * A1;

    % Time-stepping
    uu = zeros(N, K + 1);
    uu(:, 1) = u0;
    for j = 2:K + 1
        uu(:, j) = A * uu(:, j - 1);
    end

    % Exact solution
    exact = zeros(N, K + 1);
    for j = 1:K + 1
        exact(:, j) = 0.5 * (sign(x - t(j) - 0.2) - sign(x - t(j) - 0.4));
    end

    % Plot the numerical and exact solutions
    figure;
    plot(x, uu(:, K + 1), 'r', 'LineWidth', 1.5);
    hold on;
    plot(x, exact(:, K + 1), 'g--', 'LineWidth', 1.5);
    xlabel('x');
    ylabel('u(x, t)');
    title('Numerical and Exact Solutions at Final Time Step');
    legend('Numerical Solution', 'Exact Solution');
    grid on;
    hold off;

    solution = uu;
end
