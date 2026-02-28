function [exact, solution] = UP1(N, K)
    h = 2 / N;
    k = 5 / K;
    v = k / h;
    t = 0:k:5;
    x = (1:N)' * h;

    % Initialize the A1 matrix
    A1 = diag(-1 * ones(N, 1)) + diag(ones(N - 1, 1), -1);
    A1(1, N) = 1;
    A1 = v * A1;

    % Initialize the A matrix
    A = diag(ones(N, 1)) + A1;

    % Initial condition
    u0 = exp(-10 * (4 * (x - 1)).^2);
    uu = zeros(N, K + 1);
    uu(:, 1) = u0;

    % Time-stepping loop
    for j = 2:K + 1
        uu(:, j) = A * uu(:, j - 1);
    end

    solution = uu;

    % Exact solution
    exact = zeros(N, K + 1);
    for j = 1:K + 1
        exact(:, j) = exp(-10 * (4 * (x - t(j)) - 1).^2);
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
end
