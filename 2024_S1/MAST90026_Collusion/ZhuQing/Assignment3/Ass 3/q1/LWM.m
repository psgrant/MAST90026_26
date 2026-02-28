function solution = LWM(K, N)
    h = 2 / N;
    k = 5 / K;
    v = k / h;
    t = 0:k:5;

    % Define matrix A
    A = createMatrixA(N, v);

    % Initial condition
    x = (1:N)' * h;
    u0 = initialCondition(x);

    % Time-stepping
    uu = timeStepping(A, u0, K, N);

    % Exact solution
    exact = exactSolution(x, t, K, N);

    % Plot the numerical and exact solutions
    plotSolutions(x, uu, exact, K);

    solution = uu;
end

function A = createMatrixA(N, v)
    A1 = diag(-1 * ones(N-1,1), -1) + diag(ones(N-1,1), 1);
    A1(1, N) = -1;
    A1(N, 1) = 1;
    A1 = v * A1 / 2;

    A2 = diag(-2*ones(N,1)) + diag(ones(N-1,1), -1) + diag(ones(N-1,1), 1);
    A2(1, N) = 1;
    A2(N, 1) = 1;
    A2 = v * v * A2 / 2;

    A = diag(ones(N, 1)) - A1 + A2;
end

function u0 = initialCondition(x)
    u0 = abs(sign(x - 0.2)) / 2 + sign(x - 0.2) / 2 - abs(sign(x - 0.4)) / 2 - sign(x - 0.4) / 2;
end

function uu = timeStepping(A, u0, K, N)
    uu = zeros(N, K + 1);
    uu(:, 1) = u0;
    for j = 2:K + 1
        uu(:, j) = A * uu(:, j-1);
    end
end

function exact = exactSolution(x, t, K, N)
    exact = zeros(N, K + 1);
    for j = 1:K + 1
        for i = 1:N
            exact(i, j) = abs(sign(x(i) - t(j) - 0.2)) / 2 + sign(x(i) - t(j) - 0.2) / 2 - abs(sign(x(i) - t(j) - 0.4)) / 2 - sign(x(i) - t(j) - 0.4) / 2;
        end
    end
end

function plotSolutions(x, uu, exact, K)
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
