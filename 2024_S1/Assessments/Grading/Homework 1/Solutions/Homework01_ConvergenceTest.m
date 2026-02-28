% Number of simulations to run. h is halved each time
nSim = 20;
n = 2.^linspace(0, nSim - 1, nSim);

eps = 0.3;

h = zeros(1, nSim);
err2bi = zeros(1, nSim);
t2bi = zeros(1, nSim);
% run the simulations and only save the error
for i = 1:nSim
    tic
    [h(i), err2bi(i), ~, ~] = Homework01_2bi(n(i), eps);
    t2bi(i) = toc;
end

err2bii = zeros(1, nSim);
t2bii = zeros(1, nSim);
% run the simulations and only save the error
for i = 1:nSim
    tic
    [~, err2bii(i), ~, ~] = Homework01_2bii(n(i), eps);
    t2bii(i) = toc;
end

% loglog plot of the error.
subplot(1, 2, 1)
loglog(h, err2bi, h, err2bii)
title('Error vs grid size')
xlabel('h')
ylabel('||E||_2')
legend('Central', 'Central Upwind')
subplot(1, 2, 2)
loglog(t2bi, err2bi, t2bii, err2bii)
title('Error vs execution time')
xlabel('t')
ylabel('||E||_2')
legend('Central', 'Central Upwind')