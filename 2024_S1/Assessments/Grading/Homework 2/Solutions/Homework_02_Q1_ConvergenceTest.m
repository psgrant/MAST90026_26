% Number of simulations to run. h is halved each time
nSim = 24;
N = 2.^linspace(0, nSim - 1, nSim);

p = 5;
q = -1;
r = -6;
a = -3;
b = 1;
alpha = 14;
beta = 2;

h = zeros(1, nSim);
err1i = zeros(1, nSim);
t1i = zeros(1, nSim);
% run the simulations and only save the error and time taken
for i = 1:nSim
    tic
    [h(i), err1i(i), ~, ~] = Homework_02_Q1i(N(i), p, q, r, a, b, alpha, beta);
    t1i(i) = toc;
end

err1ii = zeros(1, nSim);
t1ii = zeros(1, nSim);
% run the simulations and only save the error and time taken
for i = 1:nSim
    tic
    [~, err1ii(i), ~, ~] = Homework_02_Q1ii(N(i), p, q, r, a, b, alpha, beta);
    t1ii(i) = toc;
end

% loglog plot of the error.
subplot(1, 2, 1)
loglog(h, err1i, h, err1ii)
title('Error vs grid size')
xlabel('h')
ylabel('||E||_2')
legend('First order', 'Second order')
subplot(1, 2, 2)
loglog(t1i, err1i, t1ii, err1ii)
title('Error vs execution time')
xlabel('t')
ylabel('||E||_2')
legend('First order', 'Second order')