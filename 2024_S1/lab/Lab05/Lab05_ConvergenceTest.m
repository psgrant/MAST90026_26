% Number of simulations to run. h is halved each time
nSim = 10;
n = 2.^linspace(1, nSim, nSim);

h = zeros(1, nSim);
err = zeros(1, nSim);
t = zeros(1, nSim);
% run the simulations and only save the error
for i = 1:nSim
    tic
    [h(i), err(i), ~, ~] = Lab05_FDM_Poisson_2(n(i));
    t(i) = toc;
end

% loglog plot of the error.
subplot(1, 2, 1)
loglog(h, err)
title('Error vs grid size')
xlabel('h')
ylabel('||E||_2')
subplot(1, 2, 2)
loglog(t, err)
title('Error vs execution time')
xlabel('t')
ylabel('||E||_2')