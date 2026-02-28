% Number of simulations to run. h is halved each time
nSim = 20;
n = 2.^linspace(0, nSim - 1, nSim);

h = zeros(1, nSim);
err = zeros(1, nSim);
t = zeros(1, nSim);
% run the simulations and only save the error
for i = 1:nSim
    tic
    [h(i), err(i), ~, ~] = Lab01_FDM(n(i));
    t(i) = toc;
end

% loglog plot of the error.
subplot(1, 2, 1)
loglog(h, err)
title('Error vs grid size for central finite difference scheme')
xlabel('h')
ylabel('||E||_2')
subplot(1, 2, 2)
loglog(t, err)
xlabel('t')
ylabel('||E||_2')