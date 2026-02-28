NSim = 30;
N = linspace(2, NSim+1, NSim);

u0 = 1;
u1 = 2;

t = zeros(1, NSim);
err = zeros(1, NSim);
for i = 1:NSim
    tic

    NN = N(i);
    [D, x] = cheb(NN);

    A = D^2 - D - 6*eye(NN + 1);
    A(1, 1) = 1;
    A(1, 2:NN+1) = 0;
    A(NN+1, NN+1) = 1;
    A(NN+1, 1:NN) = 0;

    b = 30*x;
    b(1) = u1;
    b(NN+1) = u0;

    u = A\b;
    err(i) = sqrt((pi/NN)*sum(sqrt(1 - x(2:end-1).^2).*(u(2:end-1) - Lab01_ExactSolution(x(2:end-1))).^2));
    t(i) = toc;
end

subplot(1, 2, 1)
loglog(N, err)
xlabel('N')
ylabel('||E||_2')
subplot(1, 2, 2)
loglog(t, err)
xlabel('t')
ylabel('||E||_2')
xlim([min(t)/10 10*max(t)])

