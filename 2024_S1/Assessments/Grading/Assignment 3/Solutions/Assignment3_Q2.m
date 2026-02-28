% We first make k much smaller than h so the error is dominated by the
% spatial discretisaion

T = 10000;
Nlist = [100, 200, 500, 1000];

h = zeros(size(Nlist));
err_CN_h = zeros(size(Nlist));
err_TR_h = zeros(size(Nlist));
for i = 1:length(Nlist)
    N = Nlist(i);
    % Solve using Crank-Nicolson
    [x, u, h(i), ~] = Assignment3_Q2_CrankNicolson(N, T);
    % Note we use a p2 norm to get around the x = 1/2 issue as discussed in
    % the marking scheme
    err_CN_h(i) = sqrt(h(i)*sum((u - exp(-pi^2)*cos(pi*x)).^2));

    % Solve using TR-BDF2 and calculate error
    [x, u, ~, ~] = Assignment3_Q2_TRBDF2(N, T);
    err_TR_h(i) = sqrt(h(i)*sum((u - exp(-pi^2)*cos(pi*x)).^2));
end
subplot(2, 2, 1)
loglog(h, err_CN_h)
xlabel('h')
ylabel('Grid 2-norm')
title('Fixed k, Crank-Nicolson')

subplot(2, 2, 2)
loglog(h, err_TR_h)
xlabel('h')
ylabel('Grid 2-norm')
title('Fixed k, TR-BDF2')



% Next we make the spatial discretisation large and vary the time step 
N = 10000;
hFixed = 1/(N - 1);
Tlist = [100, 200, 500, 1000];

k = zeros(size(Tlist));
err_CN_k = zeros(size(Tlist));
err_TR_k = zeros(size(Tlist));
for i = 1:length(Tlist)
    T = Tlist(i);
    [x, u, ~, k(i)] = Assignment3_Q2_CrankNicolson(N, T);
    % Note we use a p2 norm to get around the x = 1/2 issue as discussed in
    % the marking scheme
    err_CN_k(i) = sqrt(hFixed*sum((u - exp(-pi^2)*cos(pi*x)).^2));

    % Solve using TR-BDF2 and calculate error
    [x, u, ~, ~] = Assignment3_Q2_TRBDF2(N, T);
    err_TR_k(i) = sqrt(hFixed*sum((u - exp(-pi^2)*cos(pi*x)).^2));
end
subplot(2, 2, 3)
loglog(h, err_CN_k)
xlabel('k')
ylabel('Grid 2-norm')
title('Fixed h, Crank-Nicolson')

subplot(2, 2, 4)
loglog(h, err_CN_k)
xlabel('k')
ylabel('Grid 2-norm')
title('Fixed h, TR-BDF2')