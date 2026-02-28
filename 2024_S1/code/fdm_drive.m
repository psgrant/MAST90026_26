uexact = @(x) sin(pi*x) + x;
f = @(x) pi^2*sin(pi*x);

N = 2.^(3:6);
maxIt = length(N);
err = zeros(1,maxIt);
err2 = zeros(1,maxIt);
err3 = zeros(1,maxIt);
for k = 1:maxIt
    node = (0:1/N(k):1)';
    u = fdm(node, f, 0, 1);
    node2 = (legslb(N(k)+1)+1)/2;
    u2 = fdm(node2, f, 0, 1);
    err(k) = max(abs(u-uexact(node)));
    err2(k) = max(abs(u2-uexact(node2)));
    node3 = node;
    node3(2:end-1,1) = node3(2:end-1,1) + 1/(N(k)+1)*rand(N(k)-1,1);
    node3 = sort(node3);
    u3 = fdm(node3, f, 0, 1);
    err3(k) = max(abs(u3-uexact(node3)));
end

loglog(N, err);
hold on
loglog(N, err2);
loglog(N, err3);
loglog(N, 1./N.^2, '-*');
legend('uniform', 'nonuniform', 'randon', 'line with slope 2');