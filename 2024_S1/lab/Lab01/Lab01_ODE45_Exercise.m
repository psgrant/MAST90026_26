tspan = [0 20];
uinit = [1 0];
[t, u] = ode45(@umat, tspan, uinit);
plot(t, u(:, 1))

function dydt = umat(~, u)
    dydt = [u(2); (1-u(1).^2).*u(2) - u(1)];
end