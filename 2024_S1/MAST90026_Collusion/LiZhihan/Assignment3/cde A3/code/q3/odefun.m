function dydt = odefun(t,y,N,L)
h=L/N;
dydt = zeros(N,1);
dydt(1) = (-2*y(1)+y(2))/h^(2)+1/h^(2)+y(1)*(1-y(1));
for i=2:N-1
dydt(i) = (-2*y(i)+y(i-1)+y(i+1))/h^(2)+y(i)*(1-y(i));
end
dydt(N)=(2*y(N-1)-2*y(N))/h^(2)+y(N)*(1-y(N));
