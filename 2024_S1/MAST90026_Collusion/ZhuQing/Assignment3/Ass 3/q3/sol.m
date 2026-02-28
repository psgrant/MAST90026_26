function [t1,solution1]=sol(L,N)
h=L/N;
u=zeros(N,1);
tspan=[0 100];
x=h:h:L;
for i=1:N
    u(i)=1-abs(sign(x(i)-1))/2-sign(x(i)-1)/2;
end
[t1,solution1]=ode45(@(t,y) odefun(t,y,N,L),tspan,u);
[X,T] = meshgrid(x,t1);
figure(1)
mesh(X,T,solution1)
xlabel('x')
ylabel('T')
zlabel('solution')





