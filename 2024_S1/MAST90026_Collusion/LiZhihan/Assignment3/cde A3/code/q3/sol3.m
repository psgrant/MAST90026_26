function [t3,solution3]=sol3(L,N)
h=L/N;
x=h:h:L;
tspan=[0 10];
y0=zeros(N,1);
x0=zeros(N,1);
for i=1:N
    x0(i)=L*i/N;
    if x0(i)<1
    y0(i)=1;
    else
    y0(i)=cos(-(x0(i)-1)/5);
    end
end
[t3,solution3]=ode45(@(t,y) odefun(t,y,N,L),tspan,y0);
%plot 
[X,T] = meshgrid(x,t3);
figure(1)
mesh(X,T,solution3)
xlabel('x')
ylabel('T')
zlabel('solution')
