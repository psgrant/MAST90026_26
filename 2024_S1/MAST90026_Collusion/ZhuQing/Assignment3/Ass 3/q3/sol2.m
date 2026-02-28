function [t2,solution2]=sol2(L,N)
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
    y0(i)=exp(-(x0(i)-1)/3);
    end
end
[t2,solution2]=ode45(@(t,y) odefun(t,y,N,L),tspan,y0);
%plot 
[X,T] = meshgrid(x,t2);
figure(1)
mesh(X,T,solution2)
xlabel('x')
ylabel('T')
zlabel('solution')
