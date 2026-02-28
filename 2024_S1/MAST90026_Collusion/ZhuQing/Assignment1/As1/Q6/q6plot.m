L=[10,20,30,40,50,60];
u1=Newton(L(1));
u2=Newton(L(2));
u3=Newton(L(3));
u4=Newton(L(4));
u5=Newton(L(5));
u6=Newton(L(6));
x=0:1/100:1;
figure;

for i = 1:length(L)
    subplot(2, 3, i);
    u = Newton(L(i));
    plot(x, u);
    title(['Newton Interpolation for L(', num2str(i), ') = ', num2str(L(i))]);
end