%To Test toSelfAdjointForm
a=0;
b=1;
N=100;
x=a:1/N:b;

%Define the original non-self-adjoint functions
p = @(x) x+6;
q = @(x) -(x+1);
f = @(x) x^2+3;
pt=zeros(1,N+1);
qt=zeros(1,N+1);
ft=zeros(1,N+1);
for i=1:length(x)
    pt(i) = p(x(i));
    qt(i) = q(x(i));
    ft(i) = f(x(i));
end

[D, q, f] = toSelfAdjointForm(x, pt, qt, ft);