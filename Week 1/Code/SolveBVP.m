% Function to solve y''+x^2 y' - 4y^2 = 0.
function SolveBVP()

solinit = bvpinit([0,1],[0 0]);

sol = bvp4c(@deriv,@bcs,solinit);

x=linspace(0,1,100);
y=deval(sol,x);

plot(x,y(1,:),'b');

function dYdx = deriv(x,Y)

dYdx(1) = Y(2); 
dYdx(2) = 4*Y(1)-3*Y(2);

% boundary conditions y(0)=1, y(1)=1.
function res = bcs(ya,yb)
res = [ ya(1)-1 
        yb(1)-1];
    
% % boundary conditions y(0)=1, y(1)=1.
% function res = bcs(ya,yb)
% res = [ ya(2) 
%         yb(1)-1];