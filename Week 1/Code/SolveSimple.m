% Function to solve dydx=xy.
function SolveSimple(y0)

[x,y] = ode45(@deriv,[0,1],y0);

plot(x,y,'b');

function dydx = deriv(x,y)

dydx = x*y;
