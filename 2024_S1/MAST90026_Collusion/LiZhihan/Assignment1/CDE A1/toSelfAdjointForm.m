%To compute the coefficients in the self-adjoint form of the ODE
function [D, q, f] = toSelfAdjointForm(x, pt, qt, ft)

%Computes the cumulative integral
temp=cumtrapz(x, pt);
D=exp(temp); %Exponentiated to get D(x).
q=-qt.*D;
f=-ft.*D;



