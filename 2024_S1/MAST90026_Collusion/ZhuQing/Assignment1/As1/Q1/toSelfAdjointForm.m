function [D, q, f] = toSelfAdjointForm(x, pt, qt, ft)
temp=cumtrapz(x, pt);
D=exp(temp);
q=-qt.*D;
f=-ft.*D;




