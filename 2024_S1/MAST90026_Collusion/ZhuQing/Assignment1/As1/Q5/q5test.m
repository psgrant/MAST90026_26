node = -1:0.01:1;
N=length(node);
[pt, qt, ft] = Assignment_1_Q5_BVP_Functions(node);
[D, q, f] = toSelfAdjointForm(node, pt, qt, ft);
elem=[linspace(1,N-1,N-1);linspace(2,N,N-1)]';
u=BvpFE(node,elem,D,q,f,2,0);
plot(node,u)
