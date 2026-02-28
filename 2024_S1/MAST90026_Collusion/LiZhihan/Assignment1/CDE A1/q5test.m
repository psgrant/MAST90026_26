%using the finite element method (FEM) with linear basis functions.
%creates a vector of nodes from -1 to 1 with a step size of 0.01
node = -1:0.01:1;
N=length(node);
[pt, qt, ft] = Assignment_1_Q5_BVP_Functions(node);
[D, q, f] = toSelfAdjointForm(node, pt, qt, ft);
elem=[linspace(1,N-1,N-1);linspace(2,N,N-1)]';%Creates an element connectivity matrix elem.
u=BvpFEnew(node,elem,D,q,f,2,0);% finite element solver
plot(node,u)
