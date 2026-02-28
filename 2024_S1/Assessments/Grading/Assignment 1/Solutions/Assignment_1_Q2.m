
% N+1 is the number of Chebyshev-Gauss-Lobatto collocation points to use
N = 10;

% Calculate the collocation points and differentiation matrix
[D, x] = cheb(N);


u = BvpSpec(D, pt, qt, ft, alpha, beta);