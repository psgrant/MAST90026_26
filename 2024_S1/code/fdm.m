function u = fdm(node, f, u0, u1)
%% U = FDM(NODE, F, U0, U1) finite difference mehtod on general meshes
%
%   Copyright (C) Hailong Guo
%   08/05/2022

%% find mesh size
N = length(node);
h = diff(node);

%% Form the differential matrix
A = zeros(N-2,N-2);
A(1,1) = 2/(h(1)*h(2));
A(1,2) = -2/(h(2)*(h(1)+h(2)));
for i = 2:N-3
   A(i,i-1) = -2/(h(i)*(h(i)+h(i+1)));
   A(i,i) = 2/(h(i+1)*h(i));
   A(i,i+1) = -2/(h(i+1)*(h(i)+h(i+1)));
end
A(N-2,N-3) = -2/(h(N-2)*(h(N-2)+h(N-1)));
A(N-2,N-2) = 2/(h(N-2)*h(N-1));
rhs = f(node(2:end-1));

%% Impose boundary condition
rhs(1) = rhs(1) + 2*u0/(h(1)*(h(1)+h(2)));
rhs(N-2) = rhs(N-2) + 2*u1/(h(N-1)*(h(N-1)+h(N-2)));

%% Solve the linear system
u = A\rhs;
u = [u0; u; u1];