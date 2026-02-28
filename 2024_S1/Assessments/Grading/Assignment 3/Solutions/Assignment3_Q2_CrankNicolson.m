function [x, u, h, k] = Assignment3_Q2_CrankNicolson(N, T)
% Function takes the arguments N which is the number of grid points (including end but excluding ghost)
% and T which is the number of time points

    % Define the time step and grid size
    h = 1/(N - 1);
    k = 1/T;
    
    % List of x values
    x = linspace(0, 1, N)';
    
    % Set the initial condition and use this to define the boundary conditions
    u = cos(pi*x);
    
    % Create A matrix of the discretised system (note we don't need the vector b because it is all zeros 
    e = ones(N, 1);
    A = spdiags([e/(h^2), -2*e/(h^2), e/(h^2)], -1:1, N, N);
    A(1, 2) = 2/(h^2);
    A(N, N - 1) = 2/(h^2);
    
    % Perform the Crank-Nicolson time stepping
    for i = 1:T
        RHS = (speye(N) + k*A/2)*u;
        LHSMat = speye(N) - k*A/2;
        u = LHSMat\RHS;
    end

end