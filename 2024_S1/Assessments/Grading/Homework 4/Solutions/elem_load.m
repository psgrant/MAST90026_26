function FEj = elem_load(f, xj, xjp1)
    hj = xjp1 - xj;
    x = @(xi)xj + hj*xi;
    N1 = @(xi)1 - xi;
    N2 = @(xi)xi;
    int1 = @(xi)f(x(xi))*N1(xi);
    int2 = @(xi)f(x(xi))*N2(xi);
    xi0 = (1/2)*(1 - sqrt(3/5));

    % The quadrature points/weights
    w0 = 5/18;
    xi1 = 1/2;
    w1 = 8/18;
    xi2 = (1/2)*(1 + sqrt(3/5));
    w2 = 5/18;
    % Calculate the integrals
    FEj = hj*(w0*[int1(xi0); int2(xi0)] + w1*[int1(xi1); int2(xi1)] + w2*[int1(xi2); int2(xi2)]);
end