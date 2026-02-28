function [pt, qt, ft] = Assignment_1_Q5_BVP_Functions(node)
    assert(~(min(node) < -1) && ~(max(node) > 1), 'x values must lie on the interval [-1, 1]')
    
    xSpline = [-1; -0.9; -0.5; -0.1; 0.1; 0.2; 0.4; 1];
    ptSpline = [1; 2; 3; 2; 1; 1; 2; 1];
    qtSpline = [0; 0; 0.5; 1; 1; 0.5; 0; 0];
    ftSpline = [1; -1; -2; -1; 1; 1; -1; 1];

    pt = spline(xSpline, ptSpline, node);
    qt = spline(xSpline, qtSpline, node);
    ft = spline(xSpline, ftSpline, node);
end