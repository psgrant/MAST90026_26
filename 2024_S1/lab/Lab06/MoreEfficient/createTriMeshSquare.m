function [node, elem, bndNode, freeNode] = createTriMeshSquare(xS, xE, yS, yE, nx, ny)
    % List of x and y points
    x = linspace(xS, xE, nx);
    y = linspace(yS, yE, ny);
    % Create grid of x and y values
    [X, Y] = meshgrid(x, y);
    % Turn grid into list we can use to represent node
    X = X(:);
    Y = Y(:);
    node = [X, Y];
    % Generate the triangulation from the nodes
    elem = delaunay(node);


    eps = 1e-14;
    isLeftBnd = abs(node(:, 1) - xS) < eps;
    isRightBnd = abs(node(:, 1) - xE) < eps;
    isBottomBnd = abs(node(:, 2) - yS) < eps;
    isTopBnd = abs(node(:, 2) - yE) < eps;
    isBndNode = false(nx*ny, 1);
    isBndNode(isLeftBnd) = true;
    isBndNode(isRightBnd) = true;
    isBndNode(isBottomBnd) = true;
    isBndNode(isTopBnd) = true;
    bndNode = find(isBndNode);
    freeNode = find(~isBndNode);

end