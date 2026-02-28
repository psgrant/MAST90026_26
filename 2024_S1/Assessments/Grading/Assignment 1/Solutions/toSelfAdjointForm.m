function [D, q, f] = toSelfAdjointForm(x, pt , qt , ft)
    D = exp(cumtrapz(x, pt));
    q = -qt.*D;
    f = -ft.*D;
end