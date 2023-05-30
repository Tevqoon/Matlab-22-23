function [P_out, dP_out] = fun_vzmeti(alpha_n, X, e, P0, dP0, t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[x, y] = size(X);
P_out = P0;
dP_out = dP0;

for i = 1:y
    P_out = P_out + alpha_n(i) * exp(e(i) * t) * X(:, i);
    dP_out = dP_out + alpha_n(i) * exp(e(i) * t) * X(:, i) * e(i);
end

end

