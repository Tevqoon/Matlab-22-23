function [x] = tikhonov_svd(A, b, alfa)
%TIKHONOV_SVD Summary of this function goes here
%   Detailed explanation goes here

[U, S, V] = svd(A);

D = S ./ (S .^ 2 + alfa^2);
x = V * (D * (U' * b));

end

