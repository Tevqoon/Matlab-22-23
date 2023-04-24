function [x] = odrezan_svd(A, b, odrez)
%ODREZAN_SVD Summary of this function goes here
%   Detailed explanation goes here
[U, S, V] = svd(A);

% Odrezani singularni razcep
U2 = U(:, 1:odrez);
S2 = S(1:odrez, 1:odrez);
V2 = V(:, 1:odrez);

x = V2 * (S2 \ (U2' * b));

end

