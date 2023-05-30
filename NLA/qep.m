function [Lam,X] = qep(M, C, K)
%QEP Summary of this function goes here
%   Detailed explanation goes here

[x, y] = size(M);

A = [zeros(x,y), eye(x, y); -inv(M) * K, -inv(M) * C];
[V, D] = eig(A);

Lam  = diag(D);
X = V(1:x,:);
end

