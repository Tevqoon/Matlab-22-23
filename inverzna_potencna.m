function [e,x,k] = inverzna_potencna(A, x0, tol, N)

if nargin < 4, N = Inf; end
if nargin < 3, tol = 1e-10; end

x = x0/norm(x0);
e = x' * A * x;
k = 0;
while k < N
    y = (A - e * eye(length(x0))) \ x;
    x = y / norm(y);
    e = x' * (A - e * eye(length(x0))) * x;
    k = k + 1;
end
end