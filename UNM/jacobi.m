function [y, k] = jacobi(x0, fun, tol, N)

dx = inf;
k = 0;
n = length(x0);
x1 = zeros(1, n);

while dx >= tol && k < N
    for i = 1:n
        x1(i) = fun(x0, i);
    end
    dx = norm(x1 - x0, 'inf');
    x0 = x1;
    k = k + 1;
end

y = x1;
end
