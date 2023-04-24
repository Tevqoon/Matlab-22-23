function [y, k] = seidl(x0, fun, tol, N)

dx = inf;
k = 0;
n = length(x0);

while dx >= tol && k < N
    xstari = x0;
    for i = 1:n
        x0(i) = fun(x0, i);
    end
    dx = norm(xstari - x0, 'inf');
    k = k + 1;
end

y = x0;
end
