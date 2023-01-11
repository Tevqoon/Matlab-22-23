function [e, x] = hotelling(A, x0, tol, N, redukcije)

B = @(z) A * z;
for i = 1:redukcije
    [e,x] = potencna(B, x0, tol, N);
    B = @(z) B(z) - e * (x' * z) * x;
end

end
