function [x, X, n_final] = racionalni_newton(d, u, rho, i, n_max, eps, x0)
%RACIONALNI_NEWTON Nicla sekularne funkcije na (d_i+1, d_i)

% d je vektor diagonalcev iz odstevanja
% i je indeks desnega krajisca, imamo torej (d(i+1), d(i))
% S tem i predstavlja tudi ravno i-to niclo po velikosti

% Zaenkrat deluje le za prave intervale. ne na koncu

psi1 = @(x) rho * sum(u(1:i) .^ 2 ./ (d(1:i) - x));
dpsi1 = @(x) rho * sum(u(1:i) .^ 2 ./ (d(1:i) - x).^2);
psi2 = @(x) rho * sum(u(i+1:end) .^ 2 ./ (d(i+1:end) - x));
dpsi2 = @(x) rho * sum(u(i+1:end) .^ 2 ./ (d(i+1:end) - x).^2);

n_final = n_max;
X = zeros(n_max + 1, 1);
X(1) = x0;

for n = 1:n_max
    if n > 1 && abs(X(n) - X(n-1)) < eps
        n_final = n;
        X(n+1) = X(n);
        break;
    end
    fr1 = 1 / (d(i) - X(n));
    c1 = [fr1 1; fr1^2 0] \ [psi1(X(n)); dpsi1(X(n))];
    fr2 = 1 / (d(i+1) - X(n));
    c2 = [fr2 1; fr2^2 0] \ [psi2(X(n)); dpsi2(X(n))];
    
    h1 = @(l) c1(1) / (d(i) - l) + c1(2);
    h2 = @(l) c2(1) / (d(i+1) - l) + c2(2);
    h = @(l) 1 + h1(l) + h2(l);

    X(n+1) = fzero(h, X(n));

end

x = X(n_final + 1);

