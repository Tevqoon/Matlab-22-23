format long

disp("Prva naloga")
h = 0.1;
x1 = arrayfun(@(i) i / 10, 0:9);
y1 = [2 8 4 7 4 7 6 2 8 3];

jy = jacobi(y1, @g2, 0, 300);
disp(jy(1))

yn = newton(@f2, @Jf2, y1', 0, 3);
disp(yn(1))

disp("Druga naloga")

b = 322/101;
x2 = arrayfun(@(i) i / 10, 0:20);
f2 = @(x) cos(b*x^2 + 3/7);
y2 = arrayfun(f2, x2);

% Prvi del
p2 = @(x) polyval(polyfit(x2, y2, 4), x);
disp(p2(1))

% Drugi del
df2 = @(x) -2*b*x*sin(3/7 + b * x^2);
disp(polyval(polyfix(x2, y2, 4, 1, f2(1), 1, df2(1)), 0))

% Polyfix dela samo za least squares, ne da se mi implementirat naše
% minimizacije.

disp("Tretja naloga")

c = 46/101;
A3 = [5 6 4 2;
      6 1 c 3;
      4 c 5 7;
      2 3 7 1];
x3 = (1/2) * ones(4, 1);
e = potencna(A3, x3, 1e-2, Inf);
disp(e);

e = hotelling(A3, x3, 1e-2, Inf, 3);
disp(e)

disp("Četrta naloga")

d = 541/103;
A4 = zeros(5,5);
for i = 1:5
    for j = 1:5
        if j > i - 2
            A4(i,j) = (d - i + j)^2;
        end
    end
end

[x,y] = size(A4);
X = A4;
for i = 1:5
    s = X(x,y);
    [Q, R] = qr(X - s * eye(x,y));
    X = R * Q + s * eye(x,y);
end
mu = X(x,y);
disp(mu)

x4 = (1/2) * [1 0 0 0 0]';
x41 = inverzna_potencna(A4, x4, mu, 3, 0);

disp(norm(x41, 1))
