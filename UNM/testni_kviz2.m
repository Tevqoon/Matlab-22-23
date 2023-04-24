format long

disp("Prva naloga")

x1 = arrayfun(@(i) i * h1, 1:100);
h = 0.01;

[yj, kj] = jacobi(x1, @g1, 1e-3, 1000);
disp(yj(100))
[ys, ks] = seidl(x1, @g1, 1e-3, 1000);
disp(ys(100))
[yn, Xn, kn] = newton(@f1, @Jf1, x1', 1e-1, 1000);
disp(yn(100))

disp("Druga naloga")
f2 = @(x) x * sin(3*x);
df2 = @(x) 3 * x * cos(3 * x) + sin(3 * x);
x2 = arrayfun(@(i) i/5, 0:10);
y2 = arrayfun(f2, x2);

% Prvi del
p2 = @(x) polyval(polyfit(x2, y2, 3), x);
disp(norm(arrayfun(@(x) f2(x) - p2(x), x2), "inf"))

% Drugi del
x25 = x2(6);
y25 = y2(6);

p2f = @(x) polyval(polyfix(x2, y2, 3, x25, y25), x);
disp(norm(arrayfun(@(x) f2(x) - p2f(x), x2), "inf"))

% Tretji del
x2l = x2(1:6);
y2l = y2(1:6);
x2r = x2(7:11);
y2r = y2(7:11);
p2fl = @(x) polyval(polyfix(x2l, y2l, 3, x25, y25, x25, df2(x25)), x);
p2fr = @(x) polyval(polyfix(x2r, y2r, 3, x25, y25, x25, df2(x25)), x);

p2lr = @(x) (x <= 1) * p2fl(x) + (x > 1) * p2fr(x);
disp(norm(arrayfun(@(x) f2(x) - p2lr(x), x2), "inf"))

disp("Tretja naloga")
A3 = delsq(numgrid("C", 10));
l3 = length(A3);
x3 = ones(l3, 1)/norm(ones(l3, 1));

% Prvi del
e = potencna(A3, x3, 1e-5, Inf);
disp(e)

% Drugi del
e = 1 / potencna(inv(A3), x3, 0, 5);
disp(e)

% Tretji del
e = hotelling(A3, x3, 0, 3000, 5);
disp(e)

disp("Četrta naloga")
A4 = [4 -8 8 4;
      2 6 -2 2;
      2 -2 9 -2;
      4 8 -4 0];
B4 = hess(A4);
C4 = qr_iter(B4, 1e-8);

% Prvi del
l41 = max(diag(C4));
disp(l41)

% Drugi del
l44 = min(diag(C4));
disp(l44)

% Tretji del
x4 = (1/2) * [1 1 1 1]';
x41 = inverzna_potencna(A4, x4, l41, Inf, 1e-6);
x42 = inverzna_potencna(A4, x4, l44, Inf, 1e-6);
disp(norm(x41,"inf") + norm(x42, "inf"))
