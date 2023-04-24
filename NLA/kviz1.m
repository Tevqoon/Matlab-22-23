format long;

disp("Prva naloga")

% Prvi del
n = 10; T = 7;
A1 = gen_A(n, T);
disp(A1(3,2) + A1(7,5));

% Drugi del
n = 500; T = 7;
A2 = gen_A(n, T);

f2 = @(x) (x - 5).^2 .* (x - 3).^2;

x2 = linspace(0, T, n+1)';
x2 = x2(2:end);
% Piecewise funkcija na roke
for i = 1:n
    if 3 <= x2(i) && x2(i) <= 5
        t = x2(i);
        x2(i) = (t - 5)^2 * (t - 3)^2;
    else
        x2(i) = 0;
    end
end
b2 = A2 * x2;
disp(b2(400));

% Tretji del
x3 = odrezan_svd(A2, b2, 10);

disp(norm(x2 - x3, "inf"))

% Cetrti del
noise = (sin((1:n)' / 10) / 300);
b3 = b2 + noise;

x4 = odrezan_svd(A2, b3, 20);
disp(norm(x2 - x4, "inf"))

% Peti del
alfa = 1 / 8;

x5 = tikhonov_svd(A2, b3, alfa);
disp(norm(A2 * x5 - b3, 2));


disp("Druga naloga")
