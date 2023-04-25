format long;

disp("Prva naloga")

% Prvi del
n = 10; T = 7;
A1 = gen_A(n, T);
disp(A1(3,2) + A1(7,5));

% Drugi del
n = 500; T = 7;
A1 = gen_A(n, T);

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
b2 = A1 * x2;
disp(b2(400));

% Tretji del
x3 = odrezan_svd(A1, b2, 10);

disp(norm(x2 - x3, "inf"))

% Cetrti del
noise = (sin((1:n)' / 10) / 300);
b3 = b2 + noise;

x4 = odrezan_svd(A1, b3, 20);
disp(norm(x2 - x4, "inf"))

% Peti del
alfa = 1 / 8;

x5 = tikhonov_svd(A1, b3, alfa);
disp(norm(A1 * x5 - b3, 2));


disp("Druga naloga")

A2 = [1.2500 -0.5000 1.0000 -1.0000 0.2500 0.2500 0 0.2500;
    -0.5000 11.0000 -2.0000    2.0000   -0.5000 -0.5000 0 -0.5000;
    1.0000 -2.0000 9.0000 -4.0000 1.0000 1.0000 0 1.0000;
    -1.0000 2.0000 -4.0000 -2.0000 -1.0000 -1.0000 0 -1.0000;
    0.2500 -0.5000 1.0000 -1.0000 3.2500 0.2500 0 0.2500;
    0.2500 -0.5000 1.0000 -1.0000 0.2500 4.2500 0 0.2500;
    0 0 0 0 0 0 1.0000 0;
    0.2500 -0.5000 1.0000 -1.0000 0.2500 0.2500 0 2.2500];

[nx, ny] = size(A2);

[u, d, rho] = odstej_do_diagonalne(A2);
% Koncna norma napake
% disp(norm(A2 - diag(d) - rho * (u * u')));

% Trace
disp(sum(d));

% Sortiranje & brisanje
[d, Pn] = sort(d, "descend");

% Funkcija sort hrani permutacijo svojih elementov, ki jo lahko uporabimo
% kot zacetno permutacijo. Tu jo spremenimo v permutacijsko matriko.
P = eye(nx);
P = P(Pn, :);
u = P * u;

n_max = nx; % Gledamo, koliko slabih vrednosti imamo
i = 1;
while i <= n_max
    if abs(u(i)) < 1e-10 || (i < n_max && d(i) == d(i + 1))
        % Menjamo.
        [u(i), u(n_max)] = deal(u(n_max), u(i));
        [d(i), d(n_max)] = deal(d(n_max), d(i));
        [P(i, :), P(n_max, :)] = deal(P(n_max, :), P(i, :));
        n_max = n_max - 1;
    else
        i = i + 1;
    end
end

% Sekularna funkcija
f = @(x) 1 + sum(rho * u .^ 2 ./ (d - x));
disp(f(7));

% Racionalna Newtonova metoda
x0 = (d(1) + d(2)) / 2;
[x1, X1, n1] = racionalni_newton(d,u,rho, 1, 1, "inf", x0);
disp(x1);

% Malo pointless met opcijo za epsilon, ce mamo itak max korake glede na to
% da ne gre inicializirat infinite arraya. Bolje bi bilo shranjevati samo
% zadnji dve vrednosti direktno in ju potem menjati.
[x2, X2, n2] = racionalni_newton(d,u,rho,1, 100000, 1e-15, x0);
disp(x2);

% Lastni vektor

% z najprej inicializiramo na celo zato da lahko kasneje
% permutiramo nazaj. V trenutni verziji d in i je vse na prvih n_max
% mestih, zato ostalo brez problema ostane nicelno.
z = zeros(nx, 1); 
z(1:n_max) = (diag(d(1:n_max)) - x1 * eye(n_max)) \ u(1:n_max);

% Pomnozitev z inverzom permutacije je enak kot resevanje linearnega
% sistema, ki je bolj stabilno in optimizirano.
z1 = P \ z;
z1 = z1 / norm(z1);

disp(abs(z1(1)))

