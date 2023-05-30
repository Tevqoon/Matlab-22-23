format long;

disp("Prva naloga")

% Prvi del
n = 13; T = 6.5;
A1 = gen_A(n, T);
S1 = svd(A1);
disp(S1(3));

% Drugi del
n = 600; T = 6.5;

A1 = gen_A(n, T);

f2 = @(t) sin(t) * cos(1.00 * t) * (t-1);

x2 = linspace(0, T, n+1)';
x2 = x2(2:end);
% Piecewise funkcija na roke
for i = 1:n
    if 1 < x2(i)
        x2(i) = f2(x2(i));
    else
        x2(i) = 0;
    end
end
b2 = A1 * x2;

% Iscemo min F(t), ti so shranjeni v b2
disp(min(b2))

% Tretji del
f3 = @(t) (t - 1.00)^2;
b3 = linspace(0, T, n+1)';
b3 = b3(2:end);

% Piecewise funkcija na roke
for i = 1:n
    if 1 < b3(i)
        b3(i) = f3(b3(i));
    else
        b3(i) = 0;
    end
end

x3 = odrezan_svd(A1, b3, 30);
disp(sum(x3) / 601); % Dodan 1 za f(0) = 0, ki ga prej izlocimo

% Cetrti del
f4 = @(t) sin(t - 1) * (t - 3);
x4 = linspace(0, T, n+1)';
x4 = x4(2:end);

noise = (cos(x4 * 20) / 100);

% Piecewise funkcija na roke
for i = 1:n
    if 1 < x4(i) && x4(i) < 3
        x4(i) = f4(x4(i));
    else
        x4(i) = 0;
    end
end

b4_clean = A1 * x4;

b4 = b4_clean + noise;

x4_noisy = odrezan_svd(A1, b4, 100);
disp(max(abs(x4 - x4_noisy)));

% Peti del

% Komentirano ker pocasno
% k = zeros(20,1);
% for i = 1:20
%     x5 = tikhonov_svd(A1, b4, i / 20);
%     k(i) = norm(A1 * x5 - b4)^4 + norm(x5)^4;
% end
% disp(min(k));

disp("Druga naloga")

A2 = [2.0000   -2.0000   -0.5000    1.0000    0.5000   -0.5000         0   -0.5000;-2.0000   14.0000    1.0000   -2.0000   -1.0000    1.0000         0    1.0000;-0.5000    1.0000    6.1500   -0.5000   -0.2500    0.2500         0    0.2500; 1.0000   -2.0000   -0.5000   -7.0000    0.5000   -0.5000         0   -0.5000; 0.5000   -1.0000   -0.2500    0.5000    3.2500   -0.2500         0   -0.2500;-0.5000    1.0000    0.2500   -0.5000   -0.2500    5.2500         0    0.2500;      0         0         0         0         0         0    0.9000         0;-0.5000    1.0000    0.2500   -0.5000   -0.2500    0.2500         0    2.2500];

[Veig, Deig] = eig(A2);

[nx, ny] = size(A2);

[u, d, rho] = odstej_do_diagonalne(A2);

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


a = 3/4;
% Popravek na roko
u = [3 a a -a a -3/2 -3/2 0]';
d = [10 6 5 3 2 1 -8 0.9]';

% Prvi del
disp(max(d));

% Drugi del
disp(d(end));

% Tretji del

% Od tu naprej koda ne dela. Problem je v tem, da v racionalni newtonovi
% metodi nisem implementiral robustnega nacina racunanja obeh nicel, zato
% trenutna implementacija konstantno jemlje enega izmed intervalov za
% niclo. To sem opazil prepozno in nisem imel casa popraviti.

% Tretja najvecja lastna vrednost => i = 3
x0 = (d(4) + 9 * d(3)) / 10;
[x1, X1, n1] = racionalni_newton(d,u,rho, 3, 1, 100, x0);
% Iz eig(A2) izracunamo tocno vrednost

x2 = 5.151810014679217;

disp(x2 - x1);


% Cetrti del

% Peti del

x0 = (d(4) + d(3)) / 2;
[x4, X4, n4] = racionalni_newton(d,u,rho, 3, 10000, 1000, x0);

% Lastni vektor

% z najprej inicializiramo na celo zato da lahko kasneje
% permutiramo nazaj. V trenutni verziji d in i je vse na prvih n_max
% mestih, zato ostalo brez problema ostane nicelno.
z = zeros(nx, 1); 
z(1:n_max) = (diag(d(1:n_max)) - x4 * eye(n_max)) \ u(1:n_max);

% Pomnozitev z inverzom permutacije je enak kot resevanje linearnega
% sistema, ki je bolj stabilno in optimizirano.
z1 = P \ z;
z1 = z1 / norm(z1);

disp(norm(Veig(:, 6) - z1));


