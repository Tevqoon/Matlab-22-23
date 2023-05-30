format long;

disp("Prva naloga")
rng(1000);
A = (rand(10));
B = (rand(10));

T = @(lam) A + 1./(lam-2).*B + lam.*eye(size(A));
dT = @(lam) -1./((lam-2).^2).*B + eye(size(A));

% Prvi del
[X, Y] = meshgrid(linspace(-1, 1, 100));
Z = X + 1i*Y;
detT = arrayfun(@(z) abs(det(T(z))), Z);
contourf(X, Y, detT, [0.1, 3e-3]);
% Range je finnicky, we're trying to find a value small enough to show the
% individual zeroes yet large enough to make the ranges visible
% Higher exponents of the absolute value seem to be slightly easier to work
% with
title('|det(T(z))|');

% Drugi del
lam_init = 0.3i;
x_init = [ones(5, 1); zeros(5, 1)];
v = ones(10, 1);
[lam, x] = newtonNEP(T, dT, lam_init, x_init, v, 3);
residual = norm(T(lam) * x);

disp(residual);

% Tretji del
lam_init = 0.6i; % Choose something from the bucket of the eigenvalue
x_init = [ones(5, 1); zeros(5, 1)];
v = ones(10, 1);
[lam, x] = newtonNEP(T, dT, lam_init, x_init, v, 100);
disp(imag(lam));

disp("Druga naloga")

% Define the matrices
M = [1 0 0; 0 2 0; 0 0 3];
K = [2 -1 0; -1 2 -1; 0 -1 1];
C = zeros(size(M));

P = [1.5; 3; 4.5];
% P enacbe ne naredi nehomogene, samo zacetne pogoje da druge

% R pove zacetne polozaje v mirujocem stanju
R = [1; 2; 3]; 

% Initial conditions
v0 = [0; 0; 0];

% Prvi del

disp(norm(R, 1));

% Drugi del

[X, e] = polyeig(K, C, M);
disp(max(abs(e)));

% Tretji del

% compute the alpha_n; q0 - R gives a good coordinate system
    
alpha_n = [X; e.' .* X] \ [P - R; v0];

[q20, v20] = fun_vzmeti(alpha_n, X, e, R, v0, 20); % R daje default polozaje

%disp("Here")
%disp(q(0));

disp(q20);

% Cetrti del

disp(v20);

% Peti del

M_slon = M;
M_slon(3,3) = M_slon(3,3) + 100;
[Xs, es] = polyeig(K, C, M_slon);

alpha_s = [Xs; es.' .* Xs]  \ [q20 - R; v20];

[q30, v30] = fun_vzmeti(alpha_s, Xs, es, R, v0, 10);
disp(q30);

disp("Tretja naloga")

[x, y] = meshgrid(linspace(-1, 1, 100));
M1 = sin(cos(10*x) .* (abs(y) + 1) .^ x);

[x, y, z] = meshgrid(linspace(-1, 1, 100));
M2 = sin(cos(10*x) .* (abs(y) + 1) .^ x) .* z.^2;

% Prvi del
m = 4;
[U, S, V] = svd(M1);
Sm = S;
Sm(m+1:end,:) = 0;
Sm(:,m+1:end) = 0;
M1_approx = U * Sm * V';
diff_M1 = abs(M1 - M1_approx);
max_diff_M1 = max(diff_M1(:));

disp(max_diff_M1);

% Drugi del

% Tenzoriraj
M2_tensor = tensor(M2);

rank_m = 2;
init = {eye(size(M2, 1), rank_m), eye(size(M2, 2), rank_m), eye(size(M2, 3), rank_m)}; % Zacetne

% CP ALS with initial matrices and maximum 4 iterations
M2_approx = cp_als(M2_tensor, 2, "maxiters", 4, "init", init);
M2_app_matrix = double(M2_approx); % Back to a MD array

diff_M2 = abs(M2 - M2_app_matrix);
max_diff_M2 = max(diff_M2(:));

disp(max_diff_M2);

