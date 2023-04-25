function [u, d, rho] = odstej_do_diagonalne(A)
%ODSTEJ_DO_DIAGONALNE find u, D, rho that A2 = diag(D) - rho * (u * u')

[n, ~] = size(A);
% Assume square

% Racunanje ustreznega razcepa za matriko
minimizator_aux = @(u, D, rho) A - diag(D) - rho * (u * u');
minimizator = @(uD) minimizator_aux(uD(1:n), uD(n+1:2*n), uD(2*n+1));
u0D0 = ones(2 * n + 1,1);

% Na roke:
% u = [0.5 -1 2 -2 0.5 0.5 0 0.5];

% Optimization toolbox
ud = fsolve(minimizator, u0D0, optimset("display", "off", "algorithm", "levenberg-marquardt"));
u = ud(1:n);
d = ud(n+1:2*n);
rho = ud(2*n+1);

end

