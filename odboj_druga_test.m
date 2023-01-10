function [new_point, new_dir, index] = odboj_druga_test(point,dir)

% poli = @(r) [(sx^2 + sy^2), (2*(sx * x + sy * y)), (x^2 + y^2 - r^2)]
% Razpisana verzija

poli = @(r) [norm(dir)^2, 2*dot(point,dir), (norm(point)^2 - r^2)];
% Polinom, ki enkodira presečišče nove točke in enotske krožnice s polmerom r

roots1 = lfilter(roots(poli(1)), @(x) x > 1e-5);
coef = min(roots1);
index = 1;

new_point = coef * dir' + point;
% Nova točka je premaknjena za koeficient * smer.

% Novo smer računamo tako, da najprej dobimo kot premice skozi novo točko.
normalized = new_point / norm(new_point,2);
cost = normalized(1);
sint = normalized(2);

T = [cost -sint; sint cost] * [1 0; 0 -1] * [cost sint; -sint cost];
% Matrika zrcaljenja točke čez premico skozi izhodišče in presečišče 
new_dir = - T * dir;
% S tem dobimo vektor smeri novega žarka, obrnjen, da res dobimo odboj

end

