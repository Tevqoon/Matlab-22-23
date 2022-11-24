function [new_point, new_dir, index] = odboj_druga(point,dir)

% poli = @(r) [(sx^2 + sy^2), (2*(sx * x + sy * y)), (x^2 + y^2 - r^2)]
% Razpisana verzija

poli = @(r) [norm(dir)^2, 2*dot(point,dir), (norm(point)^2 - r^2)];
% Polinom, ki enkodira presečišče nove točke in enotske krožnice s polmerom r

roots3 = lfilter(roots(poli(3)), @(x) x > 0);
roots4 = lfilter(roots(poli(4)), @(x) x > 0);
% Zanimajo nas samo pozitivne ničle. S tem poberemo tudi kompleksne.

% Koeficient, ki ga iščemo, je najmanjša pozitivna ničla izmed obeh
% polinomov. Najprej prevermimo, če kakšen nima pozitivne.
if isempty(roots4)
    index = 1;
    coef = min(roots3);
elseif isempty(roots3)
    index = 2;
    coef = min(roots4);
else
    [coef, index] = min([min(roots3), min(roots4)]);
end

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

