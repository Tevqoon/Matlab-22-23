function [lam, x] = newtonNEP(T, dT, lam, x, v, maxKorakov)

% Vhodni podatki

% T : funkcija, katero lastno vrednost iscemo
% dT : odvod funkcije T
% lam : priblizek za lastno vrednost
% maxKorakov : maksimalno stevilo korakov metode
% x : priblizek za lastni vektor
% v : izbrani vektor, da velja v^H x_k = 1

% Izhodni podatki

% lam : priblizek za lastno vrednost
% x : priblizek za lastni vektor

 for i = 1:maxKorakov
    u = T(lam) \ dT(lam) * x;
    lam = lam - (v' * x) / (v' * u);
    x = u / (v' * u);
 end
end
