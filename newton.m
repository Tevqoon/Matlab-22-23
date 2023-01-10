function [x,X,k] = newton(f,Jf,x0,tol,N)
% Opis:
%  newton izvede Newtonovo metodo za resevanje sistema nelinearnih enacb
%
% Definicija:
%  [x,X,k] = newton(f,Jf,x0,tol,N)
%
% Vhodni podatki:
%  f    preslikava, ki doloca nelinearni sistem f(x) = 0,
%  Jf   Jacobijeva matrika preslikave f,
%  x0   zacetni priblizek (stolpec),
%  tol  toleranca absolutnega ujemanja dveh zaporednih priblizkov (privzeta
%       vrednost je 1e-10),
%  N    maksimalno stevilo korakov iteracije (privzeta vrednost je 100)
%
% Izhodni podatki:
%  x    koncni priblizek za resitev sistema f(x) = 0,
%  X    tabela vseh izracunanih priblizkov,
%  k    stevilo opravljenih korakov

if nargin < 5
    N = 100;
    if nargin < 4
        tol = 1e-10;
    end
end

X = [x0 zeros(length(x0),N)];
for k = 1:N
   deltax = -Jf(X(:,k))\f(X(:,k));
   X(:,k+1) = X(:,k) + deltax;
   if norm(deltax) < tol
       break;
   end
end

X = X(:,1:k+1);
x = X(:,end);

end