function [e,x,k] = potencna(A,x0,tol,N)
% Opis:
%  potencna izracuna priblizek za dominantni lastni par matrike A
%
% Definicija:
%  [e,x,k] = potencna(A,x0,tol,N)
%
% Vhodni podatki:
%  A    kvadratna matrika ali mnozenje z matriko A,
%  x0   zacetni priblizek za dominantni lastni vektor matrike A
%  tol  toleranca napake priblizka za dominantni lastni par (privzeta
%       vrednost je 1e-10),
%  N    maksimalno stevilo korakov iteracije (privzeta vrednost je 200)
%
% Izhodni podatki:
%  e    priblizek za dominantno lastno vrednost (Rayleighjev kvocient
%       priblizka za lastni vektor),
%  x    priblizek za dominantni lastni vektor,
%  k    stevilo opravljenih korakov iteracije

if nargin < 4, N = Inf; end
if nargin < 3, tol = 1e-10; end

if ~isa(A,'function_handle')
    if nargin < 2, x0 = rand(size(A,1),1); end
    A = @(x) A*x;
end

x = x0/norm(x0);
Ax = A(x);
e = x'*Ax;
k = 0;
while norm(Ax-e*x) >= tol && k < N
    x = Ax/norm(Ax);
    Ax = A(x);
    e = x'*Ax;
    k = k+1;
end

end