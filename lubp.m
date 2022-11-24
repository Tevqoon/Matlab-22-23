function [L,U] = lubp(A)
% Opis:
%  lubp izracuna LU razcep A = LU matrike A brez pivotiranja
%
% Definicija:
%  [L,U] = lubp(A)
%
% Vhodni podatek:
%  A    kvadratna matrika (z nesingularnimi vodilnimi podmatrikami)
%
% Izhodna podatka:
%  L    spodnje trikotna matrika z enkami na diagonali,
%  U    zgornje trikotna matrika

n = size(A,1);
L = eye(n);
for j = 1:n-1
    for i = j+1:n
        L(i,j) = A(i,j)/A(j,j);
        for k = j+1:n
            A(i,k) = A(i,k) - L(i,j)*A(j,k);
        end
    end
end
U = triu(A);

end