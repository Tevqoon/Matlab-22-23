function x = inverzna_potencna(A, x, l, N, eps)
    k = 0;
    while k < N && norm(A * x - l * x) >= eps 
        y = (A - l * eye(size(A))) \ x;
        x = y/norm(y);
        k = k + 1;
    end
end