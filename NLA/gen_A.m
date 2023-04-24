function [A] = gen_A(n,T)
    % Create a time grid with n+1 points
    t = linspace(0, T, n+1);
    
    % K(t, tau) je odvisen le od razlike t in tau, ki je ravno enaka 
    % i-ti tocki, ker je t0 = 0
    K = @(razlika) 1/(2*sqrt(pi)) .* (razlika).^(-3/2) .* exp(-1./(4.*razlika));
    razlike = K(t);

    % Bookkeeping za toeplitz, nastavimo se razlike(1), kjer drugace dobimo
    % NaN.
    razlike(1) = 0;
    vrstica = zeros(n+1, 1);
    
    % Toeplitz z nicelno vrstico nam da ravno spodnjo trikotno, ki jo
    % zelimo, na koncu se odstranimo prvo vrstico in stolpec ter pomnozimo
    % s koeficientom.
    A = toeplitz(razlike, vrstica);
    A = A(2:end, 2:end) * (T / (n+1));
end

