function r = rayleigh(x, A)
    r = (x' * A * x) / (x' * x);
end