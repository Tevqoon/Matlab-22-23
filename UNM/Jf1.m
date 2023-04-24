function A = Jf1(x)

n = length(x);
A = diag(ones(1, n - 1), -1);

for i = 1:n
    A(i, i) = -exp(-x(i)) - 1;
end

end