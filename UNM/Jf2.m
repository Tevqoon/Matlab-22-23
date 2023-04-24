function A = Jf2(x)

n = length(x);
A = zeros(n,n);

for i = 1:n
    A(i,i) = 1 + 1 / x(i)^2;
end

for i = 2:n
    A(i - 1, i) = exp(-x(i));
end