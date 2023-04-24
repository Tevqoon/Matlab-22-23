function y = f2(x)

n = length(x);
y = zeros(n, 1);

for i = 1:n
    y(i) = x(i) - g2(x, i);
end

end