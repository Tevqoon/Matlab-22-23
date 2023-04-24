function y = g1(x, i)

h = 0.01;

if i == 1
    y = (i * h)^2 + exp(-x(i));
else
    y = x(i - 1) + (i * h)^2 + exp(-x(i));
end

end