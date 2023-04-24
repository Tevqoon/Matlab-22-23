function y = g2(x, i)

a = 53/101;

if i == 10
    y = 1 / x(10) + exp(-a);
else
    y = 1 / x(i)  + exp(-x(i + 1));
end

end
