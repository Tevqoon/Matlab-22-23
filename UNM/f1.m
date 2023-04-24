function y = f1(x)
% Vzamemo seznam x-ov, vrnemo seznam y-ov.

n = length(x);
y = zeros(n, 1);

for i = 1:n
    y(i) = g1(x, i) - x(i);
end

end

