function y = funG(x,u)
    y = zeros(100,1);
    y(1) = x(1)^2 + exp(-u(1));
    for i=2:100
        y(i) = u(i) - u(i-1) + x(i)^2 + exp(-u(i));
    end
end