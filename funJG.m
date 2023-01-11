function Y = funJG(u)
    Y = zeros(100,100);
    for i=1:100
        Y(i,i) = - 1 - exp(-u(i));
    end
    Y = Y + diag(ones(99,1), -1);
end