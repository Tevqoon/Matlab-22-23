function ret = map(list, f)
    l = length(list);
    ret = zeros(1,l);
    for i = 1:l
        ret(i) = f(list(i));
    end
end