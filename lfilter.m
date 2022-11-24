function [retlist] = lfilter(list, predicate)
    retlist = [];   
    next = 1;
    for i = 1:length(list)
        if predicate(list(i))
            retlist(next) = list(i);
            next = next + 1;
        end
    end
    retlist = retlist';
end