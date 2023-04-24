function [outlist] = shift_list(list,n)
% Puts the last n values on the front of the list

len = length(list);

outlist = cat(2, list((len - n + 1):len), list(1:(len - n)));

end

