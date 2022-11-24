function [output_matrix] = K_cetrta(input_list)
% Makes a permutation matrix as wanted by the 4. naloga
len = length(input_list);
output_matrix=ones(len);

for i = 1:len
    output_matrix(:,i) = shift_list(input_list, i - 1);
end

end

