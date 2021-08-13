function  xx=decoding_pop(pop,dimension,stringlength,var_bound)

% this function decodes the binary population into the matrix xx which is a 
% population_size*dimension matrix
% var_bound is a vector that contains the bounds of each dimension 


pop_size = size(pop,1);

temp = 2.^fliplr(0:stringlength-1)/(2^stringlength-1);

for i=1:pop_size
    for j=1:dimension
      
        xx (i,j) = sum(pop(i , (j-1)*stringlength+1:j*stringlength).*temp)...
            .*(var_bound(2*j)-var_bound(2*j-1))+var_bound(2*j-1);
    end
end
end

