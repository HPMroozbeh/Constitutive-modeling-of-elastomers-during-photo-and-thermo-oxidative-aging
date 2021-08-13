function [child1, child2] = cross_over(parent1, parent2,dimension,stringlength)
% this function will do the cross over from two parents and produce two
% children
ran = floor(rand(1,dimension)*(stringlength-2))+1;

% the following will do the cross over for each variable and not for the
% whole binary string
for i=1:dimension
    child1((i-1)*stringlength+1:i*stringlength) = [parent1((i-1)*stringlength+1:...
        (i-1)*stringlength+ran(i)) parent2((i-1)*stringlength+1+ran(i):i*stringlength)];
    child2((i-1)*stringlength+1:i*stringlength) = [parent2((i-1)*stringlength+1:...
        (i-1)*stringlength+ran(i)) parent1((i-1)*stringlength+1+ran(i):i*stringlength)];
end
end
