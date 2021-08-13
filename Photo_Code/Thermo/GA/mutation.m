function rr=mutation(pop_individual, dimension, stringlength, pm)
% this function mutates a binary individual according to the probability

ran = rand(1,dimension*stringlength);
rr=pop_individual;

for i=1:dimension*stringlength
    if ran(i)<pm
        rr(i) = mod(pop_individual(i)+1,2);
    end 
end
end
