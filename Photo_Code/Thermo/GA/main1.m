% this code is based on genetic algorithm for the minimum of 
% a function called func1 
% the termination is set according to the number of generations ngen
clc;
clear;
close all;


ngen = input('enter the number of generations?');
pop_size = input('enter the size of population?');
stringlength = input('enter the string length?');
dimension = input('enter the number of variables?');
var_bound = input('enter a vector containing the bounds for each variable?');
pc = input('enter the cross over probability?');
pm = input('enter the mutation probability?');


% importing the data
Aged_data = xlsread('BPU_6030.xlsx');
Aged_data1 = xlsread('BPU_8010.xlsx');
%Aged_data2 = xlsread('BPU_8030.xlsx');
%Aged_data3 = xlsread('BPU_9510.xlsx');
Aged_data4 = xlsread('BPU_9530.xlsx');
Aged_stretch = Aged_data(:,1)+1;
Aged_stress = Aged_data(:,2);
Aged_stretch1 = Aged_data1(:,1)+1;
Aged_stress1 = Aged_data1(:,2);
% Aged_stretch2 = Aged_data2(:,1)+1;
% Aged_stress2 = Aged_data2(:,2);
%  Aged_stretch3 = Aged_data3(:,1)+1;
%  Aged_stress3 = Aged_data3(:,2);
Aged_stretch4 = Aged_data4(:,1)+1;
Aged_stress4 = Aged_data4(:,2);

% construction of the initial population
pop = round(rand(pop_size,dimension*stringlength));
new_pop = zeros(pop_size,dimension*stringlength);

x_min = zeros(ngen,dimension);

for i=1:ngen
    
   xx=decoding_pop(pop,dimension,stringlength,var_bound);
   % this part calls the objective function evaluation for the population
   % func1 is the name of the objective function
   for i1=1:pop_size
       yy(i1) = func1(xx(i1,:),Aged_stretch,Aged_stress,Aged_stretch1,Aged_stress1,Aged_stretch4,Aged_stress4);
   end
   % ss1 is the vector that contains the probability
   % distribution of the individuals in selection
   % process
   % f_min is the vector of the minimum of func1 in each generation
   
   f_min(i) = min(yy);
   kk = find(yy==min(yy));
   x_min(i,:) = xx(kk(1),:);
   ss1 = yy./sum(yy);
   
   % the Roulete selection is used for the selection process
   sum1 = 0;
   for i2=1:pop_size
       ss2(i2) = ss1(i2)+sum1;
       sum1 = ss2(i2);
   end
   
   for i3=1:floor(pop_size/2)
       aa1=rand;
       aa2=rand;
       for i4=1:pop_size
           if aa1<ss2(i4)
               jj=i4;%index of the first selected parent
               break
           end
       end
       for i5=1:pop_size
           if aa2<ss2(i5)
               kk=i5;%index of the second selected parent
               break
           end
       end
    
    
       % we do cross over on the new selected parents if the cross-over
       % probability is greater than a random number
       if rand<pc
           [ch1,ch2] = cross_over(pop(jj,1:dimension*stringlength),...
         pop(kk,1:dimension*stringlength),dimension,stringlength);
       else
           ch1 = pop(jj,1:dimension*stringlength);
           ch2 = pop(kk,1:dimension*stringlength);
       end
       
       
       % we do mutaion
       ch11 = mutation(ch1, dimension, stringlength, pm);
       ch22 = mutation(ch2, dimension, stringlength, pm);
       
       
       % we substitute the new offsprings in the new population
       new_pop(2*i3-1,:)=ch11;
       new_pop(2*i3,:)=ch22;
   end
   ii = min(yy);
   jj = find(yy==ii);
   ii1 = pop(jj,:);
   new_pop(jj,:) = ii1;
   pop = new_pop;
   
end

figure(1);
hold on
title('Objective function vs. generation number')
xlabel('Generation number')
ylabel('Objective function')
axis([1 inf 0 inf])
plot((1:ngen),f_min)

figure(2)
hold on
nos = 10;
plot(Aged_stretch,Aged_stress,'d','Color','k','LineWidth',3,'MarkerSize',9)
uuu = Aging(x_min(end,:),30*24*3600,333,Aged_stretch);
x1 = linspace(1,(length(Aged_stretch)-1)/4+1,nos);
plot(x1,uuu(1,:),'-','Color','k','LineWidth',3)
plot(Aged_stretch1,Aged_stress1,'*','Color','r','LineWidth',3,'MarkerSize',9)
uuu1 = Aging(x_min(end,:),10*24*3600,353,Aged_stretch1);
x2 = linspace(1,(length(Aged_stretch1)-1)/4+1,nos);
plot(x2,uuu1(1,:),'-','Color','r','LineWidth',3)
% plot(Aged_stretch2,Aged_stress2,'o','Color','b','LineWidth',3,'MarkerSize',9)
% uuu2 = Aging(x_min(end,:),30*24*3600,353,Aged_stretch2);
% x3 = linspace(1,(length(Aged_stretch2)-1)/4+1,nos);
% plot(x3,uuu2(1,:),'-','Color','b','LineWidth',3)
%  plot(Aged_stretch3,Aged_stress3,'p','Color','c','LineWidth',3,'MarkerSize',9)
%  uuu3 = Aging(x_min(end,:),10*24*3600,368,Aged_stretch3);
%  x4 = linspace(1,(length(Aged_stretch3)-1)/4+1,nos);
%  plot(x4,uuu3(1,:),'-','Color','c','LineWidth',3)
plot(Aged_stretch4,Aged_stress4,'s','Color','m','LineWidth',3,'MarkerSize',9)
uuu4 = Aging(x_min(end,:),30*24*3600,368,Aged_stretch4);
x5 = linspace(1,(length(Aged_stretch4)-1)/4+1,nos);
plot(x5,uuu4(1,:),'-','Color','m','LineWidth',3)
title('Stress vs. stretch')
xlabel('Stretch')
ylabel('Stress')
axis([1 inf 0 inf])



