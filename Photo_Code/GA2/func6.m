function uu4 = func6(Aged_stretch4,Aged_stress4)
% T=80 and t = 30
nos = 10;
cc4 = linspace(1,(length(Aged_stretch4)-1)/4+1,nos);
uu4 = interp1(Aged_stretch4,Aged_stress4,cc4);
end