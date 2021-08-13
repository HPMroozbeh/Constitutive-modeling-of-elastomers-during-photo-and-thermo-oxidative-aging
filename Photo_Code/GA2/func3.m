function uu1 = func3(Aged_stretch1,Aged_stress1)
% T=45 and t=60days
nos = 10;
cc1 = linspace(Aged_stretch1(1,1),(length(Aged_stretch1)-1)/4+1,nos);
uu1 = interp1(Aged_stretch1,Aged_stress1,cc1);
end