function uu3 = func5(Aged_stretch3,Aged_stress3)
% T=60 and t = 30
nos = 10;
cc3 = linspace(Aged_stretch3(1,1),(length(Aged_stretch3)-1)/4+1,nos);
uu3 = interp1(Aged_stretch3,Aged_stress3,cc3);
end