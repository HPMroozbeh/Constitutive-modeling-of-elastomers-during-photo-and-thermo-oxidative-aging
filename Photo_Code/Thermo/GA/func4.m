function uu2 = func4(Aged_stretch2,Aged_stress2)
% T=60 and t =30 days
nos = 10;
cc2 = linspace(Aged_stretch2(1,1),(length(Aged_stretch2)-1)/4+1,nos);
uu2 = interp1(Aged_stretch2,Aged_stress2,cc2);
end