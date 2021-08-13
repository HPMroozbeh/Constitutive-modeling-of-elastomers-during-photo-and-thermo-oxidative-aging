function uu = func2(Aged_stretch,Aged_stress)
% T=45 and t=10days
nos = 10;
cc = linspace(Aged_stretch(1,1),(length(Aged_stretch)-1)/4+1,nos);
uu = interp1(Aged_stretch,Aged_stress,cc);
end