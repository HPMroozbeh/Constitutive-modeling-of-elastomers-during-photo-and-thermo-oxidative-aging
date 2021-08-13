function uu5 = func7(Aged_stretch5,Aged_stress5)
% T=80 and t = 30
nos = 10;
cc5 = linspace(1,(length(Aged_stretch5)-1)/4+1,nos);
uu5 = interp1(Aged_stretch5,Aged_stress5,cc5);
end