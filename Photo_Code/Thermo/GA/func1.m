function uu=func1(xx,Aged_stretch,Aged_stress,Aged_stretch1,Aged_stress1,Aged_stretch4,Aged_stress4)
nos = 10;
u1 =  func2(Aged_stretch,Aged_stress);
u2 = func3(Aged_stretch1,Aged_stress1);
%u3 = func4(Aged_stretch2,Aged_stress2);
%u4 = func5(Aged_stretch3,Aged_stress3);
u5 = func6(Aged_stretch4,Aged_stress4);
u6 = Aging(xx,30*24*3600,333,Aged_stretch);
u7 = Aging(xx,10*24*3600,353,Aged_stretch1);
%u8 = Aging(xx,30*24*3600,353,Aged_stretch2);
%u9 = Aging(xx,10*24*3600,368,Aged_stretch3);
u10 = Aging(xx,30*24*3600,368,Aged_stretch4);

%0.2*(u3-u8).^2 + 0.2* (u4-u9).^2
u = 0.3*(u1-u6).^2 + 0.3* (u2-u7).^2 + 0.4* (u5-u10).^2;
uu = (1/nos)*sum(u);
end


