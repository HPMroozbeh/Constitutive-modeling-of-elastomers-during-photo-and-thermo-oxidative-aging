function uu=func1(xx,Aged_stretch,Aged_stress,Aged_stretch1,Aged_stress1,Aged_stretch2,Aged_stress2,Aged_stretch3,Aged_stress3)
nos = 10;
u1 =  func2(Aged_stretch,Aged_stress);
u2 = func3(Aged_stretch1,Aged_stress1);
u3 = func4(Aged_stretch2,Aged_stress2);
u4 = func5(Aged_stretch3,Aged_stress3);
u6 = Aging(xx,150*24*3600,318,Aged_stretch);
u7 = Aging(xx,30*24*3600,353,Aged_stretch1);
u8 = Aging(xx,10*24*3600,333,Aged_stretch2);
u9 = Aging(xx,10*24*3600,353,Aged_stretch3);

var = dot(u9,u9-u7);
weight = 0;
if var < 0
    weight = 10;
end

u =  0.08* (u1-u6).^2 + 0.4* (u2-u7).^2 + 0.12*(u3-u8).^2 + 0.4*(u4-u9).^2;
uu = (1/nos)*sum(u)+ weight;
end


