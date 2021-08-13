clear;
clc;

% Network peobability fitted parameters
% reletive chain length
Rbar_s = 1.37; 
% an indicator for crosslink density
q_s = 0.75;
% we won't need it in this section but can be an indicator of lost chain
beta_s = 1;
nue=1.001;
n_max=100;
N_s00=8*10^5;
% scale parameter between chi and lambda
C=0.2;
% decay function parameters
A1 = 1.4*10^-9;
A2 = 8.3*10^-9;
T1 = 5.2*10^5;
T2 = 2*10^7;
%Constants
E = 80000;
R = 8.314;
%reference temperature
temp = 353;
% desired temperature
Des_temp = 333;
%Time-temperature superposition
af = exp((E/R)*((1/temp)-(1/Des_temp)));
% area under probability function
area = ProbIntegral(nue*Rbar_s,n_max,Rbar_s,q_s);
% Micro_sphere constants
xc=[1,0,0,0.7071067812,0.7071067812,0.7071067812,0.7071067812,0,0'...
    ,0.3879073041,0.3879073041,0.3879073041,0.3879073041,0.3879073041'...
    ,0.3879073041,0.3879073041,0.3879073041,0.8360955967,0.8360955967'...
    ,0.8360955967,0.8360955967];
yc=[0,1,0,0.7071067812,-0.7071067812,0,0,0.7071067812,0.7071067812'...
    ,0.3879073041,0.3879073041,-0.3879073041,-0.3879073041,0.8360955967'...
    ,0.8360955967,-0.8360955967,-0.8360955967,0.3879073041'...
    ,0.3879073041,-0.3879073041,-0.3879073041];
zc=[0,0,1,0,0,0.7071067812,-0.7071067812,0.7071067812,-0.7071067812'...
    ,0.8360955967,-0.8360955967,0.8360955967,-0.8360955967,0.3879073041'...
    ,-0.3879073041,0.3879073041,-0.3879073041,0.3879073041,-0.3879073041'...
    ,0.3879073041,-0.3879073041];
wc=[0.0265214244,0.0265214244,0.0265214244,0.0199301476,0.0199301476'...
    ,0.0199301476,0.0199301476,0.0199301476,0.0199301476,0.0250712367'...
    ,0.0250712367,0.0250712367,0.0250712367,0.0250712367,0.0250712367'...
    ,0.0250712367,0.0250712367,0.0250712367,0.0250712367,0.0250712367'...
    0.0250712367];
syms Lambda_s Lambda_max Lambda_y
F_s=[Lambda_s,0,0;0,1/sqrt(Lambda_s),0;0,0,1/sqrt(Lambda_s)];
Lambda_D=sym(zeros(1,21));
Lambda_Dmax=zeros(1,21);
for i=1:21
    D=[xc(i),yc(i),zc(i)];
    Lambda_D(i)=sqrt(D*transpose(F_s)*F_s*transpose(D));
end
% Calculating lambda based on chi = 1.5
Lambda_rel = ((1.5-C^(1/3))/(1-C^(1/3)));
Lambda_s = Lambda_rel;
q=1;
P=zeros([3 3 1]);
T=zeros([3 3 1]);
P_pp=zeros(1,21);
I=eye(3);
num = xlsread('rel_test_at_50_and_60.xlsx');
t=num(:,1);
for j=1:length(t)
    for i=1:21
        Lambda_Dmax(i)=max(1,eval(Lambda_D(i)));
        nmin=(nue*(Lambda_Dmax(i))*Rbar_s);
     P_pp(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nue*Rbar_s,nmin,n_max,Rbar_s,q_s,beta_s)'...
            *NIntegral(nmin,n_max,Rbar_s,q_s,eval(Lambda_D(i)),area);
        P(:,:,q)=P(:,:,q)+((P_pp(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
    end
    T(:,:,q)=((1/(A1+A2))*(A1*exp(-E/(R*temp)*af*(t(j)/T1))+A2*exp(-E/(R*temp)*af*(t(j)/T2))))*(P(:,:,q)-(P(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1)));
     q=q+1;
     P(:,:,q)=0;
end
%A1*exp(-A2*exp(-E/(R*temp))*t(j))
%exp(-t(j)/(Tau_s*((log(t(j)))^pow)))
%A1*exp(-E/(R*temp)*t(j)*A2)
%A1*((R*temp*A2)/(E*(t(j)+T2)))
%((R*temp*A1)/(E*t(j)))^(1/19)
%((R*temp*(exp(-E/(R*temp)*t(j))))/(E*t(j)))^(1/19)
%((R*temp*log10(t(j)))/(E*t(j)))^(1/19)
%((R*temp*(t(j))^0.6)/(E*t(j)))^(1)
%(A1*exp(-((af*t(j))/T2)-zM2*(af*t(j))^2/2)+B*exp(-(af*t(j))/T2))
%(A1*(af*t(j))^A2+B)
%(A1*exp(-E/(R*temp)*t(j)*A2)+A3*((R*temp*A4)/(E*(t(j)))))
%(A1*exp(-E/(R*temp)*af*t(j)*A2)+A3*((R*temp)/(E*af*(t(j)))))
%(A1*exp(-E/(R*temp)*af*(t(j)/T1))+A2*exp(-E/(R*temp)*af*(t(j)/T2)))
for i=1:length(T)
stress(i)=(T(1,1,i))/T(1,1,1);
end

exp_stress = num(:,2)/num(1,2);
% error calculation
square_er = 0;
max_error = 0;
for i = 1: length(t)
    square_er = square_er + ((stress(i)-exp_stress(i))/exp_stress(i))^2;
    max_error = max(max_error , sqrt(square_er));
end

Error = sqrt(square_er/length(t))*100
max_error

plot(t,stress,'blue')
%plot(log10(t),stress,'blue')
hold on
plot(t,exp_stress,'red')
%plot(log10(t),num(:,2),'red')
