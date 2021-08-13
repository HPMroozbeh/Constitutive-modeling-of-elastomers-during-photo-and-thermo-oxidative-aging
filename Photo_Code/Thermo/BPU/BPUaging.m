clear;
clc;
format long
% aging time
t=80*24*3600;
%reference temperature
temp = 333;
% desired temperature
Des_temp = 353;
% decay function parameters
T1 = 76.8;
T2 = 160.33;
%Constants
E1 = 46302;
R = 8.314;
I = 1;
%alpha = 0.72;
alpha = 1;
%beta = 0.039;
beta = 1;
gamma = 1;
%Time-temperature superposition
af1 = exp((E1/R)*((1/temp)-(1/Des_temp)));
% Network peobability fitted parameters
% reletive chain length
Rbar_s = 4;
R_1 = 3.3;
R_2 = 0.913;
R_I = 0;
Rbar_b = (Rbar_s) - (R_1 - (R_I/2))*(1-exp(-I^alpha*(t/T2)))+ (R_2 + (R_I/2))*(1-exp(-I^beta*(t/T2)));
q_s = 1.0;
q_1 = 0.7;
q_2 = 0.36;
q_I = 0;
q_b = (q_s) - (q_1 - (q_I/2))*(1-exp(-I^alpha*(t/T2)))+ (q_2 + (q_I/2))*(1-exp(-I^beta*(t/T2)));

nue = 1.02;
n_max = 100;
N_s00= 19.04*10^6;

% calculating the area under probability functions for network s & b
area_s = ProbIntegral(nue*Rbar_s,n_max,Rbar_s,q_s);
area_b = ProbIntegral(nue*Rbar_b,n_max,Rbar_b,q_b);
% making sure that there is nosegment generation
C1=((FirstMoment(nue*Rbar_s,n_max,Rbar_s,q_s)*area_b)/(FirstMoment(nue*Rbar_b,n_max,Rbar_b,q_b)*area_s));
N_b0inf = C1*N_s00;
% scale parameter between chi and lambda
C=0.2;
% micro_sphere
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
Lambda_rel=1;
Lambda_s=1;
q=1;
P_s=zeros([3 3 1]);
T=zeros([3 3 1]);
P_pp_s=zeros(1,21);
P_bT=zeros([3 3 1]);
P_pp_bT=zeros(1,21);
P_b=zeros([3 3 1]);
P_pp_b=zeros(1,21);
cunt = 1;
% 12 is about the same amount for chi = 6
if Lambda_s>=1 && Lambda_s<=12
    while Lambda_s<12
    for i=1:21
        Lambda_Dmax(i)=max(1,eval(Lambda_D(i)));
        nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
        nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
     P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
            *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
        P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
    end
    T(:,:,q)=((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))));
     q=q+1;
     cunt = cunt+1;
     Lambda_s=Lambda_s+0.5;
     P_b(:,:,q)=0;
    end
end

% Changing stress to MPa
for i=1:(q-1)
stress(i)=T(1,1,i)/10^6;
end
Micro_strain(1) = 0;
Micro_strain(2) = 0.5;
for i = 3:(cunt-1)
Micro_strain(i) = Micro_strain(i-1)+0.5;
end
%Micro_strain = 0:0.5:11;
Macro_strain = (1-(0.2^(1/3))).*(Micro_strain+1)+(0.2^(1/3)-1);
% Experimental data point
%num = xlsread('BPU_95.xlsx');
% strain_exp=num(:,1)./100;
% stress_exp=num(:,5);
or = [0.93 0.69 0.13];
pu = [0.71 0.27 1];
gr = [ 0.39 0.83 0.07];
plot(Macro_strain,stress,'-','Color',gr,'LineWidth',3)
hold on
%plot(strain_exp,stress_exp,'o','Color',gr,'LineWidth',3,'MarkerSize',9)
xlim([0 6])
ylim([0 10])
 title('\theta = 95')
 xlabel('Strain [%]') 
 ylabel('Stress [MPa]')
%legend({'BPU Aged for 0 days'},'Location','northwest')

% error calculation
% square_er = 0;
% max_error = 0;
% 
% i =2;
% 
% for j= 1:(q-1)
% if  strain_exp(i) <= Macro_strain(j) && i <= length(strain_exp)-2
%     
%     square_er = square_er + (((stress(j-1)+(((strain_exp(i)-Macro_strain(j-1))/(Macro_strain(j)-Macro_strain(j-1)))*(stress(j)-stress(j-1))))-stress_exp(i))/stress_exp(i))^2;
%     max_error = max(max_error , sqrt(square_er));
%     i = i+1;
% end
% end
% 
% Error = sqrt(square_er/length(stress_exp))*100
% max_error*100