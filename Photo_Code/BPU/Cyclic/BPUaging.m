clear;
clc;
format long
%%
% Aging time. It is the experimental test time in second.
t = 0 * 24 * 3600;

%Reference temperature. It is the temperature you use as reference. I used
%60 but in general it is better to use highest temperature as reference.
%You choose it once and keep it constant for the whole process.
temp = 333;

% Desired temperature. You change this temperature to the temperature that
% you are using the data for.
Des_temp = 333;

%% decay function parameters

T1 = 80;
T2 = 240;

% E is the activation energy parameter. If you read my paper in polymer
% degradation journal you can know all about it. Their boundary is between
% 20000 - 130000 usually.

E1 = 31634;
E2 = 27333;

% R is constant. You don't need to worry about it.
R = 8.314;

% I is an indication of UV intensity. If we put I = 1, then we would have
% only thermo_oxidative aging condition.

I = 210.84;


% alpha and beta are relative effect of UV in determining the amount of
% cross-linking and chain scission.
alpha = 0.4;
beta = 0.18742;
gamma = 0.31738;
%%
%% Time-temperature superposition
af1 = exp((E1/R)*((1/temp)-(1/Des_temp)));
%% Network peobability fitted parameters
% reletive chain lengths
% R_s is the end to end distance of chain length of unaged material
Rbar_s = 3.39;

% R_1 and R_2 are the effect of cross linking and chain scission on R,
% respectively.
R_1 = 0.6255;
R_2 = 0.0344;

% R_I is the effect of cross linking and chain scission on R.
R_I = -1.14795;

% R_b is the end to end distance of chain length of aged network
Rbar_b = (Rbar_s) - (R_1 + (R_I/2))*(1-exp(-exp(-E2/(R*Des_temp))*I^alpha*(t/T2)))+ (R_2 - (R_I/2))*(1-exp(-exp(-E2/(R*Des_temp))*I^beta*(t/T2)));

% an indicator for crosslink density
% In my paper I have  p_s which is (1-q_s). 
% q boundary is 0< q_s or q_b <1
%since totally aged is usually harder than unaged, q_b usually is smalle
%than q_s. In general higher q means softer material.
q_s = 1.0;
q_1 = 0.0441;
q_2 = 0.01877;
q_I = -0.014916;


q_b = (q_s) - (q_1 + (q_I/2))*(1-exp(-exp(-E2/(R*Des_temp))*I^alpha*(t/T2)))+ (q_2 - (q_I/2))*(1-exp(-exp(-E2/(R*Des_temp))*I^beta*(t/T2)));
% if q_b > 1
%     q_b =1;
% end

% nue is the exact same parameter as the one in 2009 paper of Dr. dargazany
nue = 1.001;

% n_max is the highest number of chains. usually I put it a large number
% and it is not going to affect our result. Becareful large number is
% reletive to R values. if you put R_s as 60, then n_max should be at much
% higher than 100. It should be at least higher than lambda*R* nue.
n_max = 100;

%This parameter only shift verything up or down without changing their
%shape.
N_s00= 13.51*10^6;

% calculating the area under probability functions for network s & b. I
% need them for calculating the total number of chain segments in the
% network.
area_s = ProbIntegral(nue*Rbar_s,n_max,Rbar_s,q_s);
area_b = ProbIntegral(nue*Rbar_b,n_max,Rbar_b,q_b);

% making sure that there is no segment generation. Total number of segments
% of aged network is equal to total number of segments of unaged network
% using these two equations:
C1=((FirstMoment(nue*Rbar_s,n_max,Rbar_s,q_s)*area_b)/(FirstMoment(nue*Rbar_b,n_max,Rbar_b,q_b)*area_s));
N_b0inf = C1*N_s00;

% scale parameter between chi and lambda. It is the same parameter as in
% 2009 paper of Dr. dargazany
C=0.2;
phi= 1/3;



%%
%% micro_sphere
% In this section you don't need to change anything in general , this the
% core of code that compute the stress based on what parameters you gave
% it.
% These are all the 21 constants for directions of microsphere.
% 0.5 is the step size. If you don't need accuracy, you can increase it for
% faster runs. 
Step_size = 0.01;

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
syms Lambda_s Lambda_max 
F_s=[Lambda_s,0,0;0,1/sqrt(Lambda_s),0;0,0,1/sqrt(Lambda_s)];
Lambda_D=sym(zeros(1,21));
Lambda_Dmax=zeros(1,21);
Lambda_maxdir=sym(zeros(1,21));
for i=1:21
    D=[xc(i),yc(i),zc(i)];
    Lambda_D(i)=sqrt(D*transpose(F_s)*F_s*transpose(D));
end

F_smax=[Lambda_max,0,0;0,1/sqrt(Lambda_max),0;0,0,1/sqrt(Lambda_max)];
for i=1:21
    D=[xc(i),yc(i),zc(i)];
    Lambda_maxdir(i)=sqrt(D*transpose(F_smax)*F_smax*transpose(D));
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
cuntL1 = 1;
%% First Loading


%Chi = 1.337568;
%Chi = 2.096;
Chi = 2.106;
Lambda = (Chi - C^(phi))/(1 - C^(phi));

if Lambda_s>=1 && Lambda_s<=Lambda
    while Lambda_s<Lambda
    for i=1:21
        Lambda_Dmax(i)=max(1,eval(Lambda_D(i)));
        nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
        nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
        P_pp_s(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nmin_s,n_max,Rbar_s,q_s)'...
            *NIntegral(nmin_s,n_max,Rbar_s,q_s,eval(Lambda_D(i)));
        P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
            *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
        P_s(:,:,q)=P_s(:,:,q)+((P_pp_s(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
        P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
    end
    T(:,:,q)=(exp(-exp(-E1/(R*temp))*af1*I^gamma*(t/T1)))*((P_s(:,:,q)-(P_s(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))))+((1-(exp(-exp(-E1/(R*temp))*af1*I^gamma*(t/T1)))))'...
        *((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))));
     q = q + 1;
     cuntL1 = cuntL1 + 1;
     Lambda_s = Lambda_s + Step_size;
     P_s(:,:,q) = 0;
     P_b(:,:,q) = 0;
    end
end
%%
%% First Unloading
Lambda_max = Lambda;
Lambda_s = Lambda_s - Step_size;
cuntUnL1 =1;

if Lambda_s>=1 && Lambda_s<=Lambda
    while Lambda_s>=1
    for i=1:21
        Lambda_Dmax(i)=max(1,eval(Lambda_maxdir(i)));
        nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
        nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
        P_pp_s(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nmin_s,n_max,Rbar_s,q_s)'...
            *NIntegral(nmin_s,n_max,Rbar_s,q_s,eval(Lambda_D(i)));
        P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
            *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
        P_s(:,:,q)=P_s(:,:,q)+((P_pp_s(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
        P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
    end
    T(:,:,q)=(exp(-exp(-E1/(R*temp))*af1*I^gamma*(t/T1)))*((P_s(:,:,q)-(P_s(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))))+((1-(exp(-exp(-E1/(R*temp))*af1*I^gamma*(t/T1)))))'...
        *((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))));
    if T(1,1,q)<=0
        break
    end
    Lambda_s    = Lambda_s - Step_size;
    q           = q + 1;
    cuntUnL1    = cuntUnL1 + 1;
    P_s(:,:,q) = 0;
    P_b(:,:,q) = 0;
    end
end

%% Start Second loading
% Lambda_l2 = Lambda;
% Lambda_s = Lambda;
Lambda_s = Lambda_s + Step_size;
%Chi = 1.711434;
%Chi = 3.22;
Chi = 3.234;
Lambda = (Chi - C^(phi))/(1 - C^(phi));

cuntL2 = 1;

if Lambda_s>=1 && Lambda_s<=Lambda
    while Lambda_s<Lambda
    for i=1:21
        Lambda_Dmax(i)=max([1,eval(Lambda_D(i)),eval(Lambda_maxdir(i))]);
        nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
        nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
        P_pp_s(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nmin_s,n_max,Rbar_s,q_s)'...
            *NIntegral(nmin_s,n_max,Rbar_s,q_s,eval(Lambda_D(i)));
        P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
            *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
        P_s(:,:,q)=P_s(:,:,q)+((P_pp_s(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
        P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
    end
    T(:,:,q)=(exp(-exp(-E1/(R*temp))*af1*I^gamma*(t/T1)))*((P_s(:,:,q)-(P_s(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))))+((1-(exp(-exp(-E1/(R*temp))*af1*I^gamma*(t/T1)))))'...
        *((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))));
     q          = q + 1;
     cuntL2     = cuntL2 + 1;
     Lambda_s   = Lambda_s + Step_size;
     P_s(:,:,q) = 0;
     P_b(:,:,q) = 0;
    end
end

%%

%% Start Second Unloading
Lambda_max = Lambda;
Lambda_s = Lambda_s - Step_size;
cuntUnL2 =1;

if Lambda_s>=1 && Lambda_s<=Lambda
    while Lambda_s>=1
    for i=1:21
        Lambda_Dmax(i)=max(1,eval(Lambda_maxdir(i)));
        nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
        nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
        P_pp_s(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nmin_s,n_max,Rbar_s,q_s)'...
            *NIntegral(nmin_s,n_max,Rbar_s,q_s,eval(Lambda_D(i)));
        P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
            *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
        P_s(:,:,q)=P_s(:,:,q)+((P_pp_s(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
        P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
    end
    T(:,:,q)=(exp(-exp(-E1/(R*temp))*af1*I^gamma*(t/T1)))*((P_s(:,:,q)-(P_s(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))))+((1-(exp(-exp(-E1/(R*temp))*af1*I^gamma*(t/T1)))))'...
        *((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))));
    if T(1,1,q)<=0
        break
    end
    Lambda_s    = Lambda_s - Step_size;
    q           = q + 1;
    cuntUnL2    = cuntUnL2 + 1;
    P_s(:,:,q) = 0;
    P_b(:,:,q) = 0;
    end
end

%%

%% Start Third Loading
%Lambda_l3 = Lambda;
%Lambda_s = Lambda;
 Lambda_s = Lambda_s + Step_size;            
%Chi = 2.05626;
%Chi = 4.356;
Chi = 4.35;
Lambda = (Chi - C^(phi))/(1 - C^(phi));

cuntL3 = 1;

if Lambda_s>=1 && Lambda_s<=Lambda
    while Lambda_s<Lambda
    for i=1:21
        Lambda_Dmax(i)=max([1,eval(Lambda_D(i)),eval(Lambda_maxdir(i))]);
        nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
        nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
        P_pp_s(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nmin_s,n_max,Rbar_s,q_s)'...
            *NIntegral(nmin_s,n_max,Rbar_s,q_s,eval(Lambda_D(i)));
        P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
            *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
        P_s(:,:,q)=P_s(:,:,q)+((P_pp_s(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
        P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
    end
    T(:,:,q)=(exp(-exp(-E1/(R*temp))*af1*I^gamma*(t/T1)))*((P_s(:,:,q)-(P_s(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))))+((1-(exp(-exp(-E1/(R*temp))*af1*I^gamma*(t/T1)))))'...
        *((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))));
     q          = q + 1;
     cuntL3     = cuntL3 + 1;
     Lambda_s   = Lambda_s + Step_size;
     P_s(:,:,q) = 0;
     P_b(:,:,q) = 0;
    end
end

%%

%% Start Third Unloading
Lambda_max = Lambda;
Lambda_s = Lambda_s - 0.05;
cuntUnL3 =1;

if Lambda_s>=1 && Lambda_s<=Lambda
    while Lambda_s>=1
    for i=1:21
        Lambda_Dmax(i)=max(1,eval(Lambda_maxdir(i)));
        nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
        nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
        P_pp_s(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nmin_s,n_max,Rbar_s,q_s)'...
            *NIntegral(nmin_s,n_max,Rbar_s,q_s,eval(Lambda_D(i)));
        P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
            *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
        P_s(:,:,q)=P_s(:,:,q)+((P_pp_s(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
        P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
    end
    T(:,:,q)=(exp(-exp(-E1/(R*temp))*af1*I^gamma*(t/T1)))*((P_s(:,:,q)-(P_s(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))))+((1-(exp(-exp(-E1/(R*temp))*af1*I^gamma*(t/T1)))))'...
        *((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))));
    if T(1,1,q)<=0
        break
    end
    Lambda_s    = Lambda_s - Step_size;
    q           = q + 1;
    cuntUnL3    = cuntUnL3 + 1;
    P_s(:,:,q) = 0;
    P_b(:,:,q) = 0;
    end
end

%%

%% Changing stress to MPa
for i=1:(q-1)
stress(i)=T(1,1,i)/10^6;
end
Micro_strain(1) = 0;
Micro_strain(2) = Step_size;
for i = 3:(cuntL1-1)
Micro_strain(i) = Micro_strain(i-1) + Step_size;
end
for i = cuntL1:(cuntL1+cuntUnL1-2)
Micro_strain(i) = Micro_strain(i-1) - Step_size;
end

for i = (cuntL1+cuntUnL1-1):(cuntL1+cuntUnL1+cuntL2-3)
Micro_strain(i) = Micro_strain(i-1) + Step_size;
end
for i = (cuntL1+cuntUnL1+cuntL2-2):(cuntL1+cuntUnL1+cuntL2+cuntUnL2-4)
Micro_strain(i) = Micro_strain(i-1) - Step_size;
end

for i = (cuntL1+cuntUnL1+cuntL2+cuntUnL2-3):(cuntL1+cuntUnL1+cuntL2+cuntUnL2+cuntL3-5)
Micro_strain(i) = Micro_strain(i-1) + Step_size;
end
for i = (cuntL1+cuntUnL1+cuntL2+cuntUnL2+cuntL3-4):(cuntL1+cuntUnL1+cuntL2+cuntUnL2+cuntL3+cuntUnL3-6)
Micro_strain(i) = Micro_strain(i-1) - Step_size;
end

%Micro_strain = 0:0.5:11;
Macro_strain = (1-(0.2^(1/3))).*(Micro_strain+1)+(0.2^(1/3)-1);
% Experimental data point
figure(1);
hold on
 num = xlsread('Virgin.xlsx');
 strain_exp=num(:,1);
 stress_exp=num(:,2);
or = [0.93 0.69 0.13];
cr = [1 0.667 0];
pu = [0.71 0.27 1];
gr = [ 0.39 0.83 0.07];
plot(Macro_strain,stress,'-','Color','k','LineWidth',4)
 hold on
 plot((strain_exp/100),stress_exp,'d','Color','r','LineWidth',3,'MarkerSize',10)
xlim([0 3.5])
ylim([0 10])
 title('\theta = 60')
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