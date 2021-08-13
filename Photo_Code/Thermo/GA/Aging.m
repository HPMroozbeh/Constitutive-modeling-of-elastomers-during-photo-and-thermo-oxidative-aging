function u = Aging(xx,t,Des_temp,Stretch)
format long
%%
%reference temperature
temp = 333;
% desired temperature
% decay function parameters
T1 = xx(1);
T2 = xx(2);
%Constants
E1 = xx(3);
R = 8.314;
I = 1;
alpha = 1;
beta = 1;
gamma = 1;
%Time-temperature superposition
af1 = exp((E1/R)*((1/temp)-(1/Des_temp)));
% Network peobability fitted parameters
% reletive chain length
Rbar_s = 4;
R_1 = xx(4);
R_2 = xx(5);
R_I = 0;
Rbar_b = (Rbar_s) - (R_1 - (R_I/2))*(1-exp(-I^alpha*(t/T2)))+ (R_2 + (R_I/2))*(1-exp(-I^beta*(t/T2)));
% an indicator for crosslink density
q_s = 1.0;
q_1 = xx(6);
q_2 = xx(7);
q_I = 0;
q_b = (q_s) - (q_1 - (q_I/2))*(1-exp(-I^alpha*(t/T2)))+ (q_2 + (q_I/2))*(1-exp(-I^beta*(t/T2)));
if q_b > 1
    q_b =1;
end

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
phi= 1/3;
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
%% First Loading

nos = 10;
Chi = linspace(Stretch(1,1),(length(Stretch)-1)/4+1,nos);
lambda = (Chi - C^(phi))./(1 - C^(phi));


for j= 1:nos
    Lambda_s = lambda(j);
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
     q=q+1;
     cunt = cunt+1;
     P_s(:,:,q)=0;
     P_b(:,:,q)=0;
end

% Changing stress to MPa
for i=1:(q-1)
stress(i)=T(1,1,i)/10^6;
end
u(1,:) = stress;
end