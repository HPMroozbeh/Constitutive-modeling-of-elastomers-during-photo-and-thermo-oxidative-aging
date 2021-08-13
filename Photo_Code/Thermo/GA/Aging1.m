function u = Aging(xx,t,Des_temp,Stretch)
format long
% aging time
% decay function parameters
T1 = xx(1);
T2 = T1 + xx(2);
%Constants
E1 = xx(3);
E2 = E1 + xx(4);
R = 8.314;
I = 1;
alpha = 1;
beta = 1;
%reference temperature
temp = 333;
%Time-temperature superposition
af1 = exp((E1/R)*((1/temp)-(1/Des_temp)));
af2 = exp((E2/R)*((1/temp)-(1/Des_temp)));
% Network peobability fitted parameters
% reletive chain length
Rbar_s = 4;
R_2 =  xx(5);
R_1 = R_2 + xx(6);
Rbar_b = (Rbar_s) - R_1*(1-exp(-exp(-E1/(R*temp))*I^alpha*af1*(t/T1)))+ R_2*(1-exp(-exp(-E2/(R*temp))*I^beta*af2*(t/T2)));
% an indicator for crosslink density
q_s = 1.0;
q_2 = xx(7);
q_1 = q_2 + xx(8);
q_b = (q_s) - q_1*(1-exp(-exp(-E1/(R*temp))*I^alpha*af1*(t/T1)))+ q_2*(1-exp(-exp(-E2/(R*temp))*I^beta*af2*(t/T2)));
nue = 1.02;
n_max = 100;
N_s00= 19.04*10^6;
% calculating the area under probability functions for network s & b
%area_bT = ProbIntegral(nue*Rbar_bT,n_max,Rbar_bT,q_bT);
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

q=1;
T=zeros([3 3 1]);
P_b=zeros([3 3 1]);
P_pp_b=zeros(1,21);
cunt = 1;
nos = 20;
Chi = linspace(1,(length(Stretch)-1)/4+1,nos);
lambda = (Chi - C^(1/3))./(1-C^(1/3));

for j= 1:nos
    Lambda_s = lambda(j);
    for i=1:21
        Lambda_Dmax(i)=max(1,eval(Lambda_D(i)));
        nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
     P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
            *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
        P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
    end
    T(:,:,q)=((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))));
     q=q+1;
     cunt = cunt+1;
     P_b(:,:,q)=0;
end


% Changing stress to MPa
for i=1:(q-1)
stress(i)=T(1,1,i)/10^6;
end
u(1,:) = stress;
end