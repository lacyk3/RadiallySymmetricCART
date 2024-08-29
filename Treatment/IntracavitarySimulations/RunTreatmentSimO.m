function [U,V,r,t] = RunTreatmentSimO(TumorType,diameter, dose,width, T)

% Shared dimensional parameters
R0 = 1;
DC = 0.0138;
b = (4*pi*R0^3/3)*1e-9; % includes spherical volume factor
d = .8*2.25; 
j = .8*0.6; 
k = 1.2*2.019e7; 
m = 1.2*0.293;
q = 1.2*1.53e-9; 

% -- Depending on TumorType, set parameters and load initial condition
if TumorType == 1 
    a = 0.25;
    DT_star = 10^(-5);
    u_star = 0.5/b;
    if diameter == 3
    load('Type1_15mm.mat');
    elseif diameter ==1
    load('Type1_5mm.mat');
    else 
    load('Type1_10mm.mat')
    end
elseif TumorType == 2 
    a = 0.125;
    DT_star =2* 10^(-5);
    u_star = 0.1/b;
    if diameter == 3
    load('Type2_15mm.mat');
    elseif diameter ==1
    load('Type2_5mm.mat');
    else 
    load('Type2_10mm.mat')
    end
elseif TumorType == 3
    a = 0.025;
    DT_star = 10^(-4);
    u_star = 0.01/b;
    if diameter == 3
    load('Type3_15mm.mat');
    elseif diameter ==1
    load('Type3_5mm.mat');
    else 
    load('Type3_10mm.mat')
    end
elseif TumorType == 4
    a = 0.25;
    DT_star = 10^(-4);
    u_star = 0.01/b;
    if diameter == 3
    load('Type4_15mm.mat');
    elseif diameter ==1
    load('Type4_5mm.mat');
    else 
    load('Type4_10mm.mat')
    end

end
% Nondimensional parameters
alpha = d^2*j/(k*a*b^2);
beta_T = 1;
beta_C = (k*q*b)/(d^2*j);
gamma = d/a;
chi = (m*k*b^2)/(d^2*j);
zeta = d^2/(k*b^2);
s = 1.2*0.36; 
l = 1.2*1.43;

DC = DC/a/R0^2;
DT_star = DT_star/a/R0^2;

u_star = u_star*b;

%%% initialize discretization
dr = 0.015;
dt = min(dr*a/DC, dr/5);
T = T*a;
N = round(T/dt + 1);

%%% set initial conditions
% -- use larger domain for larger initial tumor spread
if diameter < 2 
    if TumorType == 3 
        adjust_r = round(2.25*length(Rinit)/Rinit(end)); 
    else
    adjust_r = round(2*length(Rinit)/Rinit(end));
    end
else
    if TumorType == 3
        adjust_r = round(3*length(Rinit)/Rinit(end)); 
    else
        adjust_r = round(2.25*length(Rinit)/Rinit(end));
    end
end

r = Rinit(1:adjust_r);
U = zeros(length(r), N);
V = U;
U(:,1) = Uinit(1:adjust_r);
 
%--intracavitary
InitC  = dose*b; %nondimensionalize dose
temp = find(Uinit >0.03); % find edge of tumor
r2 = temp(end)-1+round(length(r)*width/r(end));
r1 = temp(end)-1;
height = InitC/(sum((r(r1:r2).^2))*dr*4*pi);
V(r1:r2,1) = height;

options = optimset('Display','off', 'FinDiffType', 'central', 'UseParallel', true);
tic    
for n = 2:N
        
     if min(real(U(:,n-1)))<0
         disp('broke early')
        break 
     end

     if  U(round(0.8*size(U,1)),n-1) > 0 %|| V(end,n-1) > 10^(-3) 
         disp('domain extended')
        r = [0:dr:r(end)*1.25];
        U_new = zeros(length(r), N);
        V_new = U_new;
        U_new(1:length(U(:,n-1)),1:N) = U;
        V_new(1:length(V(:,n-1)),1:N) = V;
        U = U_new;
        V = V_new;        
     end
         
         BaneOfMyExistence = @(u_nn)TumorCARTFunc([U(:,n-1); V(:,n-1)], u_nn, dt, dr, DT_star, DC, beta_T, gamma, s,l, alpha, zeta, beta_C, chi, r,u_star);
         u_nn = fsolve(BaneOfMyExistence, [U(:,n-1); V(:,n-1)], options);
              
    U(:,n) = max(real(u_nn(1:length(u_nn)/2)),0);
    V(:,n) = max(real(u_nn(length(u_nn)/2 + 1:end)),0);
     
end 
toc

t = dt*(1:N)/a;

end
