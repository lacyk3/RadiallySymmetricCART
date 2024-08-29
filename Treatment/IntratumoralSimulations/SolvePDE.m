function [t, U,V, Tdim, Cdim] = SolvePDE(IC,params,tF)
% Dimensional parameters
R0 = 1;
DC = params(1); 
a = params(2);
b = params(3);
d = params(4); 
j = params(5); 
k = params(6);
m = params(7);
q = params(8);
DT_star = params(9);
u_star = params(10);

% Dimensionless parameters
alpha = d^2*j/(k*a*b^2);
beta_T = 1;
beta_C = (k*q*b)/(d^2*j);
gamma = d/a;
chi = (m*k*b^2)/(d^2*j);
zeta = d^2/(k*b^2);
DC = DC/a;
DT_star = DT_star/a;
u_star = u_star*b;
s = params(11); 
l = params(12);

% initializes discretization
dr = 0.015;
dt = min(dr*a/DC, dr/5);
T = tF*a; 
N = round(T/dt + 1);

% Tumor Initial condition
load(IC);
adjust_r = round(2.25*length(Rinit)/Rinit(end)); 
r = Rinit(1:adjust_r);
U = zeros(length(r), N);
V = U;
U(:,1) = Uinit(1:adjust_r);

InitC  = params(13)*b; 

width = .1;
r1 = round(length(r)/r(end)*width);
height = InitC/(sum((r(1:r1-1).^2).*exp(1 - 1./(1 - (r(1:r1 - 1)/r(r1)).^2)))*dr*4*pi);
V(1:r1 - 1,1) = exp(1 - 1./(1 - (r(1:r1 - 1)/r(r1)).^2))*height;

options = optimset('Display','off', 'FinDiffType', 'central', 'UseParallel', true);
count = 0
tic   
for n = 2:N
     t = n*dt  ;
     if min(real(U(:,n-1)))<0
         disp('broke early')
        break 
     end

     if  U(round(0.8*size(U,1)),n-1) > 0 %|| V(end,n-1) > 10^(-3) 
         disp('Domain extended')
        r = [0:dr:r(end)*1.25];
        U_new = zeros(length(r), N);
        V_new = U_new;
        U_new(1:length(U(:,n-1)),1:N) = U;
        V_new(1:length(V(:,n-1)),1:N) = V;
        U = U_new;
        V = V_new;        
     end
     
     if V(end,n-1) > 10^(-5) && count ==0
         disp('CAR T-cells at the boundary')
         count = count + 1;
     end
         BaneOfMyExistence = @(u_nn) TumorCARTFunc([U(:,n-1); V(:,n-1)], u_nn, dt, dr, DT_star, DC, beta_T, gamma, s,l, alpha, zeta, beta_C, chi, r,u_star);
         u_nn = fsolve(BaneOfMyExistence, [U(:,n-1); V(:,n-1)], options);
              
    U(:,n) = max(real(u_nn(1:length(u_nn)/2)),0);
    V(:,n) = max(real(u_nn(length(u_nn)/2 + 1:end)),0);
     
end 
toc

%%

C = []; T = []; 
for n = 1:length(U(1,:))
    C = [C sum(4*pi*V(:,n).*(r.').^2)*dr];
    T = [T sum(4*pi*U(:,n).*(r.').^2)*dr];
end
Tdim = T./b;
Cdim = C./b;
t = dt*(1:N)/a; % dimensional time for plotting
end