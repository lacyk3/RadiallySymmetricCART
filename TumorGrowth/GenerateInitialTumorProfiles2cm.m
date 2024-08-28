%% Generate Initial Condition for 4 tumor types
% Start with 1 mm tumor mass at center of domain. Allow tumor to grow until
% detectable radius reaches 20 mm. Save tumor profile at detectable radius 
% of 5 mm, 10 mm, 15 mm, and 20 mm.

clear all; close all; clc

%% -- Set tumor parameters
type = 'Type2'; % 'Type1', 'Type3', 'Type4'
name1 = strcat(type,'_5mm.mat');
name2 = strcat(type,'_10mm.mat');
name3 = strcat(type,'_15mm.mat');
name4 = strcat(type,'_20mm.mat');

% Dimensional parameters
a = 0.125; % day-1 % 0.25, 0.025, 0.25
b = (4*pi/3)*1e-9; % cm^3/cell
u_star = 0.1/b; %.1/b cell/cm^3 % 0.5, 0.01, 0.01
DT_star = 2*10^(-5);% cm^2/day %10^(-5), 10^(-4), 10^(-4)

% Dimensionless parameters
beta_T = 1;
DT_starnd = DT_star/a; 
u_star = u_star*b;
u_detection = 10^8*b; 

% set up domain
dr = 0.015;
dt = a/10;
T = 1500*a; %large maximum time
N = round(T/dt + 1);

r = [0:dr:1];
r1 = round(0.25*length(r));
U = zeros(length(r), N);
U(1:r1 - 1,1) = exp(1 - 1./(1 - (r(1:r1 - 1)/r(r1)).^4))/beta_T;

options = optimset('Display','off', 'FinDiffType', 'central', 'UseParallel', true);
tic
for n = 2:N
    detectable = U(:,n-1)>u_detection;
    tumorsize = r(detectable);
    if length(tumorsize)>1
        if tumorsize(end)> 0.5
            break 
        end
    end
    % Break if tumor density has gone negative
    if min(real(U(:,n-1)))<0
        break 
    end
    % If tumor has reached edge domain, extend the domain 
    if  U(round(0.85*size(U,1)),n-1) > 0
        r = [0:dr:r(end)*1.5];
        U_new = zeros(length(r), N);
        U_new(1:length(U(:,n-1)),1:N) = U;
        U = U_new;        
    end
    % Solve PDE
    BaneOfMyExistence = @(u_nn)TumorFunc(U(:,n-1), u_nn, dt, dr, DT_starnd,beta_T, r, u_star);
    u_nn = fsolve(BaneOfMyExistence, U(:,n-1), options); 
    % Enforce that all density is positive
    U(:,n) = max(u_nn,0);
  
end 
dt*n/a % print time when detected radius is 5 mm
toc

% save
Rinit = [0:dr:4.5];
Uinit = zeros(length(Rinit),1);
Uinit(1:length(r),1)= U(:,n-1);
total_cells = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b % print total cell count 
figure
hold on
plot([-flip(r) r].',[flip(U(:,1)); U(:,1)], 'LineWidth',1.5, 'color', [.1 .1 .1 .5])
plot([-flip(r) r].',[flip(U(:,n-1)); U(:,n-1)], 'LineWidth',1.5,'color','k')
plot([-1 1], [u_detection u_detection],'k--')
set(gca, 'FontSize', 20)
save(name1, 'Uinit', 'Rinit')

%% Continue to next size
U2 = zeros(length(r), N);
U2(:,1) = U(:,n-1);
clear U

tic
for n = 2:N
    detectable = U2(:,n-1)>u_detection;
    tumorsize = r(detectable);
    if length(tumorsize)>1
        if tumorsize(end)> 1
            break 
        end
    end
    % Break if tumor density has gone negative
    if min(real(U2(:,n-1)))<0
        break 
    end
    % If tumor has reached edge domain, extend the domain 
    if  U2(round(0.85*size(U2,1)),n-1) > 0
        r = [0:dr:r(end)*1.5];
        U_new = zeros(length(r), N);
        U_new(1:length(U2(:,n-1)),1:N) = U2;
        U2 = U_new;        
    end
    % Solve PDE
    BaneOfMyExistence = @(u_nn)TumorFunc(U2(:,n-1), u_nn, dt, dr, DT_starnd,beta_T, r, u_star);
    u_nn = fsolve(BaneOfMyExistence, U2(:,n-1), options); 
    % Enforce that all density is positive
    U2(:,n) = max(u_nn,0);
  
end 
dt*n/a
toc

Rinit = [0:dr:4.5];
Uinit = zeros(length(Rinit),1);
Uinit(1:length(r),1)= U2(:,n-1);
total_cells = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b

plot([-flip(r) r].',[flip(U2(:,1)); U2(:,1)], 'LineWidth',1.5, 'color', [.1 .1 .1 .5])
plot([-flip(r) r].',[flip(U2(:,n-1)); U2(:,n-1)], 'LineWidth',1.5,'color','k')
plot([-1 1], [u_detection u_detection],'k--')
save(name2, 'Uinit', 'Rinit')

%% Continue to next size
U3 = zeros(length(r), N);
U3(:,1) = U2(:,n-1);
clear U2

tic
for n = 2:N
    detectable = U3(:,n-1)>u_detection;
    tumorsize = r(detectable);
    if length(tumorsize)>1
        if tumorsize(end)> 1.5
            break 
        end
    end
    % Break if tumor density has gone negative
    if min(real(U3(:,n-1)))<0
        break 
    end
    % If tumor has reached edge domain, extend the domain 
    if  U3(round(0.85*size(U3,1)),n-1) > 0
        r = [0:dr:r(end)*1.5];
        U_new = zeros(length(r), N);
        U_new(1:length(U3(:,n-1)),1:N) = U3;
        U3 = U_new;        
    end
    % Solve PDE
    BaneOfMyExistence = @(u_nn)TumorFunc(U3(:,n-1), u_nn, dt, dr, DT_starnd,beta_T, r, u_star);
    u_nn = fsolve(BaneOfMyExistence, U3(:,n-1), options); 
    % Enforce that all density is positive
    U3(:,n) = max(u_nn,0);
  
end 
dt*n/a
toc

Rinit = [0:dr:5];
Uinit = zeros(length(Rinit),1);
Uinit(1:length(r),1)= U3(:,n-1);
total_cells = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b

plot([-flip(r) r].',[flip(U3(:,1)); U3(:,1)], 'LineWidth',1.5, 'color', [.1 .1 .1 .5])
plot([-flip(r) r].',[flip(U3(:,n-1)); U3(:,n-1)], 'LineWidth',1.5,'color','k')
plot([-1 1], [u_detection u_detection],'k--')

save(name3, 'Uinit', 'Rinit')

%% Continue to final size

U4 = zeros(length(r), N);
U4(:,1) = U3(:,n-1);
clear U3

tic
for n = 2:N
    detectable = U4(:,n-1)>u_detection;
    tumorsize = r(detectable);
    if length(tumorsize)>1
        if tumorsize(end)> 2
            break 
        end
    end
    % Break if tumor density has gone negative
    if min(real(U4(:,n-1)))<0
        break 
    end
    % If tumor has reached edge domain, extend the domain 
    if  U4(round(0.85*size(U4,1)),n-1) > 0
        r = [0:dr:r(end)*1.5];
        U_new = zeros(length(r), N);
        U_new(1:length(U4(:,n-1)),1:N) = U4;
        U4 = U_new;        
    end
    % Solve PDE
    BaneOfMyExistence = @(u_nn)TumorFunc(U4(:,n-1), u_nn, dt, dr, DT_starnd,beta_T, r, u_star);
    u_nn = fsolve(BaneOfMyExistence, U4(:,n-1), options); 
    % Enforce that all density is positive
    U4(:,n) = max(u_nn,0);
  
end 
dt*n/a
toc

Rinit = [0:dr:5];
Uinit = zeros(length(Rinit),1);
Uinit(1:length(r),1)= U4(:,n-1);
total_cells = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b

plot([-flip(r) r].',[flip(U4(:,1)); U4(:,1)], 'LineWidth',1.5, 'color', [.1 .1 .1 .5])
plot([-flip(r) r].',[flip(U4(:,n-1)); U4(:,n-1)], 'LineWidth',1.5,'color','k')
plot([-1 1], [u_detection u_detection],'k--')

save(name4, 'Uinit', 'Rinit')
