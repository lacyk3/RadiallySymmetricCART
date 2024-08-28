%% Build Tumor Growth Figure
clc; clear all; %close all;
%% Show tumor growth for 4 tumor types

% -- Type 1: HPLD ~165 day TVDT

%Dimensional parameters
% a = 0.25;
% b = (4*pi/3)*1e-9;
% u_star = .5;
% u_detection = 10^8*b;
% DT_star = 10^(-5);
% dr = 0.015;
% % - tumor burden at 1 cm, 1.5 cm, 2 cm detectable tumor radius
% load('Type1_20mm.mat')
% T20mm = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b;
% load('Type1_15mm.mat')
% T15mm = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b;
% load('Type1_10mm.mat')
% T10mm = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b;


% % -- Type 2: HPHD ~38 day TVDT
a = 0.125;
b = (4*pi/3)*1e-9;
u_star = .1;
u_detection = 10^8*b;
DT_star = 2*10^(-5);
dr = 0.015;
load('Type2_20mm.mat')
T20mm = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b;
load('Type2_15mm.mat')
T15mm = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b;
load('Type2_10mm.mat')
T10mm = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b;
% 
% % % -- Type 3: LPHD ~ 182 day TVDT
% a = 0.025;
% b = (4*pi/3)*1e-9;
% u_star = 0.01;
% u_detection = 10^8*b;
% DT_star = 10^(-4);
% dr = 0.015;
% load('Type3_20mm.mat')
% T20mm = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b;
% load('Type3_15mm.mat')
% T15mm = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b;
% load('Type3_10mm.mat')
% T10mm = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b;
% 
% % -- Type 4: LPLD ~ 463 day TVDT
% a = 0.25;
% b = (4*pi/3)*1e-9;
% u_star = 0.01;
% u_detection = 10^8*b;
% DT_star = 10^(-4);
% dr = 0.015;
% load('Type4_20mm.mat')
% T20mm = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b;
% load('Type4_15mm.mat')
% T15mm = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b;
% load('Type4_10mm.mat')
% T10mm = sum(4*pi*Uinit.*(Rinit.').^2)*dr/b;

% Dimensionless parameters
beta_T = 1;
DT_starnd = DT_star/a;

%%% Plot Initial condition
figure
hold on
%%%-- Initial 1 mm tumor mass
dr = 0.015;
r = [0:dr:1];
r1 = round(0.1*length(r));
U = zeros(length(r), 1);
U(1:r1 - 1,1) = exp(1 - 1./(1 - (r(1:r1 - 1)/r(r1)).^4))/beta_T;

plot([-flip(r) r].',[flip(U(:,1)); U(:,1)], 'LineWidth',1.3, 'color', [.1 .1 .1 .4])
%%%-- 1 cm tumor 
plot([-flip(Rinit) Rinit].',[flip(Uinit); Uinit], 'LineWidth',1.3,'color','k')
plot([-1 1], [u_detection u_detection],'k:', 'LineWidth',1.5)
xlim([-4.25 4.25])
yticks([0.4189 0.8378])
ylim([0 1.25])
set(gca,'FontSize',4)
set(gca,'xticklabels',[])
set(gca,'yticklabels',[])
set(gca,'LineWidth',1, 'TickLength',[0.02, 0.005])
x_width=214; y_width=86;
set(gcf, 'Position', [100 100 x_width y_width]); %

%% Simulate untreated tumor growth
% - time to progression by RECIST criteria
% - tumor volume doubling time
dr = 0.015;
dt = a/10;
T = 200*a;
N = round(T/dt + 1);
r = Rinit;
U = zeros(length(r), N);
U(:,1) = Uinit;

options = optimset('Display','off');%,  'UseParallel', true,'FinDiffType', 'central');%
flag = 0;
flag1 = 0;
for n = 2:N
    detectable = U(:,n-1)>u_detection;
    tumorsize = r(detectable);
    if length(tumorsize)>1
        if flag == 0;
            if tumorsize(end)> 2^(1/3)
                TVDT = dt*n/a
                flag = 1;
            end
        end
        if flag1 == 0;
            if tumorsize(end)> 1.2
                TTP = dt*n/a
                flag1 = 1;
            end
        end
    end
    if min(real(U(:,n-1)))<0
        break 
    end
    % If tumor has reached edge domain, extend the domain 
    if  U(end,n-1) > 0
        r = [0:dr:r(end)*1.1];
        U_new = zeros(length(r), N);
        U_new(1:length(U(:,n-1)),1:N) = U;
        U = U_new;        
    end
    % Solve PDE (u_n, u_nn, dt, dr, DT, beta_T, r, u_star)
    BaneOfMyExistence = @(u_nn)TumorFunc(U(:,n-1), u_nn, dt, dr, DT_starnd,beta_T, r, u_star);
    u_nn = fsolve(BaneOfMyExistence, U(:,n-1), options); 
    % Enforce that all density is positive
    U(:,n) = max(u_nn,0);
  
end 
%%
%al = linspace(0.05,1,28); [0.95,.425,0,al(count)]
c=parula(int16(N));
count = 0;
figure
for jj = 1:75:N
    count = count+1;
    hold on
    plot([-flip(r) r].',[flip(U(:,jj)); U(:,jj)], 'LineWidth',1.1,'Color',c(jj,:))
    axis([-4,4,0,1.25])
end

xlim([-4.25 4.25])
set(gca,'FontSize',4)
yticks([0.4189 0.8378])
set(gca,'xticklabels',[])
set(gca,'yticklabels',[])
set(gca,'LineWidth',1, 'TickLength',[0.02, 0.005])
x_width=214; y_width=86;
set(gcf, 'Position', [100 100 x_width y_width]); %


%% 

