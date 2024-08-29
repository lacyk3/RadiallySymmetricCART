%% Visualize simulation outputs to evaluate dose outcomes
clear all; close all; clc
%%
load(['InnerType1-2cm.mat'])
type = 'Inner Type 1, diameter = 2cm, a = .25, DT* = 10^{-5}';

%% Parameters in order to re-dimensionalize and include the exhausted CARTs
% -- Type 1 and 2
a = 0.25;
% -- Type 3
%a = 0.025;
% -- Type 4
%a = 0.125;
R0 = 1;
DC = 0.0138;%.0128; 
b = (4*pi*R0^3/3)*1e-9;
d = 0.8*2.25; 
j = 0.8*0.6; 
k = 1.2*2.019e7;
u_detection = 10^8*b;

alpha = d^2*j/(k*a*b^2);
q_dim = 1.2*1.53*1e-8; % dimensional exhaustion rate
q = alpha*(k*q_dim*b)/(d^2*j); % nondimensional exhaustion rate
m_dim = 0.12;  % dimensional death rate
m = alpha*(m_dim*k*b^2)/(d^2*j); %nondimensional death rate
dr = 0.015;
%% -- Plot total tumor burden, detectable diameter, and total CAR Ts over time

colors = lines(10);
figure
for idx = 1: size(SimResults,1)
U = SimResults{idx,1};
V = SimResults{idx,2};
r = SimResults{idx,3};
t = SimResults{idx,4};
dt = t(2)-t(1);
t = (t-dt)*a;
dt = dt*a;
C0 = sum(4*pi*V(:,1).*(r.').^2)*dr;
E0 =0*C0/50;
C = []; T = [];RT = []; E = []; W = zeros(size(V));
for n = 1:length(U(1,:))
    C = [C sum(4*pi*V(:,n).*(r.').^2)*dr]; % total CAR T cells
    T = [T sum(4*pi*U(:,n).*(r.').^2)*dr];  % total tumor cells
    temp = find(U(:,n)>u_detection); % where above limit of detection?
    if temp
        RT = [RT r(temp(end))]; % track tumor radius
    else 
        RT = [RT 0];
    end
    % track exhausted cells generated
    if n == 1
        Ex = zeros(length(U(:,1)),1);
    else
        Ex = exp(-m.*t(n)).*sum(q*V(:,1:n).*U(:,1:n).*exp(m.*t(1:n)).*dt,2); 
    end
    W(:,n) = Ex;
    E = [E sum(4*pi*Ex.*(r.').^2)*dr+E0*exp(-m.*t(n))]; % sum up the total exhausted cells
end
t = SimResults{idx,4};
Diameter_init = 2*RT(1);

subplot(3,1,1)
hold on
plot(t,T(1:end)/b, 'LineWidth',2,'color',colors(idx,:), 'LineStyle', '-', 'MarkerSize',1)
subplot(3,1,2)
hold on
plot(t, 2*RT/Diameter_init,'LineWidth',2,'color',colors(idx,:), 'LineStyle', '-', 'MarkerSize',1)
subplot(3,1,3)
hold on
plot(t,(C(1:end)+E(1:end))/b,'LineWidth',2,'color',colors(idx,:), 'LineStyle', '-', 'MarkerSize',1)
end

%% Add reference lines for evaluation criteria
subplot(3,1,1)
hold on
plot([0 70], [T(1) T(1)]/b,'k--') % initial tumor burden
plot([28 28], [0 2.5*T(1)/b],'k')
plot([42 42], [0 2.5*T(1)/b],'k')
plot([56 56], [0 2.5*T(1)/b],'k')
ylabel('Tumor Cells')
title(type)
xlim([0,70])
ylim([0 2.5*10^9])%

subplot(3,1,2)
hold on
ylabel('Tumor Diameter')
plot([0 70], [1.2 1.2],'k--') % progression
plot([0 70], [0.7 0.7],'k--') % partial response
plot([28 28], [0 2],'k') % week 4
plot([42 42], [0 2],'k') % week 6
plot([56 56], [0 2],'k') % week 8

xlim([0,70])
ylim([0 2])

subplot(3,1,3)
ylabel('CAR T-Cells')
xlabel('days post treatment')
xlim([0,70])
