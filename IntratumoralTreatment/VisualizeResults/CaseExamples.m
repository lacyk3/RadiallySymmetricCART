%% Treatment Outcome Examples
clear all;  clc; close all
%%

% -- Load tumor outcome without treatment
load('Type4_untreated.mat')

% -- Calculate totals over time
dr = 0.015;
a = 0.25;
b = (4*pi/3)*1e-9;
u_detection = 10^8*b;


 T = [];RT = []; 
for n = 1:length(U(1,:))
    T = [T sum(4*pi*U(:,n).*(r.').^2)*dr]; 
    temp = find(U(:,n)>u_detection);
    if temp
        RT = [RT r(temp(end))];
    else 
        RT = [RT 0];
    end
end
Diameter_init = 2*RT(1);

% -- plot
fig1 = figure;
left_color = [0 0 0];
right_color = 1.1*[0    0.6    0.75];
set(fig1,'defaultAxesColorOrder',[left_color; right_color]);
hold on
xlim([0 56])
ylim([0 2.5e9])

% -- Load tumor outcome with treatment
load(['OuterType4-2cm.mat']);
idx = 2; % pick dose of interest
dr = 0.015;
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

U = SimResults{idx,1};
V = SimResults{idx,2};
r = SimResults{idx,3};
t = SimResults{idx,4};
dt = t(2)-t(1);
t = (t-dt)*a;
dt = dt*a;
C0 = sum(4*pi*V(:,1).*(r.').^2)*dr;
E0 = 0; % no exhausted cells in initial dose
C = []; T = [];RT = []; E = []; W = zeros(size(V));
detectableindex = [];
for n = 1:length(U(1,:))
    % effector CAR T cells
    C = [C sum(4*pi*V(:,n).*(r.').^2)*dr];
    % tumor cells
    T = [T sum(4*pi*U(:,n).*(r.').^2)*dr];
    % detectable tumor radius
    temp = find(U(:,n)>u_detection);
    if temp
        RT = [RT r(temp(end))];
        detectableindex = [detectableindex temp(end)];
    else 
        RT = [RT 0];
        detectableindex = [detectableindex 1];
    end
    % exhausted CAR T cells
    if n == 1
        Ex = zeros(length(U(:,1)),1);
    else
        Ex = exp(-m.*t(n)).*sum(q*V(:,1:n).*U(:,1:n).*exp(m.*t(1:n)).*dt,2);
    end
    W(:,n) = Ex;
    E = [E sum(4*pi*Ex.*(r.').^2)*dr+E0*exp(-m.*t(n))];
end
Diameter_init = 2*RT(1);
t = SimResults{idx,4};

%%
plot(t,T(1:end)/b,'k', 'LineWidth',2, 'LineStyle', '-', 'MarkerSize',1, 'color', left_color)

yyaxis right
plot(t,(C(1:end)+E(1:end))/b,'LineWidth',2, 'LineStyle', '-', 'MarkerSize',1, 'color', right_color)
hold on
ylim([0 2.5e9])
x_width=1000; y_width=345;
set(gcf, 'Position', [250 250 x_width/2 y_width/1.3]); %
set(gca, 'FontSize',14)
xlabel('time since treatment (days)')
%%
figure
%ylabel('Detectable tumor diameter')
hold on
plot(t, 2*RT/Diameter_init,'k','LineWidth',2, 'LineStyle', '-', 'MarkerSize',1)
plot([0 60], [1.2 1.2],'k-') % progression
plot([0 60], [0.7 0.7],'k-') % partial response
x_width=1000/3; y_width=345/3;
axis([0 56 0 1.5])
xticks([0 25 50])
yticks([0 1])
set(gcf, 'Position', [250 250 x_width/2 y_width/1.3]); %
set(gca, 'FontSize',14)

%%
fig1 = figure;
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);
left_color = [0 0 0];
right_color = [.01    0.8    0.5];
set(fig1,'defaultAxesColorOrder',[left_color; right_color]);
dT = 950;
js = [1,100,250,600,2500, 5000];
for j = 1:6
    jj = js(j);
    di = detectableindex(jj);
    Ui = U(1:di,jj);
subplot(2,3,j)
area([-flip(r(1:di)) r(1:di)].',[flip(Ui); Ui],'EdgeColor','none','FaceColor','k','FaceAlpha',.1)
hold on
yyaxis left
plot([-flip(r) r].',[flip(U(:,jj)); U(:,jj)], 'LineWidth',1.5)
plot([-2.5 2.5],[0.4189  .4189],'--', 'LineWidth',1.5);

if j == 5 
    xlabel('radial position')
end

if j == 1 || j == 2 || j == 3
    set(gca,'xticklabels',[])
end

if j == 2 || j == 3 || j == 5 || j == 6
    set(gca,'yticklabels',[])
end

axis([-3 3 0 1.2])
yyaxis right

plot([-flip(r) r].',[flip(V(:,jj)); V(:,jj)],'LineWidth',1.5, 'color',right_color)
set(gca, 'FontSize',14)
if j == 1 || j == 2 || j == 4 || j == 5
    set(gca,'yticklabels',[])
end
axis([-2.5 2.5 0 1])

title([num2str(round(t(jj))) ' days'],'FontWeight','Normal', 'FontSize',14)
set(gca,'FontSize',14)
end

x_width=1000; y_width=345;
set(gcf, 'Position', [250 250 x_width/2 y_width/1.3]); %
subplot(2,3,1)
set(gca, 'Position', [.1, 0.625 , 0.24, 0.25]);
subplot(2,3,2)
set(gca, 'Position', [.38, 0.625 , 0.24, 0.25]);
subplot(2,3,3)
set(gca, 'Position', [.66, 0.625 , 0.24, 0.25]);
subplot(2,3,4)
set(gca, 'Position', [.1, 0.225 , 0.24, 0.25]);
subplot(2,3,5)
set(gca, 'Position', [.38, 0.225 , 0.24, 0.25]);
subplot(2,3,6)
set(gca, 'Position', [.66, 0.225 , 0.24, 0.25]);
