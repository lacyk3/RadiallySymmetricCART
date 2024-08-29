%% CAR Sensitivity Analysis Figures
%clear all; close all; clc;
param = 'q'
load(strcat(param,'Sims.mat'))
%% Plot Tumor and CAR T cell trajectories
figure
subplot(2,1,1)
set(gca, 'FontSize', 18)
grid on
xlim([0 30])
ylim([0 15*10^8])
title(param)
ylabel('Tumor')
hold on
subplot(2,1,2)
set(gca, 'FontSize', 18)
ylabel('CAR T cells')
xlabel('Days post treatment')
xlim([0 30])
ylim([0 10*10^8])
grid on
hold on

C=jet(size(Results,1));
C = flipud(C);
labels = {};
for count = [1:4]
params = Results{count,1} ;
labels = [labels,num2str(params(8),2)];
t = Results{count,2};
U = Results{count,3};
V = Results{count,4};
Tdim = Results{count,5}; 
Cdim = Results{count,6}; 
[minval,min_idx] = min(max(Tdim,1000));
[maxval,max_idx] = max(Cdim);
subplot(2,1,1)
plot(t, Tdim, 'color',C(count,:), 'LineWidth',3)
plot(t(min_idx), Tdim(min_idx),'.', 'color',C(count,:), 'MarkerSize',20)
subplot(2,1,2)
plot(t, Cdim, 'color',C(count,:), 'LineWidth',3)
plot(t(max_idx), Cdim(max_idx),'k.','color',C(count,:), 'MarkerSize',20)
end

pos = [238   169   334   628]
set(gcf, 'Position', pos)

%% Tumor nadir/time against CAR T cell peak

figure
subplot(2,1,1)
set(gca, 'FontSize', 18)
hold on
grid on
xlabel('CAR T cell peak')
ylabel('tumor nadir')
subplot(2,1,2)
set(gca, 'FontSize', 18)
hold on
grid on
xlabel('CAR T cell peak')
ylabel('time tumor nadir')


for count = [1:4]
params = Results{count,1} ;
t = Results{count,2};
U = Results{count,3};
V = Results{count,4};
Tdim = Results{count,5}; 
Cdim = Results{count,6}; 
[minval,min_idx] = min(max(Tdim,1000));
[maxval,max_idx] = max(Cdim);
%AUC = trapz(t,Cdim);
subplot(2,1,1)
plot(maxval,minval,'.','MarkerSize',20,'color',C(count,:))
subplot(2,1,2)
plot(maxval,t(min_idx),'.','MarkerSize',20,'color',C(count,:))
end

subplot(2,1,1)
legend(labels)
pos = [238   169   334   628];
set(gcf, 'Position', pos)