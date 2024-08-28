%% Build heat maps of tumor growth for Figure S1 
clear all;  clc; close all;
%% Set up
% Load relevant simulations
a = '025';
load(strcat('a',a,'.mat'))
% Carrying capacity for re-dimensionalization
b = (4*pi/3)*1e-9; % cm^3/cell

%% Tumor Burden at detection 
figure
colormap(jet)%flipud(parula))
surf(U_star_grid,DT_star_grid,TBsatDetection_grid/b, 'EdgeColor', 'none', 'FaceAlpha',0.5)
hold on
title(strcat('a = 0.',a))
% -- for a = 0.25
%plot3([0.5, 0.01], [1e-5, 1e-4], [3e9, 3e9], 'ko','MarkerFaceColor','w', 'MarkerSize',10, 'LineWidth',2)
% -- for a = 0.125
%plot3([0.1], [2e-5], [3e9], 'ko','MarkerFaceColor','w', 'MarkerSize',10, 'LineWidth',2)
% -- for a = 0.025
plot3([0.01], [1e-4], [3e9], 'ko','MarkerFaceColor','w', 'MarkerSize',10, 'LineWidth',2)
set(gca, 'FontSize', 18)
grid off
xlim([0 .9])
ylim([5e-6 1.4e-4])
caxis([7*10^8 1.2*10^9])

colorbar
for val = [5 6 7 8 9 10 11 12 13 14 15]*10^8
contour(U_star_grid,DT_star_grid,TBsatDetection_grid/b,[val val], 'LineWidth',3)
end

set(gca, 'View',[-0.6000   90.0000])
set(gca,'TickDir','out');
ax = gca;
k = 0.02;
ax.TickLength = [k,k];
ax.LineWidth = 2;
pos = [900   424   275   275];
set(gcf, 'Position', pos)


%% tumor volumne doubling time 
figure
colormap(flipud(jet))%flipud(parula))
surf(U_star_grid,DT_star_grid,TVDTs_grid, 'EdgeColor', 'none', 'FaceAlpha',0.5)
hold on
title(strcat('a = 0.',a))

% -- for a = 0.25
%plot3([0.5, 0.01], [1e-5, 1e-4], [500, 500], 'ko','MarkerFaceColor','w', 'MarkerSize',10, 'LineWidth',2)
% -- for a = 0.125
%plot3([0.1], [2e-5], [500], 'ko','MarkerFaceColor','w', 'MarkerSize',10, 'LineWidth',2)
% -- for a = 0.025
plot3([0.01], [1e-4], [500], 'ko','MarkerFaceColor','w', 'MarkerSize',10, 'LineWidth',2)

set(gca, 'FontSize', 18)
grid off
xlim([0 .9])
ylim([5e-6 1.4e-4])
caxis([21 500])

%colorbar
for val = [10 20 30 40 50 60 70 80 90 100 150 200 250 300 350 400 450 500]
contour(U_star_grid,DT_star_grid,TVDTs_grid,[val val], 'LineWidth',3)
end

set(gca, 'View',[-0.6000   90.0000])
set(gca,'TickDir','out');
ax = gca;
k = 0.02;
ax.TickLength = [k,k];
ax.LineWidth = 2;
pos = [900   424   275   275];
set(gcf, 'Position', pos)

