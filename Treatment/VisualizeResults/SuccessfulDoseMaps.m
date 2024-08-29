% Generate visualization of responses based on number of effector cells
clear all; close all; clc
%% 
total_doses = 10^8*(0:.01:10);
percent_effector = .1*(0:.01:10);

[X,Y] = meshgrid(total_doses, percent_effector);

effective_doses = X.*Y;
response_map = effective_doses;

% %%-- Inner T1 
% RECIST
% SD_lowlim = 0*10^7;
% PR_lowlim =4.3*10^8;%3
% CR_lowlim = 4.5*10^8;
% type = 'Inner Type I';
% choi
% SD_lowlim = 0*10^7;
% PR_lowlim =2.5*10^8;%3
% CR_lowlim = 4.5*10^8;
% type = 'Inner Type I';

% Inner T2
% RECIST
SD_lowlim = 2.2e8;
PR_lowlim = 2.2e8 ;
CR_lowlim = 2.5e8;
type = 'Inner Type II';
%choi
% SD_lowlim = 1.1e8;
% PR_lowlim = 1.1e8 ;
% CR_lowlim = 5e8;
% type = 'Inner Type II';

% Inner T3
% RECIST
% SD_lowlim = 7e7;
% PR_lowlim = 1.3e8;
% CR_lowlim = 1.3*10^8;
% type = 'Inner Type III';
%choi
% SD_lowlim = 1e8;% 1e8;% 
% PR_lowlim = 1e8;%1*10^8;
% CR_lowlim = 1.3*10^8;
% type = 'Inner Type III';

% % Inner T4
% RECIST
% SD_lowlim = 2.5e8;%2.3e8
% PR_lowlim = 2.5e8;
% CR_lowlim = 2.5e8;
% type = 'Inner Type IV';
%choi
% SD_lowlim = 2.4e8;%2.3e8
% PR_lowlim = 2.4e8;
% CR_lowlim = 2.5e8;
% type = 'Inner Type IV';

% Outer T1
% SD_lowlim = 0*10^7;
% PR_lowlim = 8.8e8;%-----%9.5*10^8;%9.8
% CR_lowlim = 10*10^8;
% type = 'Outer Type I';

% % Outer T2
% SD_lowlim = 2e8;
% PR_lowlim = 8.5e8;%-----%9.6e8;%9.5
% CR_lowlim = 10*10^8;
% type = 'Outer Type II';

% Outer T3
% SD_lowlim = 3*10^8;%8e7
% PR_lowlim = 10*10^8;
% CR_lowlim = 10*10^8;
% type = 'Outer Type III';

% Outer T4
% SD_lowlim = 9.6e8; %2.5*10^7;
% PR_lowlim = 10*10^8;
% CR_lowlim = 10*10^8;
% type = 'Outer Type IV';


for j = 1:length(percent_effector)
    for jj = 1:length(total_doses)
        if effective_doses(j,jj) < SD_lowlim
            response_map(j,jj) = 4 ;
        elseif effective_doses(j,jj) < PR_lowlim
            response_map(j,jj) = 3;
        elseif effective_doses(j,jj) <= CR_lowlim 
            response_map(j,jj) = 2;
        else
            response_map(j,jj) = 1;
        end
    end
end

figure
cmap = colormap;
N = 4;
x = round(linspace(1, size(cmap, 1), N));
cmap = cmap(x, :);
colormap(cmap)
image(total_doses, flipud(percent_effector),response_map)
set(gca,'YDir','normal') 
grid on
set(gca, 'FontSize',15)
title(type)
xlabel([])
ylabel([])
yticks([.1 .2 .3 .4 .5 .6 .7 .8 .9])
yticklabels([])
xticks(1e9*[.1 .2 .3 .4 .5 .6 .7 .8 .9])
xticklabels([])
pos = [238   424   316   306];
set(gcf, 'Position', pos)