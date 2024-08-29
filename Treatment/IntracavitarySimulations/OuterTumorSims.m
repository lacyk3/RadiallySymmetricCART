%% Intracavitary injection Sims
clear all; close all; clc
%% Type 1: High proliferation, low diffusivity
doses = [1e7 5e7 1e8 5e8 1e9];
times = [84 84 84 84 84]; % duration of simulation 
SimResults = cell(length(doses),4);
count = 0;
diameter = 2; % detectable diameter at time of treatment in cm
width = 0.15; % width of CAR T cell layer for initial condition

for j = 1:numel(doses) 
d = doses(j);
T = times(j);
count = count+1 % counter to watch progress
[U,V,r,t] = RunTreatmentSimO(1,diameter,d, width, T);
SimResults{count,1} = U;
SimResults{count,2} = V;
SimResults{count,3} = r;
SimResults{count,4} = t;
end

save('OuterType1-2cm.mat','SimResults')
clear all
%% Type 2: Medium proliferation, medium diffusivity
doses = [1e7 5e7 1e8 5e8 1e9];
times = [84 84 84 84 84]; % duration of simulation 
SimResults = cell(length(doses),4);
count = 0;
diameter = 2;
width = 0.15;
for j = 1:numel(doses)
	d = doses(j);
	T = times(j);
count = count+1
[U,V,r,t] = RunTreatmentSimO(2,diameter,d, width, T);
SimResults{count,1} = U;
SimResults{count,2} = V;
SimResults{count,3} = r;
SimResults{count,4} = t;
%
end

save('OuterType2-2cm.mat','SimResults')
clear all
%% Type 3: Low proliferation, high diffusivity
doses = [1e7 5e7 1e8 5e8 1e9];
times = [84 84 84 84 84]; % duration of simulation 
SimResults = cell(length(doses),4);
count = 1;

diameter = 2;
width = 0.15;
for j = 1:numel(doses)
	d = doses(j);
	T = times(j);
[U,V,r,t] = RunTreatmentSimO(3,diameter,d, width, T);
SimResults{count,1} = U;
SimResults{count,2} = V;
SimResults{count,3} = r;
SimResults{count,4} = t;
count = count+1
end

save('OuterType3-2cm.mat','SimResults')
clear all
%% Type 4: Low proliferation, low diffusivity
doses = [1e7 5e7 1e8 5e8 1e9];
times = [84 84 84 84 84]; % duration of simulation 
SimResults = cell(length(doses),4);
count = 1

diameter = 2;
width = 0.15;
for j = 1:numel(doses)
	d = doses(j);
	T = times(j);
[U,V,r,t] = RunTreatmentSimO(4,diameter,d, width, T);
SimResults{count,1} = U;
SimResults{count,2} = V;
SimResults{count,3} = r;
SimResults{count,4} = t;
count = count+1
end

save('OuterType4-2cm.mat','SimResults')
clear all