%% Intratumoral injection simulations

clear all; close all; clc
%% Type 1: High proliferation, low diffusivity
diameter = 2; % detectable diameter at time of treatment in cm
doses = [1e7 5e7 1e8 5e8 1e9];
times = [84 84 84 84 84]; % duration for each dose
SimResults = cell(length(doses),4);
 
for j = 1:length(doses)
 d = doses(j) % loop over dose levels in doses
 T = times(j); % loop over durations levels in times
 [U,V,r,t] = RunTreatmentSim(1,diameter,d,T);
 SimResults(j,1) = {U}; % save tumor density
 SimResults(j,2) = {V}; % save CART density
 SimResults(j,3) = {r}; % save spatial domain
 SimResults(j,4) = {t}; % save time points
 end
 
save('InnerType1-2cm.mat','SimResults')
clear all
%% Type 2: medium proliferation, medium diffusivity
diameter = 2; % detectable diameter at time of treatment in cm
doses = [1e7 5e7 1e8 5e8 1e9];
times = [84 84 84 84 84]; % duration for each dose
SimResults = cell(length(doses),4);

for j = 1:length(doses)
d = doses(j)
T = times(j);
[U,V,r,t] = RunTreatmentSim(2,diameter,d,T);
SimResults(j,1) = {U};
SimResults(j,2) = {V};
SimResults(j,3) = {r};
SimResults(j,4) = {t};
end

save('InnerType2-2cm.mat','SimResults')
clear all
%% Type 3: low proliferation, high diffusivity
diameter = 2; % detectable diameter at time of treatment in cm
doses = [1e7 5e7 1e8 5e8 1e9];
times =[84 84 100 100 100];
SimResults = cell(length(doses),4);
 
for j = 1:length(doses)
d = doses(j)
T=times(j);
[U,V,r,t] = RunTreatmentSim(3,diameter,d,T);
SimResults(j,1) = {U};
SimResults(j,2) = {V};
SimResults(j,3) = {r};
SimResults(j,4) = {t};
end
 
save('InnerType3-2cm.mat','SimResults')
clear all
%% Type 4: high proliferation, high diffusivity
diameter = 2; % detectable diameter at time of treatment in cm
doses = [1e7 5e7 1e8 5e8 1e9];
times =[70 70 70 84 84];
SimResults = cell(length(doses),4);

for j = 1:length(doses)
d = doses(j)
T = times(j);
[U,V,r,t] = RunTreatmentSim(4,diameter,d,T);
SimResults(j,1) = {U};
SimResults(j,2) = {V};
SimResults(j,3) = {r};
SimResults(j,4) = {t};
end

save('InnerType4-2cm.mat','SimResults')
clear all

