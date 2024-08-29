%% Sensitivity Analysis on CAR parameters
clear all; close all; clc

%% Loop through parameters

% 2 -- l
clear all
DC = 0.0138;
b = (4*pi/3)*1e-9;
d = .8*2.25; 
j = .8*0.6; 
k = 1.2*2.019e7;
m = 1.2*0.293;
q = 1.2*1.53e-9;
a = 0.25;
DT_star = 10^(-5);
u_star = 0.5/b;
s = 1.2*0.36; 
lref = 1.2*1.43;

dose = 5*10^8;


Results = cell(11,6);
count = 1;
tF = 42;
IC = 'Type2_10mm.mat';
for l = [0.1,0.5,1, 1.1, 1.2, 1.4, 1.6, 1.8, 2,4,8];
params = [DC,a,b,d,j,k,m,q,DT_star,u_star,s,l,dose];
[t, U, V, Tdim, Cdim] = SolvePDE(IC,params,tF);
Results{count,1} = params;
Results{count,2} = t;
Results{count,3} = U;
Results{count,4} = V;
Results{count,5} = Tdim;
Results{count,6} = Cdim;
count = count +1
end
save('lSims.mat', 'Results')
%% 3 -- m
clear all
DC = 0.0138;
b = (4*pi/3)*1e-9;
d = .8*2.25; 
j = .8*0.6; 
k = 1.2*2.019e7;
mref = 1.2*0.293;
q = 1.2*1.53e-9;
a = 0.25;
DT_star = 10^(-5);
u_star = 0.5/b;
s = 1.2*0.36; 
l = 1.2*1.43;

dose = 5*10^8;

Results = cell(8,6);
count = 1;
tF = 42;
IC = 'Type2_10mm.mat';
for m =[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8];
params = [DC,a,b,d,j,k,m,q,DT_star,u_star,s,l,dose];
[t, U, V, Tdim, Cdim] = SolvePDE(IC,params,tF);
Results{count,1} = params;
Results{count,2} = t;
Results{count,3} = U;
Results{count,4} = V;
Results{count,5} = Tdim;
Results{count,6} = Cdim;
count = count +1
end
save('mSims.mat', 'Results')
%% 4 -- j

DC = 0.0138;
b = (4*pi/3)*1e-9;
d = .8*2.25; 
jref = .8*0.6; 
k = 1.2*2.019e7;
m = 1.2*0.293;
q = 1.2*1.53e-9;
a = 0.25;
DT_star = 10^(-5);
u_star = 0.5/b;
s = 1.2*0.36; 
l = 1.2*1.43;

dose = 5*10^8;
Results = cell(9,6);
count = 1;
tF = 42;
IC = 'Type2_10mm.mat';
for j = [5e-2, 1e-1,2e-1,3e-1,4e-1, 5e-1, 6e-1, 7e-1,8e-1, 9e-1,1] ;
params = [DC,a,b,d,j,k,m,q,DT_star,u_star,s,l,dose];
[t, U, V, Tdim, Cdim] = SolvePDE(IC,params,tF);
Results{count,1} = params;
Results{count,2} = t;
Results{count,3} = U;
Results{count,4} = V;
Results{count,5} = Tdim;
Results{count,6} = Cdim;
count = count +1
end
save('jSims.mat', 'Results')

%% 5 -- q
clear all
DC = 0.0138;
b = (4*pi/3)*1e-9;
d = .8*2.25; 
j = .8*0.6; 
k = 1.2*2.019e7;
m = 1.2*0.293;
q = 1.2*1.53e-9;
a = 0.25;
DT_star = 10^(-5);
u_star = 0.5/b;
s = 1.2*0.36; 
l = 1.2*1.43;

dose = 5*10^8;

Results = cell(9,6);
count = 1;
tF = 42;
IC = 'Type2_10mm.mat';
for  q = [1e-11,5e-11, 1e-10, 5e-10, 1e-9, 5e-9,1e-8,5e-8];
params = [DC,a,b,d,j,k,m,q,DT_star,u_star,s,l,dose];
[t, U, V, Tdim, Cdim] = SolvePDE(IC,params,tF);
Results{count,1} = params;
Results{count,2} = t;
Results{count,3} = U;
Results{count,4} = V;
Results{count,5} = Tdim;
Results{count,6} = Cdim;
count = count +1
end
save('qSims.mat', 'Results')

%% 6 -- d
clear all
DC= 0.0138;
b = (4*pi/3)*1e-9;
d = .8*2.25; 
j = .8*0.6; 
k = 1.2*2.019e7;
m = 1.2*0.293;
q = 1.2*1.53e-9;
a = 0.25;
DT_star = 10^(-5);
u_star = 0.5/b;
s = 1.2*0.36; 
l = 1.2*1.43;

dose = 5*10^8;

Results = cell(14,6);
count = 1;
tF = 42;
IC = 'Type2_10mm.mat';
for  d = [.0125, .125,0.25, 0.5, 0.75, 1.25, 2.5, 5, 8, 10,12, 25,50,100];
params = [DC,a,b,d,j,k,m,q,DT_star,u_star,s,l,dose];
[t, U, V, Tdim, Cdim] = SolvePDE(IC,params,tF);
Results{count,1} = params;
Results{count,2} = t;
Results{count,3} = U;
Results{count,4} = V;
Results{count,5} = Tdim;
Results{count,6} = Cdim;
count = count +1
end
save('dSims.mat', 'Results')

%% 7 -- k
clear all
DC = 0.0138;
b = (4*pi/3)*1e-9;
d = .8*2.25; 
j = .8*0.6; 
k = 1.2*2.019e7;
m = 1.2*0.293;
q = 1.2*1.53e-9;
a = 0.25;
DT_star = 10^(-5);
u_star = 0.5/b;
s = 1.2*0.36; 
l = 1.2*1.43;

dose = 5*10^8;

Results = cell(8,6);
count = 1;
tF = 42;
IC = 'Type2_10mm.mat';
for  k = [1e4,1e5,1e6,1e7,1e8, 1e9,1e10,1e11];
params = [DC,a,b,d,j,k,m,q,DT_star,u_star,s,l,dose];
[t, U, V, Tdim, Cdim] = SolvePDE(IC,params,tF);
Results{count,1} = params;
Results{count,2} = t;
Results{count,3} = U;
Results{count,4} = V;
Results{count,5} = Tdim;
Results{count,6} = Cdim;
count = count +1
end
save('kSims.mat', 'Results')
%% 8 -- s
clear all
DC = 0.0138;
b = (4*pi/3)*1e-9;
d = .8*2.25; 
j = .8*0.6; 
k = 1.2*2.019e7;
m = 1.2*0.293;
q = 1.2*1.53e-9;
a = 0.25;
DT_star = 10^(-5);
u_star = 0.5/b;
s = 1.2*0.36; 
l = 1.2*1.43;

dose = 5*10^8;

Results = cell(9,6);
count = 1;
tF = 42;
IC = 'Type2_10mm.mat';
for  s = [.001,.01, 0.1, 0.25, 0.5, 1, 2, 4,10];
params = [DC,a,b,d,j,k,m,q,DT_star,u_star,s,l,dose];
[t, U, V, Tdim, Cdim] = SolvePDE(IC,params,tF);
Results{count,1} = params;
Results{count,2} = t;
Results{count,3} = U;
Results{count,4} = V;
Results{count,5} = Tdim;
Results{count,6} = Cdim;
count = count +1
end
save('sSims.mat', 'Results')

%% 1 -- DC
clear all
DCref = 0.0138;
b = (4*pi/3)*1e-9;
d = .8*2.25; 
j = .8*0.6; 
k = 1.2*2.019e7;
m = 1.2*0.293;
q = 1.2*1.53e-9;
a = 0.25;
DT_star = 10^(-5);
u_star = 0.5/b;
s = 1.2*0.36; 
l = 1.2*1.43;

dose = 5*10^8;
Results = cell(11,6);
count = 1;
tF = 42;
IC = 'Type2_10mm.mat';
for DC = [1e-4 5e-4 1e-3 2.5e-5 5e-3 7.5e-3 1e-2 2.5e-2 5e-2 7.5e-2 1e-1];
params = [DC,a,b,d,j,k,m,q,DT_star,u_star,s,l,dose];
[t, U, V, Tdim, Cdim] = SolvePDE(IC,params,tF);
Results{count,1} = params;
Results{count,2} = t;
Results{count,3} = U;
Results{count,4} = V;
Results{count,5} = Tdim;
Results{count,6} = Cdim;
count = count +1
end
save('DCSims.mat', 'Results')