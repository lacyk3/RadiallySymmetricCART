%% ---Simulations for Tumor growth maps
clear all; close all; clc
%% Set up 
% -- Build grid of parameter values
u_stars = [.01, .05, .1, .2, .4, .6, .8, .9,.005];
DT_stars = [1e-5, 2e-5, 4e-5, 6e-5, 8e-5, ...
            1e-4, 1.2e-4, 1.4e-4,5e-6];

[U_star_grid, DT_star_grid] = meshgrid(u_stars, DT_stars); 

U_star_line = reshape(U_star_grid, 1, numel(U_star_grid)); 
DT_star_line = reshape(DT_star_grid, 1, numel(DT_star_grid)); 
% -- Vector to save values
TVDTs = zeros(1, numel(U_star_line)); % tumor volume doubling time
TTDs = zeros(1, numel(U_star_line)); % time to detection
TBsatDetection = zeros(1, numel(U_star_line)); % tumor burden at detection
indices = string(1:numel(U_star_line));


parfor i= 1:numel(U_star_line)
    b = (4*pi/3)*1e-9; % cm^3/cell
    u_detection = 10^8*b; 
    beta_T = 1;
    
    dr = 0.015;
    r = [0:dr:1];
    r1 = round(0.5*length(r));
    Uinit =  exp(1 - 1./(1 - (r(1:r1 - 1)/r(r1)).^4))/beta_T;

    %% --- a = 0.25
    % set up domain
    a = 0.25;
    dt = a/10;
    T = 1000*a;
    N = round(T/dt + 1);

    u_star = U_star_line(i); DT_starnd = DT_star_line(i)/a;
    U = zeros(length(r), N);
    U(1:r1 - 1,1) = Uinit;
    m=matfile(sprintf('output%d.mat', i),'writable',true)
    m.a = a;
    m.u_star = u_star;
    m.DT_starnd = DT_starnd;
    options = optimset('Display','off', 'FinDiffType', 'central', 'UseParallel', false);
flag = 0;
tic
    for n = 2:N
        detectable = U(:,n-1)>u_detection;
        tumorsize = r(detectable);
        
        if length(tumorsize)>1
            % when detectable tumor reaches 1 cm, record time and tumor
            % burden
            if flag == 0
                if tumorsize(end) >= 1
                    TTD = dt*n/a
                    TTDs(i) = TTD;
                    m.TTD = dt*n/a;
                    TBatDetection = sum(4*pi*U(:,n-1).*(r.').^2)*dr  
                    TBsatDetection(i) = TBatDetection;
                    m.TBatDetection = sum(4*pi*U(:,n-1).*(r.').^2)*dr;
                    flag = 1
                end
            end
        % After that, when detectable tumor doubles in volume,
        % record time and end
            if tumorsize(end)>= 2^(1/3)
                TVDT = dt*n/a - TTD
                TVDTs(i) = TVDT;
                m.TVDT = TVDT;
                toc
                break
            end
        end
        % If tumor density has negative values, end the sim
        if min(real(U(:,n-1)))<0
            break 
        end
        % If tumor has reached edge domain, extend the domain 
        if  U(end,n-1) > 0
            r = [0:dr:r(end)*1.2];
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
    m.U = U;
    m.r = r;
end

TVDTs_grid = reshape(TVDTs, size(U_star_grid));
TTDs_grid = reshape(TTDs, size(U_star_grid));
TBsatDetection_grid = reshape(TBSatDetection, size(U_star_grid));

save('a25.mat','TVDTs_grid','TTDs_grid','TBsatDetection_grid','U_star_grid','DT_star_grid')

%% --- a = 0.125
%a = 0.125;
%save('a125.mat','TVDTs_grid','TTDs_grid','TBsatDetection_grid','U_star_grid','DT_star_grid')
%% --- a = 0.025
%a = 0.025;
%save('a025.mat','TVDTs_grid','TTDs_grid','TBsatDetection_grid','U_star_grid','DT_star_grid')