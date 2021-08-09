%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation Project
% Day1_Setup.m
% Yongseok Kim - Indiana University
% 2021 Summer Summer School on Structural Estimation in Corporate Finance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%% Set up grids and discretization and payoffs')
disp(' ')
addpath('build/function')
tic

%%%%%%%%%% set up grid for endogenous capitals

kgrid = linspace(log(kmin),log(kmax),knum);
kgrid = exp(kgrid)';

%%%%%%%%%% set up loglinear grid & discretization for exogenous profitability process

[zgrid, pr_mat_z] = tauchen(0, sig, rho, znstdev, znum);
zgrid = exp(zgrid)'; %convert to levels, mnum x 1
pr_mat_z(pr_mat_z<0) = 0; %fix rounding issue leading to negatives in tauchen() function

%%%%%%%%%% set up index of the state space
grid_val = zeros(statenum,2); %mnum*pnum x 2, with (i,1) = value of k, (i,2) = value of z
grid_ind = zeros(statenum,2); %mnum*pnum x 2, with (i,1) = index of k, (i,2) = index of z

%insert values
grid_val(:,1) = kron(ones(znum,1),kgrid);
grid_val(:,2) = kron(zgrid,ones(knum,1));

%insert indexes
grid_ind(:,1) = kron(ones(znum,1),(1:knum)');
grid_ind(:,2) = kron((1:znum)',ones(knum,1));
%%%%%%%%%% set up payoff matrices

%Emat is statenum x knum, with (i,j) static payoff to state i, policy j
Estarmat = repmat(grid_val(:,2),1,knum) .* (repmat(grid_val(:,1),1,knum).^alpha) ...
    - repmat(kgrid',statenum,1) ...
    + (1-delta) * repmat(grid_val(:,1),1,knum);

Emat = Estarmat;
Emat(Estarmat<0) = Estarmat(Estarmat<0) * (1 + lambda);

toc
disp(' ')