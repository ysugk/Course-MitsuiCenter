%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation Project
% Day1_Wrapper.m
% Yongseok Kim - Indiana University
% 2021 Summer Summer School on Structural Estimation in Corporate Finance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global parameters:
% model: beta, rho, sig, lambda
% solution: knum, kmin, kmax, znum, znstdev, statenum, soltol, maxit

%%%%%%%%%% set up the grids, discretization of marg cost, payoff matrices
run('Day2_Setup.m');
%%%%%%%%%% do the VFI loop
run('Day2_VFI_MPQ.m')
