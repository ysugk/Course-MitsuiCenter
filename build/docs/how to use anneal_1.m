% This snippet of code explains how to implement anneal_1.m


% To make simulated annealing search more thoroughly (and slowly), try
% increasing the values of MaxConsRej, MaxSuccess, MaxTries, and change
% CoolSched factor from 0.8 to some higher number that's still < 1.
% Setting the values below is a bit of an art, and it depends on how many
% parameters you're estimating.
InitTemp=1000;

anneal_opts=struct(...
    'CoolSched',@(T) (0.85*T), ...
    'Generator',@(x) (x+jumpsizes.*randn(Nparams,1)), ...
    'InitTemp', InitTemp, ...
    'MaxConsRej', 15, ...
    'MaxSuccess', 15, ...
    'MaxTries', 25, ...
    'StopTemp',0.00001 , ...
    'StopVal',-Inf, ...
    'Verbosity',2);

%  Search for parameter values 

disp('Finding SMM estimates')
disp('     Simulated annealing:')

tic

% First use simulated annealing to get close to the global minimum

[ParamsEst_temp,FVAL]=...
    anneal_1(@(Params)Score(Params,M_emp,WeightMat,ParamsOther,NFirms,NYearsPerFirm,NSim),...
    ParamsEst_guess,jumpsizes,anneal_opts);










