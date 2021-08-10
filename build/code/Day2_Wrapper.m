%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation Project
% Day2_Wrapper.m
% Yongseok Kim - Indiana University
% 2021 Summer Summer School on Structural Estimation in Corporate Finance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('~/Courseworks/KelleyPhd/MitsuiCenter/Project/')
addpath('build/code/')
addpath('build/function/')
clear all; close all; clc;

disp('%%%%%%%%%%%%% Solving the investment problem')
disp(' ')

global beta rho sig lambda knum znum znstdev statenum ...
    soltol maxit ...
    firmnum Terg Tsim Ttot zinit kinit

%%%%%%%%%% set fixed model parameters
beta = 0.96;
rho = 0.75;
sig = 0.30;
lambda = 0.05;

%%%%%%%%%% set solution parameters
%grid parameters
knum = 1000; %number of capital grid points
znum = 35; %number of profitability shock grid points
znstdev = 3; %number of standard deviations to cover around steady state profitability
statenum = znum*knum; %total number of state grid points

%solution parameters
soltol = 1e-7; %tolerance on model solution
maxit = 1000; %max iterations on model solution

%simulation parameters
firmnum = 5000; %number of firms in simulation
Terg = 50; %number of burn-in periods per firm
Tsim = 15; %number of periods for useful simulation
Ttot = Terg + Tsim+1; %total number of periods in simulation
rng(345891); %set random seed for simulation draws
zinit = floor(znum/2); %initial point for m in simulation
kinit = floor(knum/2); %inital point for p_{-1} in simulation

%%%%%%%%% compute data moments
load('build/input/RealData.mat')
prof = RealData(:,3);
inv = RealData(:,4);

global datamom nummom
datamom = Day2_compute_moments(prof, inv);
nummom = length(datamom);

%%%%%%%%%% minimize objective function value
InitTemp=1000;
jumpsizes = [0.01, 0.01];
Nparams = 2;

anneal_opts=struct(...
    'CoolSched',@(T) (0.85*T), ...
    'Generator',@(x) (x+jumpsizes.*randn(Nparams,1)), ...
    'InitTemp', InitTemp, ...
    'MaxConsRej', 100, ...
    'MaxSuccess', 30, ...
    'MaxTries', 20, ...
    'StopTemp',0.00001 , ...
    'StopVal',-Inf, ...
    'Verbosity',2);

%  Search for parameter values 

disp('Finding SMM estimates')
disp('     Simulated annealing:')

tic

% First use simulated annealing to get close to the global minimum
[params, loss] = ...
    anneal_1(@(params)Day2_SMM(params), [0.5, 0.05], jumpsizes, anneal_opts);
