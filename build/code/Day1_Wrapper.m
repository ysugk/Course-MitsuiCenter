%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation Project
% Day1_Wrapper.m
% Yongseok Kim - Indiana University
% 2021 Summer Summer School on Structural Estimation in Corporate Finance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('~/Courseworks/KelleyPhd/MitsuiCenter/Project/')
addpath('build/code/')
clear all; close all; clc;

disp('%%%%%%%%%%%%% Solving the investment problem')
disp(' ')

%%%%%%%%%% set model parameters

alpha = 0.5;
delta = 0.05;
beta = 0.96;
rho = 0.75;
sig = 0.30;
lambda = 0.05;

%%%%%%%%%% set solution parameters

%grid parameters
knum = 1000; %number of capital grid points
kmin = 1; %lowest grid point for capital
kmax = 250; %highest grid point for capital
znum = 35; %number of profitability shock grid points
znstdev = 3; %number of standard deviations to cover around steady state profitability
statenum = znum*knum; %total number of state grid points

%solution parameters
soltol = 1e-7; %tolerance on model solution
maxit = 1000; %max iterations on model solution

%simulation parameters
%simulation parameters
firmnum = 5000; %number of firms in simulation
Terg = 50; %number of burn-in periods per firm
Tsim = 15; %number of periods for useful simulation
Ttot = Terg + Tsim+1; %total number of periods in simulation
rng(345891); %set random seed for simulation draws
zinit = floor(znum/2); %initial point for m in simulation
kinit = floor(knum/2); %inital point for p_{-1} in simulation

%some formatting stuff
lwidnum = 2;
fsizenum = 12;

%%%%%%%%%% set up the grids, discretization of marg cost, payoff matrices
run('Day1_Setup.m');
%%%%%%%%%% implement simple VFI
run('Day1_VFI.m');
%%%%%%%%%% simulate the model
run('Day1_Simulate.m');