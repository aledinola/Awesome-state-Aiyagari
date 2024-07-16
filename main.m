%% Replication of Appendix E.2 of GKKOC "Use it or Lose it", QJE 2023
clear
clc
close all
format long g
% This is the folder where the VFI toolkit files are saved
folder1 = 'C:\Users\aledi\Documents\GitHub\VFIToolkit-matlab';
%folder2 = fullfile('..','tools');
addpath(genpath(folder1))

%% Set flags
do_GE = 0;     % 0=No general equilibrium,1=general eqm
m_toolkit = 2; % Toolkit method: 1=default, 2=correlated shocks

%% Set economic parameters
% Mostly taken from Table II of GKKOC (2023)
Params.alpha = 0.4;   % Capital share of output
Params.delta = 0.05;  % Depreciation rate
Params.crra  = 4.0;   % CRRA utility function
Params.gamma = 0.445; % Consumption share in utility
Params.beta  = 0.924; % Discount factor
Params.pen   = 1;     % Pension benefits

% Demographics
Params.prob_ret   = 0.022;
Params.prob_death = 0.066;

% Initial guess for K/L ratio
Params.K_to_L = 14.04; % 12 gives err>0
r = Params.alpha*Params.K_to_L^(Params.alpha-1)-Params.delta;
w = (1-Params.alpha)*Params.K_to_L^Params.alpha;

%% Set grids and shocks

% - Grid for labor efficiency shock e
% GKKOC (2023) state that the value for e1,e2,e3 are taken from Table 5 of
% Castaneda et al. (2003)
e_awesome = 15;
e_grid = [1.00, 3.15, 9.78, e_awesome]';
n_e    = length(e_grid);

% - Transition matrix PI(e,e'), for labor efficiency shock e. GKKOC (2023)
% state that PI(e,e') is taken from Table 4 Castaneda et al. (2003)
pi_e(1,:) = [96.24 1.14 .39 .006]/100;
pi_e(2,:) = [3.07 94.33 .37 .000]/100;
pi_e(3,:) = [1.50 .43 95.82 .020]/100;
pi_e(4,:) = [10.66 .49 6.11 80.51]/100;
pi_e = pi_e./sum(pi_e,2);

% - Initial distribution G_e, from which newborn draw e
aux = pi_e^1000;
G_e = aux(1,:)';
G_e = G_e/sum(G_e);

% - Set grid for age
n_age = 2;
age_grid = [1,2]'; % 1= YOUNG, 2=OLD

% Grid for assets
a_curve = 3.0;
a_min   = 0.0;
a_max   = 250;
n_a     = 601;
a_grid  = a_min+(a_max-a_min)*(linspace(0,1,n_a).^a_curve)';

%grid for labor
n_d    = 51;
d_grid = linspace(0.001,0.999,n_d)';

if m_toolkit==1
    tic
    [p_eqm,V,Policy,StationaryDist,AggVars,AllStats,TopWealthShares] = ...
        solve_toolkit_default(Params,e_grid,age_grid,G_e,a_grid,d_grid,n_e,...
        n_age,n_a,n_d,pi_e,do_GE);
    toc
elseif m_toolkit==2
    tic
    [p_eqm,V,Policy,StationaryDist,AggVars,AllStats,TopWealthShares] = ...
        solve_toolkit_CORR_SHOCKS(Params,e_grid,age_grid,G_e,a_grid,d_grid,n_e,...
        n_age,n_a,n_d,pi_e,do_GE);
    toc
else
    error('m_toolkit out of bounds')
end

disp('p_eqm')
disp(p_eqm)
