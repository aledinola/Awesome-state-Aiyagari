%% Replication of Appendix E.2 of GKKOC "Use it or Lose it", QJE 2023
clear;clc;close all
format long g
addpath(genpath('C:\Users\aledi\Documents\GitHub\VFIToolkit-matlab'))

%% Set flags
do_GE = 1;

%% Set grid sizes
n_a = 601; % assets
n_d = 51;  % labor supply

%% Set parameters to calibrate
Params.r   = 0.0306151766296052; % Initial guess for interest rate
Params.pen = 0.517906812023212;  % Pension benefits
Params.e_awesome = 265;          % Awesome productivity level

% Calibration targets
%   pen calibrated to match P_to_Y
%   e_awesome calibrated to match wealth_top1
Params.P_to_Y_data      = 0.049; % Pensions to GDP equal to 4.9%
Params.wealth_top1_data = 0.30; % Wealth share top 1% should be 30%

%% Set toolkit options

% --- Value functions options
vfoptions=struct(); 
vfoptions.lowmemory     = 0;
vfoptions.verbose       = 0;
vfoptions.tolerance     = 1e-7;
vfoptions.maxiter       = 500;
vfoptions.howards       = 80; 
vfoptions.maxhowards    = 500;
vfoptions.howardsgreedy = 0;
vfoptions.gridinterplayer = 0;
vfoptions.ngridinterp     = 15;
%vfoptions.divideandconquer = 0;

% Distribution options
simoptions=struct(); % Use default options for solving for stationary distribution
simoptions.verbose = 0;
simoptions.tolerance = 1e-8;
simoptions.gridinterplayer = vfoptions.gridinterplayer;
simoptions.ngridinterp     = vfoptions.ngridinterp;

% Heteroagentoptions
heteroagentoptions = struct();
heteroagentoptions.verbose=1; % verbose means that you want it to give you feedback on what is going on
heteroagentoptions.toleranceGEprices=1e-6; % default is 1e-4
heteroagentoptions.toleranceGEcondns=1e-6; % default is 1e-4
heteroagentoptions.fminalgo = 1;  % 0=fzero, 1=fminsearch, 8=lsqnonlin 
heteroagentoptions.maxiter = 2;

%% Set economic parameters

% Source: GKKOC (2023), Table II.
Params.sigma = 4.0;  % CRRA utility function
Params.gamma = 0.445; % Consumption share in utility
Params.alpha = 0.4;   % Capital share of output
Params.delta = 0.05;  % Depreciation rate

% Source: Castañeda, Díaz-Giménez, and Ríos-Rull, 2003, Table 3.
Params.beta       = 0.924; % Discount factor
Params.prob_ret   = 0.022; % prob. worker --> retired
Params.prob_death = 0.066; % prob. retired --> worker
Params.phi_2      = 0.525; % Shifting probability

%% Set grids and shocks

% - Grid for labor efficiency shock e
% GKKOC (2023) state that the value for e1,e2,e3 are taken from Table 5 of
% Castaneda et al. (2003)

e_grid = [1.00, 3.15, 9.78, Params.e_awesome]';
n_e    = length(e_grid);

% - Transition matrix PI(e,e'), for labor efficiency shock e. GKKOC (2023)
% state that PI(e,e') is taken from Table 4 Castaneda et al. (2003)
pi_e(1,:) = [96.24 1.14 .39 .006]/100;
pi_e(2,:) = [3.07 94.33 .37 .000]/100;
pi_e(3,:) = [1.50 .43 95.82 .020]/100;
pi_e(4,:) = [10.66 .49 6.11 80.51]/100;
pi_e = pi_e./sum(pi_e,2);

% - Invariant or stationary distribution of PI(e,e'), called gamma_e
aux = pi_e^1000;
gamma_e = aux(1,:)';
gamma_e = gamma_e/sum(gamma_e);

% - Initial distribution G_e, from which newborn draw e
G_e = funCDR(gamma_e,Params.phi_2);
%G_e = gamma_e;

% - Set grid for age
n_age = 2;
age_grid = [1,2]'; % 1= YOUNG, 2=OLD
pi_age = [1-Params.prob_ret, Params.prob_ret
          Params.prob_death, 1-Params.prob_death];
auzz = pi_age^1000;
prob_age = auzz(1,:)'; % marginal distribution of age
mass_retired = prob_age(2); % mass of agents with age=2 (i.e. retirees)

% Grid for assets
a_curve = 3.0;
a_min   = 0.0;
a_max   = 2000;

a_grid  = a_min+(a_max-a_min)*(linspace(0,1,n_a).^a_curve)';

%grid for labor
d_grid = linspace(0,0.999,n_d)';

%% Recast model in toolkit notation
% The state variables are (a,e,age), age=Young,Retired, and if age=R, e 
% does not matter. So the total number of exogenous states is n_e+1: All 
% the values of e when age=Y and then an arbitrary value of e when age=R.
% Set n_z = [n_e+1,1], size(z_grid)=[n_e+1,2], size(pi_z)=[n_e+1,n_e+1]
% z_grid = [e_1,Young
%           e_2,Young
%           ...
%           e_ne,Young
%           NaN,Retired]

n_z    = [n_e+1,1];
z_grid = [e_grid,age_grid(1)*ones(n_e,1)
          0,     age_grid(2)];

pi_z = [(1-Params.prob_ret)*pi_e, Params.prob_ret*ones(n_e,1);
        Params.prob_death*G_e',   (1-Params.prob_death)];

% Check size of pi_z
if ~isequal(size(pi_z),[n_e+1,n_e+1])
    error('pi_z has wrong size!')
end
if n_z(1)~=(n_e+1)
    error('Num. of grid points for exog shock must be n_e+1')
end

%% Set return function
DiscountFactorParamNames={'beta'};

ReturnFn = @(d,aprime,a,e,age,r,alpha,delta,pen,gamma,sigma) ...
    f_ReturnFn(d,aprime,a,e,age,r,alpha,delta,pen,gamma,sigma);

%% Functions to be evaluated
FnsToEvaluate.K = @(d,aprime,a,e,age) a; 
FnsToEvaluate.L = @(d,aprime,a,e,age) e*d*(age==1); %Only young supply labor
FnsToEvaluate.Pensions = @(d,aprime,a,e,age,pen) pen*(age==2); 
% If you are retired (age=2) you earn pension pen (otherwise it is zero).

%% GE conditions
% Parameters calibrated internally to match targets:
% r --> capital market residual equal to 0 (real general equilibrium)
% pen --> ratio of transfers to Y equal to 0.049
% e4 --> share of wealth held by top 1% equal to 0.30

% Declare GE conditions (or calibration targets)
GeneralEqmEqns.CapitalMarket = @(K,L,r,alpha,delta) r-alpha*(K/L)^(alpha-1)+delta;
GeneralEqmEqns.target_pen = @(K,L,Pensions,alpha,P_to_Y_data) Pensions-P_to_Y_data*(K^alpha*L^(1-alpha));
GeneralEqmEqns.target_wealth_top1 = @(wealth_top1,wealth_top1_data) wealth_top1-wealth_top1_data;

% Declare names of GE parameters
GEPriceParamNames={'r','pen','e_awesome'};
heteroagentoptions.constrainpositive = {'pen'};
heteroagentoptions.multiGEweights = [1,1,1];
heteroagentoptions.CustomModelStats = @(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_z,heteroagentoptions,vfoptions,simoptions) ...
    fun_custom_stats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_z,heteroagentoptions,vfoptions,simoptions);

%% Solve model
disp('===================================================================')
disp('TOOLKIT')

if do_GE==1
    disp('SOLVING FOR GENERAL EQUILIBRIUM')
    [p_eqm,~,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
    disp('GeneralEqmCondn:')
    disp(GeneralEqmCondn)
    if heteroagentoptions.maxiter>0
    disp('Updating the GE parameters with values found')
    for ii=1:numel(GEPriceParamNames)
        name_ii = GEPriceParamNames{ii};
        Params.(name_ii) = p_eqm.(name_ii); 
    end
    end
else
    disp('PARTIAL EQUILIBRIUM')
end

%% Value function iteration

tic
[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
time_vfi = toc;
fprintf('Time to compute VFI:  %f \n',time_vfi)

%% Distribution

tic
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);
time_dist = toc;
fprintf('Time to compute Dist:  %f \n',time_dist)

%% Model moments

% Additional functions to be evaluated
FnsToEvaluate.Income = @(d,aprime,a,e,age,r,alpha,delta,pen) f_Income(d,aprime,a,e,age,r,alpha,delta,pen);
FnsToEvaluate.Consumption = @(d,aprime,a,e,age,r,alpha,delta,pen) f_Consumption(d,aprime,a,e,age,r,alpha,delta,pen);

tic
AllStats=EvalFnOnAgentDist_AllStats_Case1(StationaryDist,Policy,FnsToEvaluate,...
    Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
time_allstats = toc;
fprintf('Time to compute AllStats:  %f \n',time_allstats)

TopWealthShares = 1-AllStats.K.LorenzCurve([80,95,99]); % Need the 20,5,and 1 top shares

% --- Obtain policy functions in values, from indexes
PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_z,d_grid,a_grid,vfoptions);

% Wage consistent with equilibrium r
Params.w = fun_prices(Params.r,Params.alpha,Params.delta);
KK = AllStats.K.Mean;
LL = AllStats.L.Mean;
CC = AllStats.Consumption.Mean;
YY = KK^Params.alpha*LL^(1-Params.alpha);
PP = Params.pen*mass_retired;
walras = abs(CC+Params.delta*KK-YY-PP);

% Recompute GE conditions
GE_CapitalMarket = GeneralEqmEqns.CapitalMarket(KK,LL,Params.r,Params.alpha,Params.delta);
GE_wealth_top1 = GeneralEqmEqns.target_wealth_top1(TopWealthShares(3),Params.wealth_top1_data);

%% Look at results
disp('===================================================================')

disp('==================================================================')
disp('PRICES')
fprintf('r   : %f \n',Params.r)
fprintf('w   : %f \n',Params.w)
disp('GE CONDITIONS')
fprintf('Capital Market : %f \n',GE_CapitalMarket)
fprintf('Walras resid   : %f \n',walras)
disp('CALIBRATION TARGETS')
fprintf('e_awesome : %f \n',abs(TopWealthShares(3)-Params.wealth_top1_data))
fprintf('pension   : %f \n',abs(PP/YY-Params.P_to_Y_data))
disp('QUANTITIES')
fprintf('K/Y        : %f \n',AllStats.K.Mean/YY)
fprintf('K/L        : %f \n',AllStats.K.Mean/AllStats.L.Mean)
fprintf('K          : %f \n',AllStats.K.Mean)
fprintf('L          : %f \n',AllStats.L.Mean)
fprintf('Y          : %f \n',YY)
fprintf('Pensions/Y : %f \n',PP/YY)
%fprintf('Pensions/Y : %f \n',pen_to_Y)
disp('INEQUALITY')
fprintf('Gini consumption    : %f \n',AllStats.Consumption.Gini)
fprintf('Gini income         : %f \n',AllStats.Income.Gini)
fprintf('Gini wealth         : %f \n',AllStats.K.Gini)
fprintf('Share wealth top 20 : %f \n',TopWealthShares(1))
fprintf('Share wealth top 5  : %f \n',TopWealthShares(2))
fprintf('Share wealth top 1  : %f \n',TopWealthShares(3))

% Plot policy for assets
pol_d    = squeeze(PolicyValues(1,:,:));
pol_aprime = squeeze(PolicyValues(2,:,:));
StationaryDist_a = sum(StationaryDist,2);

% Policy for consumption
pol_c = zeros(n_a,n_z(1));
for z_c=1:n_z(1)
    e_val = z_grid(z_c,1);
    age_val = z_grid(z_c,2);
    pol_c(:,z_c) = f_Consumption(pol_d(:,z_c),pol_aprime(:,z_c),a_grid,e_val,age_val,...
        Params.r,Params.alpha,Params.delta,Params.pen);
end

figure
plot(a_grid,a_grid,'--')
hold on
plot(a_grid,pol_aprime(:,1))
hold on
plot(a_grid,pol_aprime(:,n_e))
legend('e_1','e_{ne}','Location','best')
title('Policy function for assets')

figure
plot(a_grid,pol_d(:,1))
hold on
plot(a_grid,pol_d(:,n_e))
legend('e_1','e_{ne}','Location','best')
title('Policy function for labor supply')

figure
plot(a_grid,pol_c(:,1))
hold on
plot(a_grid,pol_c(:,n_e))
legend('e_1','e_{ne}','Location','best')
title('Policy function for consumption')

figure
plot(a_grid,StationaryDist_a,'LineWidth',2)
title('Marginal distribution of assets')