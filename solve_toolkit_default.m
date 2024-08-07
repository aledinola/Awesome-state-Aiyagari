function [V,Policy,StationaryDist,AggVars,AllStats,TopWealthShares] = ...
    solve_toolkit_default(Params,e_grid,age_grid,G_e,a_grid,d_grid,n_e,...
    n_age,n_a,n_d,pi_e)

%% Solve model using toolkit with inefficient method
% z_grid = [e_grid;age_grid]
% If e_grid has J points, z_grid will have 2*J points (all combos)
% This is a waste since once retired, e is not a state anymore

verbose = 1;

% - Set grid and transition prob for exogenous state z
n_z    = [n_e,n_age];
z_grid = [e_grid;age_grid];

pi_z = [(1-Params.prob_ret)*pi_e,          Params.prob_ret*eye(n_e);
        Params.prob_death*ones(n_e,1)*G_e', (1-Params.prob_death)*eye(n_e)];


%% Setup Return function
DiscountFactorParamNames={'beta'};

ReturnFn = @(d,aprime,a,e,age,K_to_L,alpha,delta,pen,gamma,crra) ...
    f_ReturnFn(d,aprime,a,e,age,K_to_L,alpha,delta,pen,gamma,crra);

%% Set Functions to evaluate 
% Functions to be evaluated
FnsToEvaluate.K = @(d,aprime,a,e,age) a; 
FnsToEvaluate.L = @(d,aprime,a,e,age) e*d*(age==1); %Only young supply labor
FnsToEvaluate.Income = @(d,aprime,a,e,age,K_to_L,alpha,delta,pen) IncomeFn(d,aprime,a,e,age,K_to_L,alpha,delta,pen);
FnsToEvaluate.Wealth = @(d,aprime,a,e,age) a;
FnsToEvaluate.Consumption = @(d,aprime,a,e,age,K_to_L,alpha,delta,pen) ConsumptionFn(d,aprime,a,e,age,K_to_L,alpha,delta,pen);
FnsToEvaluate.Pensions = @(d,aprime,a,e,age,pen) pen*(age==2); 
% If you are retired (age=2) you earn pension pen (otherwise it is zero).

%%% GE conditions
% General equilibrium condition(s)
% Note: Inputs can be any parameter, price, or aggregate of the FnsToEvaluate
%GeneralEqmEqns.CapitalMarket = @(K_to_L,K,L) K_to_L-K/L;

% Declare names of GE parameters
%GEPriceParamNames={'K_to_L'};

%% VFI

vfoptions=struct(); % Use default options for solving the value function (and policy fn)
%vfoptions.lowmemory = 1;
vfoptions.verbose = verbose;

[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

%% Plot value and policy functions
disp('Size of V:')
disp([n_a,n_e,n_age])
disp(size(V))

disp('Size of Policy:')
disp([2,n_a,n_e,n_age])
disp(size(Policy))

Policy_aprime = squeeze(gather(Policy(2,:,:,:))); % (a,e,age)
Policy_aprime_val = a_grid(Policy_aprime);
Policy_d = squeeze(gather(Policy(1,:,:,:))); % (a,e,age)
Policy_d_val = d_grid(Policy_d);

figure
plot(a_grid,a_grid,'--','LineWidth',2)
hold on
plot(a_grid,Policy_aprime_val(:,1,1),'LineWidth',2)
hold on
plot(a_grid,Policy_aprime_val(:,n_e,1),'LineWidth',2)
legend('45 line','e_1','e_4')
title('Policy a'' for YOUNG')
xlabel('Current-period assets, a')
ylabel('Next-period assets, a'' ')

figure
plot(a_grid,Policy_d_val(:,1,1),'LineWidth',2)
hold on
plot(a_grid,Policy_d_val(:,2,1),'LineWidth',2)
hold on
plot(a_grid,Policy_d_val(:,3,1),'LineWidth',2)
hold on
plot(a_grid,Policy_d_val(:,n_e,1),'LineWidth',2)
legend('e_1','e_2','e_3','e_4')
title('Policy labor, l, for YOUNG')
xlabel('Current-period assets, a')
ylabel('Labor supply, l')

figure
plot(a_grid,a_grid,'--','LineWidth',2)
hold on
plot(a_grid,Policy_aprime_val(:,1,2),'LineWidth',2)
hold on
plot(a_grid,Policy_aprime_val(:,n_e,2),'LineWidth',2)
legend('45 line','e_1','e_4')
title('Policy a'' for RETIRED')
xlabel('Current-period assets, a')
ylabel('Next-period assets, a'' ')

figure
plot(a_grid,Policy_d_val(:,1,2),'LineWidth',2)
hold on
plot(a_grid,Policy_d_val(:,2,2),'LineWidth',2)
hold on
plot(a_grid,Policy_d_val(:,3,2),'LineWidth',2)
hold on
plot(a_grid,Policy_d_val(:,n_e,2),'LineWidth',2)
legend('e_1','e_2','e_3','e_4')
title('Policy labor, l, for RETIRED')
xlabel('Current-period assets, a')
ylabel('Labor supply, l')

%% Stationary Distribution
simoptions = struct(); 
simoptions.verbose = verbose;
simoptions.tolerance = 1e-9;
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);

disp('Size of StationaryDist:')
disp([n_a,n_e,n_age])
disp(size(StationaryDist))

% Marginal PDF over assets
mu_a = gather(sum(StationaryDist,[2,3]));

figure
plot(a_grid,mu_a,'LineWidth',2)
title('Stationary distribution of assets, PDF')
xlabel('Current-period assets, a')

figure
plot(a_grid,cumsum(mu_a),'LineWidth',2)
title('Stationary distribution of assets, CDF')
xlabel('Current-period assets, a')

%% Aggregate variables

AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist,Policy,FnsToEvaluate,...
    Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid);

AllStats=EvalFnOnAgentDist_AllStats_Case1(StationaryDist,Policy,FnsToEvaluate,...
    Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

TopWealthShares=100*(1-AllStats.Wealth.LorenzCurve([80,95,99])); % Need the 20,5,and 1 top shares

K_to_L = AggVars.K.Mean/AggVars.L.Mean;
r = Params.alpha*K_to_L^(Params.alpha-1)-Params.delta;
w = (1-Params.alpha)*K_to_L^Params.alpha;
YY = AggVars.K.Mean^Params.alpha*AggVars.L.Mean^(1-Params.alpha);

% Ratio of pension benefits to GDP
pen_to_Y = AggVars.Pensions.Mean/YY;

%% Display results

disp('==================================================================')
disp('PRICES')
fprintf('r   : %f \n',r)
fprintf('w   : %f \n',w)
disp('QUANTITIES')
fprintf('K          : %f \n',AggVars.K.Mean)
fprintf('L          : %f \n',AggVars.L.Mean)
fprintf('K/L        : %f \n',AggVars.K.Mean/AggVars.L.Mean)
fprintf('K/Y        : %f \n',AggVars.K.Mean/YY)
fprintf('Pensions/Y : %f \n',pen_to_Y)
disp('INEQUALITY')
fprintf('Gini consumption    : %f \n',AllStats.Consumption.Gini)
fprintf('Gini income         : %f \n',AllStats.Income.Gini)
fprintf('Gini wealth         : %f \n',AllStats.Wealth.Gini)
fprintf('Share wealth top 20 : %f \n',TopWealthShares(1))
fprintf('Share wealth top 5  : %f \n',TopWealthShares(2))
fprintf('Share wealth top 1  : %f \n',TopWealthShares(3))

end %end function