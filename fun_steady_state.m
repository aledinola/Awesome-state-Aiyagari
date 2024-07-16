function [KN_bar,V,Policy,Dist,ED,agg,mom] = fun_steady_state(KN_0,par)

% PURPOSE:
% fun_steady_state computes the steady-state of the model and
% calculates model moments and aggregate quantities.
% INPUTS:
% KN_0: Initial condition for equilibrium price or capital-labor ratio in
%       corporate sector; either a scalar or a 1*2 vector.
%
% OUTPUTS:
% K_L:    capital-labor ratio in equilibrium
% V:      value function, dim: (na,nz)
% Policy: Structure with policy functions
% Dist:   stationary distribution, dim: (na,nz)
%

if numel(KN_0)~=1 && numel(KN_0)~=2
    error('KN_0 must be either a scalar or a two-element vector')
end

switch par.do_GE
    case 0
        % Partial equilibrium
        [ED,V,Policy,Dist,agg] = excess_demand(KN_0,par);
        KN_bar = KN_0;
    case 1
        % Solve general equilibrium conditions:
        % TODO: try fsolve, fminsearch, etc
        options = optimset('Display','iter','TolX',1e-5);
        [KN_bar,ED,flag] = fzero(@(KN) excess_demand(KN,par),KN_0,options);
        if flag<0
            error('FZERO failed, steady-state not found!')
        end
        fprintf(" \n")
        fprintf("Equilibrium found! \n")
        fprintf("Capital-labor ratio = %f \n",KN_bar)
        fprintf("Error (ED)          = %f \n",ED)
        
        % Now that the equilibrium prices (actually, the equilibrium K/N) have been find,
        % compute V and Policy
        % NOTE: Policy is a structure, V and Dist are multidim arrays
        [~,V,Policy,Dist,agg] = excess_demand(KN_bar,par);
        
    otherwise
        error('do_GE is not specified correctly')
end

% Compute some model moments
mom = compute_targets(KN_bar,Policy,Dist,par);

end %end function "compute_steady_state"

