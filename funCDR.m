function [G_e] = funCDR(gamma_e,phi_2)
% This function reqeights the probability of income shocks
% Given invariant distribution of income shocks gamma_e, we shift
% probability mass towards the first column, obtaining a new distribution
% called G_e. Newborn agents draw e from G_e.

n_e = length(gamma_e);
G_e = zeros(n_e,1);

G_e(1) = gamma_e(1)+phi_2*gamma_e(2)+phi_2^2*gamma_e(3)+phi_2^3*gamma_e(4);
G_e(2) = (1-phi_2)*(gamma_e(2)+phi_2*gamma_e(3)+phi_2^2*gamma_e(4));
G_e(3) = (1-phi_2)*(gamma_e(3)+phi_2*gamma_e(4));
G_e(4) = (1-phi_2)*gamma_e(4);

end