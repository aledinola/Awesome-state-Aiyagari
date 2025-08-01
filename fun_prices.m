function [w,K_to_L] = fun_prices(r,alpha,delta)
% fun_prices computes w and K/L given interest rate

K_to_L = (alpha/(r+delta))^(1/(1-alpha));
w = (1-alpha)*K_to_L^alpha;

end