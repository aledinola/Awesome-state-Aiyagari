function c = ConsumptionFn(d,aprime,a,e,age,K_to_L,alpha,delta,pen)
% INPUTS
%   d:       Hours worked
%   aprime:  Next-period's assets
%   a:       Current period assets
%   e:       Labor efficiency shock
%   age:     Age of individual: young or old
% TOOLKIT NOTATION
% (d,aprime,a,z), where z = [e;age]

r = alpha*K_to_L^(alpha-1)-delta;
w = (1-alpha)*K_to_L^alpha;

income = (w*e*d)*(age==1)+pen*(age==2)+r*a;

c = income+a-aprime; % Budget Constraint

end