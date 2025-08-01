function c = f_Consumption(d,aprime,a,e,age,r,alpha,delta,pen)
% f_Consumption computes consumption from the budget constraint
% INPUTS
%   d:       Hours worked
%   aprime:  Next-period's assets
%   a:       Current period assets
%   e:       Labor efficiency shock
%   age:     Age of individual: young or old
% TOOLKIT NOTATION
% (d,aprime,a,z), where z = [e,age]


w = fun_prices(r,alpha,delta);

income = (w*e*d)*(age==1)+pen*(age==2)+r*a;

c = income+a-aprime; % Budget Constraint

end %end function