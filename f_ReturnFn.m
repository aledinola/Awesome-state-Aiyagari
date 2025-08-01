function F = f_ReturnFn(d,aprime,a,e,age,r,alpha,delta,pen,gamma,sigma)
% INPUTS
%   d:       Hours worked
%   aprime:  Next-period's assets
%   a:       Current period assets
%   e:       Labor efficiency shock
%   age:     Age of individual: young or old
% TOOLKIT NOTATION
% (d,aprime,a,z), where z = [e,age]

F = -inf;

w = fun_prices(r,alpha,delta);

income = (w*e*d)*(age==1)+pen*(age==2)+r*a;

c = income+a-aprime; % Budget Constraint

if c>0
    % NOTE: 0<l<1 is already built into the grid
    % WARNING: this will not work if crra=1
    inside = (c^gamma)*((1-d)^(1-gamma));
    F = inside^(1-sigma)/(1-sigma);
end

end