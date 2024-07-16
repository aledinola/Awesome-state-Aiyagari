function F = f_ReturnFn(d,aprime,a,e,age,K_to_L,alpha,delta,pen,gamma,crra)
% INPUTS
%   d:       Hours worked
%   aprime:  Next-period's assets
%   a:       Current period assets
%   e:       Labor efficiency shock
%   age:     Age of individual: young or old
% TOOLKIT NOTATION
% (d,aprime,a,z), where z = [e;age]

F = -inf;

r = alpha*K_to_L^(alpha-1)-delta;
w = (1-alpha)*K_to_L^alpha;

income = (w*e*d)*(age==1)+pen*(age==2)+r*a;

c = income+a-aprime; % Budget Constraint

if c>0
    % NOTE: 0<l<1 is already built into the grid
    % WARNING: this will not work if crra=1
    inside = (c^gamma)*((1-d)^(1-gamma));
    F = inside^(1-crra)/(1-crra);
end

end