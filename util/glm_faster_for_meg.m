function [cope, varcope, coape, dof, pe] = glm_faster_for_meg(y,x,cs,con)

% [cope, varcope, coape, dof]=glm_faster_for_meg(y,x,cs,con)
%
% Fast way to compute the GLM, y=x*beta + e, when pseudo inverse of x has
% yet to be computed. 
% cs is the list of contrasts
% y is the data (num_observations x 1)
% x is design matrix (num_observations x num_regressors)
% con is 0 to compute all contrasts, or an index of the contrast to be
% computed in cs
%
% cope is the contrast of the parameter estimates = c'*beta
% varcope is the variance of the cope
% coape is the contrast of the absolute values of the parameter estimates = c'*abs(beta)
% dof is the degrees of freedom
%
% This code is about 3 times faster than 
%   pinvx = pinv(x); pinvxtx = pinv(x'*x); 
%   glm_fast_for_meg(y,x,pinvxtx,pinvx,cs,con);

% GC, based on MW, 2015

% input parsing
assert(size(y,1) == size(x,1), ...
       [mfilename ':SizeError'], ...
       'Number of rows in x and y must be the same. \n');
   
if nargin < 3,
    do_contrasts = 0;
else
    do_contrasts = 1;
end

% main code
xtxC = chol(x'*x);
pe   = xtxC \ (xtxC' \ (x'*y));
r    = y - x*pe;
vr   = diag(r' * r / (size(y,1)-size(x,2)));

if do_contrasts
    if iscell(cs)
        if con > 0,
          c = cs{con}';
        else
          c = cell2mat(cs);
        end%if
    else
        c = cs;
    end%if

    varcope = diag(c' * (xtxC \ (xtxC' \ (c * vr))));
    cope    = c' * pe;
    coape   = c' * abs(pe);

else
    cope    = pe;
    R       = inv(xtxC); % fast and accurate for upper diagonal
    varcope = R * R' * vr;
end%if

dof = length(y)-size(x,2);

end%glm