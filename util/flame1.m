function [gam, beta, covgam] = flame1(y,z,S,do_fixed_effects)

% [gam, beta, covgam] = flame1(y,z,S)
%
% solves GLM y=z*gam+e where e~N(0, beta+diag(S)) using FLAME1

if(nargin<4)
    do_fixed_effects=0;
end;

if(do_fixed_effects)
    beta=0;
    
    iU=diag(1./S);
else
    logbeta = solveforlogbetanew(1e-10,y,z,S);

    beta = exp(logbeta);
    
    iU=diag(1./(S+exp(logbeta)));
end;

covgam = pinv(z'*iU*z);

gam = covgam*z'*iU*y;
