function [cope, varcope, coape, dof, pe]=glm_fast_for_meg(y,x,pinvxtx,pinvx,cs,con)

% [cope, varcope, coape, dof]=glm_fast_for_meg(y,x,pinvxtx,pinvx,cs,con)
%
% Fast way to compute the GLM, y=x*beta + e, where the pinvx= pinv(x) and
% pinvxtx=pinv(x'*x) are precomputed and passed in
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

% y = demean(randn(100,1))+s'*0.1+s2'*2;
% t=[1:100];s=demean(sin(2*pi*t/20));s2=demean(sin(2*pi*t/30));x=[s; s2]';
% c=[1 1]';x=demean(x);[tstat, pe]=glm(y,x,c)
% t=[1:100];s=demean(sin(2*pi*t/20));s2=demean(sin(2*pi*t/30));x=[s+s2]'; c=[1]';x=demean(x);[tstat, pe]=glm(y,x,c)

if(size(y,1)~=size(x,1))
    disp('size(y):');
    size(y)
    disp('size(x):');
    size(x)
    
    error('size(y,1)~=size(x,1)');
end;

do_contrasts=1;
if(nargin<5),
    do_contrasts=0;
end;

pe=pinvx*y;

r=y-x*pe;

vr=diag(r'*r/(size(y,1)-size(x,2)));

if(do_contrasts)
    if(iscell(cs))
        if(con>0),
          c=cs{con}';
        else
          c=cell2mat(cs);
        end;
    else
        c=cs;
    end;

    varcope=diag(c'*pinvxtx*c*vr);
    cope=c'*pe;
    coape=c'*abs(pe);

else,
    cope=pe;
    varcope=pinvxtx*vr;
end;

dof=length(y)-size(x,2);