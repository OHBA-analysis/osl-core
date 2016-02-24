function [ent] = dirichlet_entropy(alpha)

% [ent] = dirichlet_entropy(alpha)
%
% comupte entropy for a dirichlet dist.
% MWW 2016
    
NK=length(alpha);
alpha_0=sum(alpha);
B=prod(gamma(alpha))/gamma(sum(alpha));
ent= log(B) + (alpha_0-NK)*digamma(alpha_0);
ent= ent - sum((alpha-1).*digamma(alpha));

