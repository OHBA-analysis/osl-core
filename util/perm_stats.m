function [corrp tstats] = perm_stats( c, nP )

% [corrp tstats] = perm_stats( c, nP )
%
% Only works for a group average at the moment to do permutation tests 
% (using sign flipping).
%
% c is data where:
% nS = size(c,1); %number of 'subjects' (recording sites)
% nF = size(c,2); %number of observations
%
%
% nP is num of permutations


if nargin<3
  nP = 500;
end

nS = size(c,1); %number of 'subjects' (recording sites)
nF = size(c,2); %number of obs

%build up a null distribution of the maximum t-stat for each
%permuatation i:
nulldist=zeros(nP,1);

for i = 1:nP

    dm = ((rand(nS,1)>0.5)-0.5)*2;
  
    [cg,vg,tg] = ols(c,dm,1);

    
    nulldist(i)=max(tg);
  
end

%run a one sample t-test on data and compare cluster size to
%permutation data
dm = ones(nS,1);
[cg,vg,tg] = ols(c,dm,1);

for kk = 1:nF %loop over obs
  cp(kk) = mean(tg(kk)>nulldist);
%  if cp(kk)~=0,    keyboard;end
end

corrp = cp;
tstats = tg;

end

