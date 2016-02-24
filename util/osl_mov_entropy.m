function [ents,ts,ents_empirical] = osl_mov_entropy(data,t,winsize,overlap)

% [ents,ts,ents_empirical] = osl_mov_entropy(data,t,winsize,overlap)
% -----------------------------------------------------------------
% data    - data (num_time_points x num_variables)
% t       - vector of sample times
% winsize - window to average over (samples)
% overlap - amount of segment overlapping (e.g. 0.75 = 75% overlap)
% -----------------------------------------------------------------
% MWW 2016
    
if nargin < 4
  overlap = 0.75;
end

if isempty(t) 
  t = 1:size(data,1);
end

NK=size(data,2);

tbf      = buffer(t,winsize,round(winsize*overlap),'nodelay');
for kk=1:NK
    tmp       = buffer(data(:,kk),winsize,round(winsize*overlap),'nodelay');

    if kk==1
        bf=zeros([NK,size(tmp)]);
    end
    
    bf(kk,:,:)=tmp;
end

ents=zeros(size(bf,3)-1,1);
ents_empirical=zeros(size(bf,3)-1,1);
ts=zeros(size(bf,3)-1,1);

for ww=1:(size(bf,3)-1)
    
    alpha=mean(bf(:,:,ww),2);
    alpha=alpha/sum(alpha);
    alpha=max(alpha,0.01);
    
    % entropy = disorder
    % compute entropy for a dirichlet distribution
    ents(ww) = dirichlet_entropy(alpha);
    
    % compute binned empirical entropy
    ents_empirical(ww) = - sum(alpha .* log(alpha));

    ts(ww)=mean(tbf(:,ww),1);

end


