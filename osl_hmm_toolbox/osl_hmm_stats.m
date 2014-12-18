function hmmstats = osl_hmm_stats(hmm,varargin)
% Compute summary statistics from HMM statepath. 
%
% HMMSTATS = OSL_HMM_STATS(HMM) returns a number of statistics from the
% HMM state path (fractional occupancy, number of occurences, mean life
% time, mean interval length, life times, interval lengths). These 
% statistics are returned in seconds if hmm.fsample exists, otherwise they 
% are returned in samples.
%
% HMMSTATS = OSL_HMM_STATS(HMM,'do_plots') also plots these statistics
%
% Adam Baker 2013

if ~isfield(hmm,'fsample') || isempty(hmm.fsample)
    fs = 1;
else
    fs = hmm.fsample;
end

hmmstats = struct;
for k = 1:hmm.K%length(unique(hmm.statepath))
    
    lifetimes = diff(logical2epoch(hmm.statepath==k),[],2)./fs; lifetimes = lifetimes(lifetimes~=0);
    intervals = diff(logical2epoch(hmm.statepath~=k),[],2)./fs; intervals = intervals(intervals~=0);
    
    hmmstats(k).FractionalOccupancy =  sum(hmm.statepath==k) / length(hmm.statepath);
    hmmstats(k).nOccurrences = length(lifetimes);
    hmmstats(k).MeanLifeTime = mean(lifetimes);
    hmmstats(k).MeanIntervalLength = mean(intervals);
    hmmstats(k).LifeTimes = lifetimes;
    hmmstats(k).Intervals = intervals;
    
end


if any(strcmpi(varargin,'do_plots'))
    
    fnames = fieldnames(hmmstats);
    for f = 1:length(fnames)
        
        sp = ABsubplot(length(fnames),f);     

        if numel(cat(1,hmmstats.(fnames{f})))==hmm.K
            bar(sp,[hmmstats.(fnames{f})])    
        else
            vp_data = vertcat(hmmstats.(fnames{f}));
            vp_grps = cell2mat(arrayfun(@(x) x*ones(1,numel(hmmstats(x).(fnames{f}))),1:hmm.K,'uniformoutput',0));
            ABviolinplot(sp,vp_data,vp_grps,50)
        end
        
        xlabel(sp,'State #');
        title(sp,fnames(f));
                
    end
    
end

end


function ev = logical2epoch(l,t)

l = l(:)';

if nargin < 2
    t = 1:length(l);
end

onsets  = t(diff([0 l]) ==  1);
offsets = t(diff([l 0]) == -1);

ev = [onsets; offsets]';

end

function sp = ABsubplot(N,n)
i = ceil(sqrt(N));
j = ceil(N/i);
if nargin > 1
    sp = subplot(i,j,n);
end
end




function ABviolinplot(h,data,groups,nbins)

cla(h), hold(h,'on')
groups = groups(:);

for n = unique(groups)'
  
  Y = data(groups==n);
  
  [h,bins] = hist(Y,nbins);
  h = 0.4*h./max(h);
  h = smooth(h,5);
  h = h(:)';
  patch([n-h n+h(end:-1:1)],[bins bins(end:-1:1)],[0.4 0.4 1],'EdgeColor','k')
  
  plot(n,mean(Y),'xk','markersize',10)
  
end

end



