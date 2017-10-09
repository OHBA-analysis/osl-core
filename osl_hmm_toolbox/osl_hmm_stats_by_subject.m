function [hmmstats, hmmstats_subj] = osl_hmm_stats_by_subject(hmm,varargin)
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


nsubs=max(hmm.subj_inds);

hmmstats_subj = struct;
for subnum = 1:nsubs,

    disp(['Computing for subj num ' num2str(subnum)]);
    for k = 1:hmm.K%length(unique(hmm.statepath))

        hmm_sub = hmm; 
        hmm_sub.statepath = hmm.statepath(hmm.subj_inds==subnum); 

        lifetimes = diff(logical2epoch(hmm_sub.statepath==k),[],2)./fs; lifetimes = lifetimes(lifetimes~=0);
        intervals = diff(logical2epoch(hmm_sub.statepath~=k),[],2)./fs; intervals = intervals(intervals~=0);

        hmmstats_subj.FractionalOccupancy(k,subnum) =  sum(hmm_sub.statepath==k) / length(hmm_sub.statepath);
        hmmstats_subj.nOccurrences(k,subnum) = length(lifetimes);
        hmmstats_subj.MeanLifeTime(k,subnum) = mean(lifetimes);
        hmmstats_subj.MeanIntervalLength(k,subnum) = mean(intervals);
        hmmstats_subj.LifeTimes{k,subnum} = lifetimes;
        hmmstats_subj.Intervals{k,subnum} = intervals;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ev = logical2epoch(l,t)

l = l(:)';

if nargin < 2
    t = 1:length(l);
end

onsets  = t(diff([0 l]) ==  1);
offsets = t(diff([l 0]) == -1);

ev = [onsets; offsets]';







