function [data] = osl_baseline_correct(S)

% [data] = osl_baseline_correct(S)
%
% S.data as [nchannels x ntpts x ntrials x nfreq]
% S.time as a vector or tpts for 2nd dimension of S.data
% S.time_window as [start end] in secs
%
% Output data will be [nchannels x ntpts x ntrials x nfreq]

baseline_time_indices = find(S.time>S.time_window(1) & S.time<S.time_window(2));
baseline_amp = mean(S.data(:,baseline_time_indices,:,:),2);

data = S.data-repmat(baseline_amp,[1 length(S.time) 1 1]);

end

