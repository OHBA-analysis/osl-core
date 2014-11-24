function data = osl_convert_trialwise_to_fieldtrip(trlwise)
% GJW and MP

% convert to fieldtrip format
data            = [];
data.trial      = cell(1,size(trlwise.dat,2));
data.time       = cell(1,size(trlwise.dat,2));

% sort trials back into original order - NB, will be missing bad trials!!
[sortedtrials,sortedinds] = sort(trlwise.triallist);
for iTrial = 1:size(trlwise.dat,2)
    data.trial{iTrial} = permute(trlwise.dat(:,sortedinds(iTrial),:),[1,3,2]);
    data.time{iTrial}  = trlwise.times;
end

data.label   = trlwise.labels;
data.fsample = trlwise.fsample;
data.dimord  = 'chan_time';
data.triallist = sortedtrials;