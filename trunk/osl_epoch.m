function D_epoched = osl_epoch(S)
% D = osl_epoch( S )
%
% Epochs the data using spm_eeg_epochs and additionally removes any trials
% that overlap with epochs previously marked as bad using OPT or OSLview
%
% S                 - struct containing arguments to pass to spm_eeg_epochs(S)
%
% S.bad_event_type  - string containing event type to remove 
%                     (default artefact_OSL')
%                     OR leave empty to keep all trials
%
% Adam Baker 2015


if ~isfield(S,'bad_event_type')
    S.bad_event_type = 'artefact_OSL';
else
    S.bad_event_type = [];
end

D_epoched = spm_eeg_epochs(S);

if ~isempty(S.bad_event_type)
    
    Badtrials = false(1,D_epoched.ntrials);
    for trl = 1:D_epoched.ntrials
        tpts = D_epoched.trialonset(trl) + D_epoched.time;
        bs = badsamples(D,':',round(tpts*D.fsample),1);
        Badtrials(trl) = any(all(bs));
    end
        
end

D_epoched.save;

end

