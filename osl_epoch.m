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

D = spm_eeg_load(S.D);

if ~isfield(S,'bad_event_type')
    S.bad_event_type = 'artefact_OSL';
end

D_epoched = spm_eeg_epochs(S);

D_epoched.epochinfo=S.epochinfo; % store epoch info inside D object
D_epoched.epochinfo.time_continuous=S.D.time;

if ~isempty(S.bad_event_type)
    
    Badtrials = false(1,D_epoched.ntrials);
    for trl = 1:D_epoched.ntrials
        tpts = D_epoched.trialonset(trl) + D_epoched.time;
        bs = badsamples(D,':',round(tpts*D.fsample),1);
        Badtrials(trl) = any(all(bs));
    end
    
    D_epoched = badtrials(D_epoched,find(Badtrials),1);
    
end

D_epoched.save;

end
