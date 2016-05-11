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

if isfield(S,'epochinfo')
    D_epoched.epochinfo = S.epochinfo; % store epoch info inside D object
else
    D_epoched.epochinfo = S;
end

D_epoched.epochinfo.time_continuous = D.time;

if ~isempty(S.bad_event_type)
    
    Badtrials = false(1,D_epoched.ntrials);
    for trl = 1:D_epoched.ntrials
        tpts = D_epoched.trialonset(trl) + [0:length(D_epoched.time)-1]/D.fsample;
        bs = badsamples(D,':',round(tpts*D.fsample),1);
        Badtrials(trl) = any(all(bs));
    end
    
    if 0
        %%debug stuff:
        trl=1
        tpts = D_epoched.trialonset(trl) + [0:length(D_epoched.time)-1]/D.fsample;
        inds=round(tpts*D.fsample);

        figure;
        plot(D.time(inds),D(1,inds,1));ho;
        plot(D.time(inds),D_epoched(1,:,trl),'r--');
    end
    
    D_epoched = badtrials(D_epoched,find(Badtrials),1);
    
end

D_epoched.save;

end
