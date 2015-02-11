function D = osl_reject_bad_epoch_trials(S)
% Removes trials that overlap with epochs marked as bad by OSLview or OPT
% D = osl_reject_bad_epoch_trials(S)
%
% S.D                - SPM MEG object filename
% S.bad_event_type   - string containing artefact to reject 
%                      (default 'artefact_OSL')
%
% Adam Baker 2015

if ~isfield(S,'bad_event_type')
    S.bad_epoch_type = 'artefact_OSL';
end

D = spm_eeg_load(S.D);

Events = D.events;

Badtrials = [D.badtrials,find(arrayfun(@(ev) any(strcmp({Events{ev}.type},S.bad_epoch_type)),1:D.ntrials))];

D = badtrials(D,Badtrials,1);

D.save

end

