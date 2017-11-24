function D_epoched = osl_epoch(S)
% D = osl_epoch( S )
%
% Epochs the data using spm_eeg_epochs and additionally removes any trials
% that overlap with epochs previously marked as bad using OPT or OSLview
%
% S                 - struct containing arguments to pass to spm_eeg_epochs(S)
%
% S.mark_bad_trials - Use continuous bad samples (e.g. 'artefact_OSL' events)
%                     to mark affected trials as bad (default = True)
%
% Romesh Abeysuriya 2017
% Adam Baker 2015

if ischar(S.D)
    S.D = spm_eeg_load(S.D);
end

if ~(S.D.montage('getindex')==0)
    warning('Montage is being switched to zero for epoching')
    S.D = S.D.montage('switch',0);
end

fprintf(1,'** Saving MEEG to disk for epoching **\n');
S.D.save;
    
if ~isfield(S,'mark_bad_trials')
    S.mark_bad_trials = true;
end

D_epoched = spm_eeg_epochs(S);

if isfield(S,'epochinfo')
    D_epoched.epochinfo = S.epochinfo; % store epoch info inside D object
else
    D_epoched.epochinfo = S;
end

D_epoched.epochinfo.time_continuous = S.D.time;

if S.mark_bad_trials
    Badtrials = squeeze(~all(good_samples(D_epoched,D_epoched.indchantype('MEEG','GOOD')),2)); % This is True if a trial contains a bad sample in any imaging modality
    fprintf('%d of %d trials contained bad samples and were marked bad\n',sum(Badtrials),length(Badtrials));
    D_epoched = D_epoched.badtrials(find(Badtrials),1);
end

fprintf(1,'** Saving epoched MEEG to disk **\n');
D_epoched.save;

