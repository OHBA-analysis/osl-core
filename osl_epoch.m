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

good_meeg_chans=D_epoched.indchantype('MEEG','GOOD');
tmp=good_samples(D_epoched,good_meeg_chans);


if S.mark_bad_trials
    if false

        Badtrials = squeeze(~all(good_samples(D_epoched,D_epoched.indchantype('MEEG','GOOD')),2)); % This is True if a trial contains any bad sample in any good meg/eeg sensor

    else
        Badtrials=zeros(D_epoched.ntrials,1);
        bad_samps=~good_samples(S.D,S.D.indchantype('MEEG','GOOD'));

        % note that the trials in epochinfo can be different to those in
        % D_epoched, due to trials that run of the start or the end of the
        % continuous data
        trial_onsets_trl=D_epoched.epochinfo.trl(:, 1)./D_epoched.fsample; % in secs

        for ee=1:D_epoched.ntrials
            % find trial index in epochinfo.trl by matching trial onsets
            ind=find(trial_onsets_trl==D_epoched.trialonset(ee));
            
            if any(bad_samps(D_epoched.epochinfo.trl(ind,1):D_epoched.epochinfo.trl(ind,2)))
                Badtrials(ee)=1;
            end
        end

    end
    fprintf('%d of %d trials contained bad samples and were marked bad\n',sum(Badtrials),length(Badtrials));
    D_epoched = D_epoched.badtrials(find(Badtrials),1);

end

fprintf(1,'** Saving epoched MEEG to disk **\n');
D_epoched.save;

