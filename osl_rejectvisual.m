function D = osl_rejectvisual(D,time_range);
% D = osl_rejectvisual(S)
%
% Visual artefact rejection script
%
% Input
% - D : SPM MEEG object
% - time_range : Latency argument for spm_eeg_ft_artefact_visual
%
% RA 2017 
% MWW and LH

%% spit out some info about what has already been rejected in this file

if ~isempty(D.badchannels)
    disp('Following channels are already marked as bad:');
    lab = chanlabels(D);
    disp(lab(badchannels(D)));
else
    disp('No bad channels currently marked.');
end

if length(D.badtrials)>0
    disp([num2str(length(find(D.badtrials))) ' trials of ' num2str(D.ntrials) ' have already been marked as bad.']);
else
    disp('No bad trials currently identified.');
end

%% first, look for EOG artefacts

if ~isempty(find(ismember(D.chantype,'EOG')))
    disp('%%% 1. EOG %%%')
    
    S = [];
    S.D = D;
    S.method = 'summary';
    S.latency = time_range*1000;
    S.channel = D.chanlabels(find(ismember(D.chantype,'EOG')));
    S.metric = 'kurtosis';
    S.save = false;
    D = spm_eeg_ft_artefact_visual(S)
else
    warning('No EOG channel found.');
end

close all;

%% next, magnetometer artefacts

if ~isempty(find(ismember(D.chantype,'MEGMAG')))
    disp('%%% 2. MAGNETOMETERS %%%')
    
    S = [];
    S.D = D;
    S.method = 'summary';
    S.channel = D.chanlabels(find(ismember(D.chantype,'MEGMAG')));
    S.latency = time_range*1000;
    S.metric = 'var';
    S.save = false;

    D = spm_eeg_ft_artefact_visual(S)
else
    warning('No magnetometers found.');
end

close all;

%% next, gradiometer artefacts

if ~isempty(find(ismember(D.chantype,'MEGPLANAR')))
    disp('%%% 3. PLANAR GRADIOMETERS %%%')
    
    S = [];
    S.D = D;
    S.method = 'summary';
    S.channel = D.chanlabels(find(ismember(D.chantype,'MEGPLANAR')));
    S.latency = time_range*1000;
    S.metric = 'var';
    S.save = false;

    D = spm_eeg_ft_artefact_visual(S)
else
    warning('No planar gradiometers found.');
end

close all;

%% next, gradiometer artefacts

if ~isempty(find(ismember(D.chantype,'MEGGRAD')))
    disp('%%% 4. AXIAL GRADIOMETERS %%%')
    
    S = [];
    S.D = D;
    S.method = 'summary';
    S.channel = D.chanlabels(find(ismember(D.chantype,'MEGGRAD')));
    S.latency = time_range*1000;
    S.metric = 'var';
    S.save = false;

    D = spm_eeg_ft_artefact_visual(S)
else
    warning('No axial gradiometers found.');
end

close all;

%% next, EEG artefacts

if ~isempty(find(ismember(D.chantype,'EEG')))
    disp('%%% 5. EEG %%%')
    
    S = [];
    S.D = D;
    S.method = 'summary';
    S.channel = D.chanlabels(find(ismember(D.chantype,'EEG')));
    S.latency = time_range*1000;
    S.metric = 'var';
    S.save = false;

    D = spm_eeg_ft_artefact_visual(S);
else
    warning('No EEG channels found.');
end

close all;

%% spit out some info about what has now been rejected in this file

if ~isempty(D.badchannels)
    disp('Following channels are now marked as bad:');
    lab = chanlabels(D);
    disp(lab(badchannels(D)));
else
    disp('No bad channels marked.');
end

if length(D.badtrials)>0
    disp([num2str(length(D.badtrials)) ' trials of ' num2str(D.ntrials) ' are marked as bad.']);
else
    disp('No bad trials identified.');
end
