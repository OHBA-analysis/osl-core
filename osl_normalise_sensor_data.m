function [Dnew pcadim normalisation] = osl_normalise_sensor_data(S)

% [D pcadim] = osl_normalise_sensor_data(S)

try D=S.D; catch, error('S.D must be specified'); end
try datatype=S.datatype; catch, error('S.datatype must be specified'); end % 'ctf' or 'neuromag' or 'eeg'

try do_plots=S.do_plots; catch, do_plots=false; end
try normalise_method=S.normalise_method; catch, normalise_method='mean_eig'; end
try conditions=S.conditions; catch, conditions=[]; conditions{1}='all'; end
try time_range=S.time_range; catch, time_range=[D.time(1) D.time(end)]; end
try force_pca_dim=S.force_pca_dim; catch, force_pca_dim=false; end

try 
    modalities=S.modalities; 
    S = ft_checkopt(S,'modalities','cell',{{'MEGMAG'},{'MEGPLANAR'},{'MEGMAG';'MEGPLANAR'}});      
catch
    switch datatype
        case 'neuromag'
            modalities={'MEGPLANAR', 'MEGMAG'}; 
        case 'ctf' 
            modalities={'MEG'}; 
        case 'eeg' 
            modalities={'EEG'}; 
        otherwise
            error('datatype not set properly. \n'); % should not be here
    end
    warning(['S.modalities not set. Will set to default']); 
end % modalities to include

switch datatype
    case 'neuromag'
        try pca_dim=S.pca_dim; catch, pca_dim=50; end
    case 'ctf' 
        try pca_dim=S.pca_dim; catch, pca_dim=260; end
    case 'eeg' 
        try pca_dim=S.pca_dim; catch, pca_dim=-1; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Establish time windows of interest

if isempty(time_range)
    time_range = [D.time(1) D.time(end)];
end

if strcmp(modalities{1},'EEG')
    modality_meeg='EEG';
else
    modality_meeg='MEGANY';
end

if(D.ntrials==1)
    goodsamples = good_samples(D,D.indchantype(modality_meeg));
else
    goodsamples = true(1,D.nsamples); 
end

samples_of_interest=zeros(1,D.nsamples);
for i=1:size(time_range,1)
    samples_of_interest(D.indsample(time_range(i, 1)):D.indsample(time_range(i, 2)))=1;
end

samples2use = samples_of_interest & goodsamples;
woi=[D.time(find(diff([0 samples2use])==1))' D.time(find(diff([samples2use 0])==-1))'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Establish trials

if strcmp(conditions{1},'all')
    trials = D.indtrial(D.condlist,'good');
else
    
    trials = D.indtrial(conditions,'good');
    if isempty(trials)
        error('No trials matched the selection, check the specified condition labels');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalise

if strcmp(normalise_method,'zscore')
    % Z-score samples within channels and trials
    [dirname, fname, ext] = fileparts(D.fullfile);
    Dnew = D.clone( fullfile(dirname,['M' fname ext]) );

    dat = D(:,:,:);
    for jj = 1:size(dat,3)
        dat(:,:,jj) = zscore(dat(:,:,jj),[],2);
    end
    Dnew(:,:,:) = dat;

    chanind=indchantype(Dnew,modality_meeg,'good');
    pcadim=length(chanind);
    normalisation=1;

elseif ~strcmp(normalise_method,'none')

    S2.D = D;
    S2.samples2use=samples2use;
    S2.trials=trials;
    S2.do_plots=do_plots;
    S2.modalities=modalities;
    S2.pca_dim=pca_dim;
    S2.normalise_method=normalise_method;
    S2.force_pca_dim=force_pca_dim;
    
    %[ Dnew pcadims pcadim ] = normalise_sensor_data( S2 );

    [ Dnew pcadim tmp norm_vec normalisation fig_handles fig_names] = normalise_sensor_data( S2 );

end

