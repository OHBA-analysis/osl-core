%% osl_africa.m
%
% AfRICA - ArteFact Rejection using Independent Component Analysis
%
% Syntax: [bad_components, fig_handles, fig_names, fig_titles] = identify_artefactual_components_manual(S)
% S needs to contain be produced inside OSL_AFRICA.m
%
% OPTIONAL inputs to compute various metrics:
%   -  ident.do_mains:       compute metric based on mains frequency [0 or 1]
%   -  ident.do_kurt:        compute metric based on extreme kurtosis [0 or 1]
%   -  ident.do_cardiac:     compute metric based on similarity to cardiac cycle [0 or 1]
%   -  ident.artefact_chans: compute metric based on correlation with external channels [e.g. {'ECG','EOG'}]
%   -  ident.do_reject:      only compute & save metrics & topos (do not launch GUI)
% Output:
%   - bad_components: A list of the components identified as bad.
% 
% Adam Baker 2014
%
% adam.baker@ohba.ox.ac.uk
%
%
% TODO: set default processing for ECG, EOG (blink), EOG (saccade)
%


function [bad_components, fig_handles, fig_names, fig_titles] = identify_artefactual_components_manual(S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set-Up

% First check if metrics have already been computed:
try 
    metrics = S.ica_res.metrics;
catch
    metrics = [];
end


if ~isfield(S,'ident'); 
    error('Need to specify ident func variables');
end


if ~isfield(S.ident,'do_reject') || S.ident.do_reject == 1;
    do_reject = 1;
else
    do_reject = 0;
end


if isfield(S.ident,'do_mains') && S.ident.do_mains == 1;
    do_mains = 1;
    if ~isfield(S.ident,'mains_frequency')
        S.ident.mains_frequency = 50;
    end
else
    do_mains = 0;
end
   
if isfield(S.ident,'do_kurt') && S.ident.do_kurt == 1
    do_kurt = 1;
else
    do_kurt = 0;
end

if isfield(S.ident,'artefact_chans')
    try
       S.ident.artefact_chans = cellstr(S.ident.artefact_chans);
    catch
        error('S.ident.artefact_chans should be strings');             
    end    
else
    S.ident.artefact_chans = [];
end


if isfield(S.ident,'do_cardiac') && S.ident.do_cardiac
    if isfield(S.ident,'artefact_chans') && any(strcmpi(S.ident.artefact_chans,'ECG'))
        do_cardiac_autocorrelation = 0;
    else
        do_cardiac_autocorrelation = 1;
    end
else
    do_cardiac_autocorrelation = 0;
end


D = spm_eeg_load(S.fname);

fig_handles = [];
fig_names   = [];
fig_titles  = [];

%% Select only good data for classification

samples_of_interest = false(1,size(S.ica_res.tc,2));

% Remove bad trials/any trial structure
c=1;
for i=1:D.ntrials
    if ~ismember(i,D.badtrials)
        samples_of_interest(:,c:c+D.nsamples-1)=true;
    end
    c=c+D.nsamples;
end

if strcmp(S.modality,'EEG')   % changed by DM
    chan_inds=setdiff(find(any([strcmp(D.chantype,'EEG')],1)),D.badchannels);
    map_inds(find(any([strcmp(D.chantype,'EEG')],1))) = 1:numel(find(any([strcmp(D.chantype,'EEG')],1)));
else
    chan_inds=setdiff(find(any([strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1)),D.badchannels);
    map_inds(find(any([strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1))) = 1:numel(find(any([strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1)));
end

% Remove bad segments
if D.ntrials==1;
    t = D.time;
    badsections = false(1,D.nsamples);
    Events = D.events;
    if ~isempty(Events)
        Events = Events(strcmp({Events.type},'BadEpoch'));
        for ev = 1:numel(Events)
            badsections = badsections | t >= Events(ev).time & t < (Events(ev).time+Events(ev).duration);
        end
    end
    samples_of_interest(1,badsections)=false;
end

%%

sm = S.ica_res.sm;
tc = S.ica_res.tc;

sm = sm(map_inds(chan_inds),:);

num_ics = S.ica_res.ica_params.num_ics;

tc_fft = fft(demean(tc(:,samples_of_interest),2),[],2);

abs_ft  = abs(tc_fft);
freq_ax = 0:(D.fsample/2)/(floor(length(abs_ft)/2)-1):(D.fsample/2);

minhz = 0.1;
minhz_ind = min(find(freq_ax>minhz));

spec = abs_ft(1:num_ics,1:floor(length(abs_ft)/2));
spec(:,1:minhz_ind) = 0;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                           Mains frequency                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if do_mains && ~isfield(metrics,'mains')
    
    spec_n = normalise(spec,2);
    
    inds = freq_ax > S.ident.mains_frequency - 1 & freq_ax < S.ident.mains_frequency + 1;
    metrics.mains.value = max(spec_n(:,inds),[],2);
    %metrics.mains.label = 'mains';
elseif isfield(S.ident,'do_mains') && ~S.ident.do_mains && isfield(metrics,'mains')
    metrics = rmfield(metrics,'mains');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                               Kurtosis                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if do_kurt && ~isfield(metrics,'kurtosis')

    kurt = kurtosis(tc(:,samples_of_interest),[],2);
    kurt_t = abs(demean(boxcox1(kurt))); % make distribution more normal to better balance low/high kurtosis
    
    metrics.kurtosis.value = kurt_t(:);
    %metrics.kurtosis.label = 'kurtosis';
elseif isfield(S.ident,'do_kurt') && ~S.ident.do_kurt && isfield(metrics,'kurtosis')
    metrics = rmfield(metrics,'kurtosis');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                               Cardiac                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if do_cardiac_autocorrelation && ~isfield(metrics,'cardiac_autocorrelation')
    
  disp('No ECG artefact channels specified - detecting cardiac components without ECG')
  
  bpm_range = [50 100]./60;
  maxlags = 2*D.fsample;

  [~,lags] = xcorr(tc(1,samples_of_interest),maxlags,'coeff');
  lags = lags./D.fsample;
  lags_range = lags>bpm_range(1) & lags<bpm_range(2);
  
  ac_max = zeros(1,num_ics);
  for ic = 1:num_ics
    tc_bp = bandpass(tc(ic,samples_of_interest),[0 48],D.fsample);
    ac = xcorr(tc_bp,maxlags,'coeff');
    ac_max(ic) = max(ac(lags_range));
  end
  
  
  metrics.cardiac_autocorrelation.value = ac_max(:);
  %metrics.cardiac_autocorrelation.label = 'cardiac_autocorrelation';
  
elseif isfield(S.ident,'do_cardiac') && ~S.ident.do_cardiac && isfield(metrics,'cardiac_autocorrelation')
    metrics = rmfield(metrics,'cardiac_autocorrelation');
  
end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            Detect artefact channels related components              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(S.ident.artefact_chans) && ~isempty(S.ident.artefact_chans)
    
    for artefact_chantype = unique(S.ident.artefact_chans)
        if ~isfield(metrics,artefact_chantype)
                   
            artefact_data = D(find(strcmp(D.chantype,artefact_chantype)),find(samples_of_interest));
            artefact_data = reshape(artefact_data,size(artefact_data,1),[]);

            % Bandpass filter
            for ac = 1:size(artefact_data,1)
                artefact_data(ac,:) = bandpass(artefact_data(ac,:),[0.1 48],D.fsample);
            end
            
            ac_corr = zeros(1,num_ics);
            for ic = 1:num_ics
                tc_bp = bandpass(tc(ic,samples_of_interest),[0.1 48],D.fsample);
                %ac_corr(ic) = sum(abs(corr(tc_bp',artefact_data')));
                for ac = 1:size(artefact_data,1)
                    ac_corr(ic) = ac_corr(ic) + abs(nanmedian(osl_movcorr(tc_bp',artefact_data(ac,:)',D.fsample*10,0)));
                end                    
            end
            ac_corr = ac_corr ./size(artefact_data,1);

            
            metrics.(char(artefact_chantype)).('value') = ac_corr(:);
        end
    end
    
end


% Precompute topographies
try
    topos = S.ica_res.topos; 
catch
    topos = [];
    h = figure('visible','off');
    modalities = unique(D.chantype(find(strncmpi(S.modality,D.chantype,3))));
    for m = 1:numel(modalities)
        disp(['Precomputing sensor topographies for modality ' modalities{m}]);
        topos_m = component_topoplot(D,sm,modalities(m));
        topos = [topos handle2struct(topos_m)];
    end
    
    close(h)
end



% Save topos & metrics
if isfield(S,'ica_file')
    S.ica_res.topos   = topos;
    S.ica_res.metrics = metrics;
    save(S.ica_file,'S');
    msg = sprintf('\n%s%s\n%','Saving metrics and topos to ', S.ica_file);
    fprintf(msg);
end

% Add variance as default metric;
metrics.variance.value = var(sm);
metrics = orderfields(metrics,[find(strcmp(fieldnames(metrics),'variance'));find(~strcmp(fieldnames(metrics),'variance'))]);


tc(:,~samples_of_interest) = nan;


if do_reject
    if isfield(S.ica_res,'bad_components')
        bad_components = S.ica_res.bad_components;
    else
        bad_components = [];
    end
    bad_components = identify_artefactual_components_manual_gui(D,tc,topos,metrics,bad_components);
else
    bad_components = [];
end

end





% function to produce component topoplot
function topos = component_topoplot(D,comp,modality)
cfg  = [];
data = [];

global OSLDIR

if strcmp(modality,'EEG')   % changed by DM
    chan_inds=setdiff(find(any([strcmp(D.chantype,'EEG')],1)),D.badchannels);
    map_inds(find(any([strcmp(D.chantype,'EEG')],1))) = 1:numel(find(any([strcmp(D.chantype,'EEG')],1)));
else
    chan_inds=setdiff(find(any([strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1)),D.badchannels);
    map_inds(find(any([strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1))) = 1:numel(find(any([strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1)));
end

comp_full(map_inds(chan_inds),:)=comp;
comp_full(D.badchannels) = 0;
comp = comp_full;

if (strcmp(modality,'MEGMAG')),
    cfg.channel     = {'MEGMAG'};
    cfg.layout      = [OSLDIR '/layouts/neuromag306mag.lay'];
    chanind = strmatch(cfg.channel, D.chantype);
    comp2view = zeros(102,size(comp,2));
    for ind=1:length(chanind),
        indp=(ind-1)*3+3;
        comp2view(ind,:)=comp(indp,:);
    end
    
elseif (strcmp(modality,'MEGPLANAR')),
    cfg.channel     = {'MEGMAG'};
    cfg.layout      = [OSLDIR '/layouts/neuromag306mag.lay'];
    chanind = strmatch(cfg.channel, D.chantype);
    comp2view = zeros(102,size(comp,2));
    for ind=1:length(chanind),
        indp=(ind-1)*3+1;
        comp2view(ind,:)=comp(indp)+comp(indp+1,:); % is this okay??
    end
    
    
elseif (strcmp(modality,'MEGGRAD')),
    cfg.channel     = {'MEG'};
    cfg.layout      = [OSLDIR '/layouts/CTF275.lay'];
    chanind = strmatch(cfg.channel, D.chantype);
    comp2view=comp;
    
    
    % elseif (strcmp(modality,'EEG')),
    %     cfg.channel     = {'EEG'};
    %     A=sensors(D,'EEG');
    %     pnt=A.chanpos;
    %     elec = [];
    %     elec.pnt = pnt(chan_inds,:);
    %     elec.label = A.label(chan_inds);
    %     cfg.elec  = elec;
    %     %cfg.layout = ft_prepare_layout(cfg);
    %     cfg.channel={'all'};
    %
    %     comp2view=comp;
    %     data.dimord='chan_comp';
    %     data.topo=comp2view;
    %     data.topolabel = A.label(chan_inds);
    %     data.label = A.label(chan_inds);
    %     data.time = {1};
else
    
    error('Unsupported modality');
end

data.dimord = 'chan_comp';
data.topo   = comp2view;
data.topolabel = D.chanlabels(chanind);
data.time = {1};

cfg.component = 1:setdiff(size(comp),length(chan_inds));
cfg.interactive = 'no';
cfg.comment     = 'no';
cfg.title     = modality{:};

cfg.layout = ft_prepare_layout(cfg);

[~] = evalc('ft_topoplotIC(cfg,data);');
topos = get(gcf,'children');
end

