%% osl_africa.m
%
% AfRICA - ArteFact Rejection using Independent Component Analysis
%
% Syntax: fname_out=osl_africa(S)
% S needs to contain:
%   -  fname: full name of input SPM object (inc. path if not in present
%      working directory) e.g. S.fname = '/home/data/spm8_mydata'.
%    - ica_file: user specified file name for the the ICA decomposition. If
%      ica_file does not exist then AfRICA will save the ICA decomposition
%      here. If it does then AfRICA will attempt to read the existing
%      result in.
%
% OPTIONAL inputs:
%   -  do_plots: set to 1 to output summary plots of artefact components
%      (in identify_artefactual_components_auto.m) DEFAULT = 0.
%   -  ica_params: a matlab structure containing the following fields:
%           - num_ics: Number of indpendent components sought: DEFAULT =
%           150.
%           - last_eig: order of PCA prewhitening: DEFAULT = num_ics.
%           - nonlinearity: DEFAULT = 'tanh'
%           - stabilisation: DEFAULT = 'on'
%           - max_iterations: DEFAULT = 1000
%           consult fastica documentation for more information on ica_params.
%    - used_maxfilter: set to 1 if using Maxfiltered data. DEFAULT = 0;

%    - ident.func: function handle for user-specified function for selecting
%      bad components. Default @IDENTIFY_ARTEFACTUAL_COMPONENTS_MANUAL. Use
%      this as a template ofr custom functions.
%    - ident: contains all settings for ident.func
% outputs
%   - fname_out: the name of the output SPM object. The prefix 'A' is
%     attached to any objects that have been cleaned up with AfRICA. 
%     e.g. fname_out = '/home/data/Aspm8_mydata.dat'
%
%
% AfRICA ignores bad trials defined by D.badtrials and bad periods defined by
% 'BadEpoch' events using OSLVIEW. Do not reset these periods to "good"
% after AfRICA.
%
% Written by Henry Luckhoo and Adam Baker
% Maintained by adam.baker@ohba.ox.ac.uk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_out, fig_handles, fig_names, fig_titles, S]=osl_africa(S)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    D=spm_eeg_load(S.fname);
    S.D=D;
catch
    error('Can not find SPM file from S.fname')
end

if ~isfield(S,'ica_file') || isempty(S.ica_file)
    error('No file specified in S.ica_file for saving AFRICA results')
end


[dir,nam]=fileparts(S.fname);
fname_out=[dir '/A' nam '.dat'];

if isfield(S,'ica_file')
  [pathstr, name, ext] = fileparts(S.ica_file);
    if isempty(pathstr)
        pathstr = dir;
    end
    if isempty(ext)
        ext = '.mat';
    end
    S.ica_file = fullfile(pathstr,[name ext]);      
end

if not(isfield(S,'modality'))   % added by DM
    S.modality='MEG';
end

if not(isfield(S,'do_plots'))   % added by DM
    S.do_plots=0;
end

if ~isfield(S,'to_do');
    S.to_do = [1 1 1]; 
end

if ~isfield(S,'ident')
    S.ident = [];
end

if ~isfield(S.ident,'func')
    S.ident.func = @identify_artefactual_components_manual;
end


if isfield(S,'logfile') && S.logfile == 1
  logdir = [dir '/Africa_logs/'];
  if ~isdir(logdir), mkdir(logdir); end
  if exist([logdir nam '_log.txt'],'file')
    unix(['rm ' logdir nam '_log.txt']);
  end
  	diary([logdir nam '_log.txt']);
end


if ~isfield(S,'used_maxfilter'); 
    S.used_maxfilter = 0;
end

fig_handles = []; 
fig_names   = [];
fig_titles  = [];


%%%%%%%%%%%%%%%%%%%%%%%%% PERFORM SENSOR SPACE ICA %%%%%%%%%%%%%%%%%%%%%%%%

if S.to_do(1)
    S.ica_res = perform_sensorspace_ica(S);
    if isfield(S,'ica_file')
        save(S.ica_file,'S');
        msg = sprintf('\n%s%s\n%','Saving ICA results to ', S.ica_file);
        fprintf(msg);
    else
        msg = sprintf('\n%s\n%','Results not being saved.');
        fprintf(msg);
    end
elseif any(S.to_do) % Load from file
    if isfield(S,'ica_file') && exist(S.ica_file,'file')==2
        msg = sprintf('\n%s%s\n%','Loading previous ICA results from ', S.ica_file);
        fprintf(msg);
        icafile = load(S.ica_file);
        if isfield(S,'ident');    icafile.S.ident    = S.ident;    end
        if isfield(S,'do_plots'); icafile.S.do_plots = S.do_plots; end
        if isfield(S,'to_do');    icafile.S.to_do    = S.to_do;    end
        S = icafile.S; clear icafile
    end 
end


% AB comment: what's with the inconsistency between ica_file and ica_res?
% Surely we only need ica_file (containing S.ica_res)?! I'm going to sort 
% this out, if it breaks anything then too bad!








%%%%%%%%%%%%%%%%%%%%%%%%%% CLASSIFY BAD COMPONENTS %%%%%%%%%%%%%%%%%%%%%%%%

if S.to_do(2)
    % Allow user-specified functions to be used to identify bad components.

    
    if isfield(S,'ica_res')
        if(S.do_plots)
            [bad_components, fig_handles, fig_names, fig_titles] = feval(S.ident.func,S);
        else
            bad_components = feval(S.ident.func,S);
        end
        
        if isfield(S,'ica_file')
            icafile = load(S.ica_file);
            S.ica_res = icafile.S.ica_res; clear icafile
        end
        
        S.ica_res.bad_components = bad_components;
        
        msg = sprintf('\n%s%s\n%','Saving bad component selection to ', S.ica_file);
        fprintf(msg);
        save(S.ica_file,'S');

    else
        error('ICA decomposition needs to be run. Set S.to_do(1) = 1 and try again');
    end
end



%%%%%%%%%%%%%%%%%% REMOVE BAD COMPONENTS FROM THE DATA %%%%%%%%%%%%%%%%%%%%

if S.to_do(3)
    if ~isfield(S,'ica_res')
        try
            load(S.ica_file);
        catch
            error('Unable to load ICA decomposition from disk. Try rerunning all AfRICA stages');
        end
    end
    fname_out = remove_bad_components(S);
else
    fname_out = [];
end

if isfield(S,'logfile') && S.logfile == 1,
    diary off;
end

end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUB-FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% perform_sensorspace_ica - subfunction to decompose meg data using fastICA

function ica_res = perform_sensorspace_ica(S)

D = spm_eeg_load(S.fname);

% 1.) Select good channels
if strcmp(S.modality,'EEG')
    chan_inds=setdiff(find(any(strcmp(D.chantype,'EEG'),1)),D.badchannels);
else
    chan_inds=setdiff(find(any([strcmp(D.chantype,'MEGMAG'); strcmp(D.chantype,'MEGPLANAR'); strcmp(D.chantype,'MEGGRAD')],1)),D.badchannels);
end

meg_dat=D(chan_inds,:,:);

% 2.) Remove bad trials/any trial structure
meg_dat = meg_dat(:,:,setdiff(1:D.ntrials, D.badtrials));
meg_dat = reshape(meg_dat,size(meg_dat,1),[]);

% 3.) Remove bad segments
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
    meg_dat(:,badsections)=[];
end

%%%%%%%%%%%%%%%%%%%% APPLY MAXFILTER SPECIFIC SETTINGS %%%%%%%%%%%%%%%%%%%%

if isfield(S,'used_maxfilter') && S.used_maxfilter
    num_ics_default = 62;
    mag_cutoff      = 62;
    plan_cutoff     = 62;
else
    num_ics_default = 150;
    mag_cutoff  = sum(strcmp(D.chantype(chan_inds),'MEGMAG'))    - 5;
    plan_cutoff = sum(strcmp(D.chantype(chan_inds),'MEGPLANAR')) - 5;
end

%%%%%%%%%%%%%%%%%%%%%%%%% SET FASTICA PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(S,'ica_params')
    if isfield(S.ica_params,'num_ics');         num_ics       = S.ica_params.num_ics;         else num_ics       = num_ics_default; end
    if isfield(S.ica_params,'last_eig');        last_eig      = S.ica_params.last_eig;        else last_eig      = num_ics_default; end
    if isfield(S.ica_params,'nonlinearity');    nonlinearity  = S.ica_params.nonlinearity;    else nonlinearity  = 'tanh';          end
    if isfield(S.ica_params,'approach');        ica_approach  = S.ica_params.approach;        else ica_approach  = 'symm';          end
    if isfield(S.ica_params,'stabilization');   stabilization = S.ica_params.stabilization;   else stabilization = 'on';            end
    if isfield(S.ica_params,'max_iterations');  max_iter      = S.ica_params.max_iterations;  else max_iter      = 1000;            end
else
    num_ics      = num_ics_default;
    last_eig     = num_ics_default;
    nonlinearity = 'tanh';
    ica_approach = 'symm';
    stabilization = 'on';
    max_iter     = 1000;
end

num_ics  = min(num_ics,size(meg_dat,1));  % added by DM
last_eig = min(last_eig,size(meg_dat,1)); % added by DM


%%%%%%%%%%%%%%%%%%%%  MINIMUM EIGENVALUE NORMALISATION %%%%%%%%%%%%%%%%%%%%

if strcmp(S.modality,'EEG')  % added by DM
    norm_vec = max(abs(meg_dat(:)))/1000*ones(size(meg_dat,1),1);
else
    norm_vec = ones(numel(chan_inds),1);
    if any(strcmp(D.chantype,'MEGMAG')) && any(strcmp(D.chantype,'MEGPLANAR'))
        mag_min_eig = svd(cov(meg_dat(strcmp(D.chantype(chan_inds),'MEGMAG'),:)')); 
        mag_min_eig = mean(mag_min_eig(mag_cutoff-2:mag_cutoff));
        
        plan_min_eig = svd(cov(meg_dat(strcmp(D.chantype(chan_inds),'MEGPLANAR'),:)')); 
        plan_min_eig = mean(plan_min_eig(plan_cutoff-2:plan_cutoff));
        
        norm_vec(strcmp(D.chantype(chan_inds),'MEGMAG'))    = mag_min_eig; 
        norm_vec(strcmp(D.chantype(chan_inds),'MEGPLANAR')) = plan_min_eig;
    else
        norm_vec = norm_vec*min(svd(cov(meg_dat(strcmp(D.chantype(chan_inds),'MEGGRAD'),:)')));
    end
    norm_vec = sqrt(norm_vec);
    
end

dat_in = meg_dat ./ repmat(norm_vec,1,size(meg_dat,2));

if S.do_plots
    figure; 
    semilogy(svd(cov(meg_dat')));
    ho;
    semilogy(svd(cov(dat_in')),'r--');
    title('Raw and normalised eigen spectra'); legend('Raw', 'Normalised');
end

%%%%%%%%%%%%%%%%%% AUTOMATIC DIMENSIONALITY ESTIMATION %%%%%%%%%%%%%%%%%%%%
if 0 == num_ics,
    num_ics = spm_pca_order(dat_in.');
end%if
if 0 == last_eig,
    last_eig = num_ics;
end%if

%%%%%%%%%%%%%%%%%%%%%%%%%%%% ICA DECOMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tc,sm,~] = fastica(dat_in,                            ...
                    'g',                nonlinearity,  ...
                    'lastEig',          last_eig,      ...
                    'numOfIC',          num_ics,       ...
                    'approach',         ica_approach,  ...
                    'stabilization',    stabilization, ...
                    'maxNumIterations', max_iter); % changed by DM

if num_ics ~= size(tc,1)
    msg = sprintf('\n%s%d%s%d%s\n%','Data dimensionality insufficient to support ', num_ics, ' components. Number of components has been reduced to ', size(tc,1), '.');
    fprintf(msg);
    num_ics = size(tc,1);
end

ica_res = [];
ica_res.ica_params.num_ics      = num_ics;
ica_res.ica_params.last_eig     = last_eig;
ica_res.ica_params.nonlinearity = nonlinearity;
ica_res.D                       = S.fname;
ica_res.tc                      = tc;
ica_res.sm                      = sm .* repmat(norm_vec,1,num_ics);

%%%%%%%%%%%%%%%%%%%%%%%%% ESTIMATE MISSING CHANNELS %%%%%%%%%%%%%%%%%%%%%%% 

if strcmp(S.modality,'EEG')  % added by DM
    sm_full = zeros(numel(find(any(strcmp(D.chantype,'EEG'),1))),ica_res.ica_params.num_ics);
    map_inds(find(any(strcmp(D.chantype,'EEG'),1))) = 1:numel(find(any(strcmp(D.chantype,'EEG'),1)));
else
    sm_full = zeros(numel(find(any([strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1))),ica_res.ica_params.num_ics);
    map_inds(find(any([strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1))) = 1:numel(find(any([strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1)));
end
sm_full(map_inds(chan_inds),:) = ica_res.sm;


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
else
    badsections = false(1,D.ntrials*D.nsamples);
end

if ~isempty(D.badchannels)
    excluded_data = D(find(D.badchannels),:,setdiff(1:D.ntrials, D.badtrials));
    excluded_data = reshape(excluded_data,size(excluded_data,1),[]); 
    excluded_data(:,badsections)=[];
    sm_full(D.badchannels,:) = excluded_data*pinv(ica_res.tc);
end

ica_res.sm = sm_full;

%%%%%%%%%%%%%%%%%%%%%%%%%% ESTIMATE MISSING EPOCHS %%%%%%%%%%%%%%%%%%%%%%%% 

excluded_timepoints = false(1,D.ntrials*D.nsamples);
for i=1:D.ntrials
    if ismember(i,D.badtrials)
        excluded_timepoints(1+(i-1)*D.nsamples:i*D.nsamples) = 1;
    end
end
excluded_timepoints(badsections)=true;

if strcmp(S.modality,'EEG')  % added by DM
    excluded_data = D(find(any([strcmp(D.chantype,'EEG')],1)),:,:);
else
    excluded_data = D(find(any([strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1)),:,:);
end
excluded_data = reshape(excluded_data,size(excluded_data,1),[]); 

excluded_data = excluded_data(:,excluded_timepoints);

tc_full = zeros(ica_res.ica_params.num_ics,D.ntrials*D.nsamples);
tc_full(:,~excluded_timepoints) = ica_res.tc;

tc_full(:,excluded_timepoints) = (excluded_data'*pinv(sm_full'))';
ica_res.tc = tc_full;

end

%% remove_bad_components - function to subtract the bad components from the
%% MEG data via spm_eeg_montage.m

function res = remove_bad_components(S) % AB - amended montage to prevent channels shuffling. Should mean that we don't need to save the raw sensor labels as should be unchanged.

%%%%%%%%%%%%%%%%%%%%%%%% LOAD AND PREPARE MEG DATA %%%%%%%%%%%%%%%%%%%%%%%%

D = spm_eeg_load(S.fname);

if strcmp(S.modality,'EEG')   % changed by DM
    chan_inds = find(any(strcmp(D.chantype,'EEG'),1));
else
    chan_inds = find(any([strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1));	
end

badchannels    = D.badchannels;
bad_components = unique(S.ica_res.bad_components);
meg_dat        = D(chan_inds,:,:);
meg_dat        = reshape(meg_dat,size(meg_dat,1),[]);

%%%%%%%%%%%%%%%%% SAVE A COPY OF THE TRA MATRIX IF NEEDED %%%%%%%%%%%%%%%%%

if strcmp(S.modality,'MEG')   % added by DM
    D = save_raw_tra_to_D(D);
    use_montage = 1;
elseif strcmp(S.modality,'EEG')
    use_montage = 0;
end

%%%%%%%%%%%%%%%%%%% REMOVE BAD COMPONENTS USING MONTAGE %%%%%%%%%%%%%%%%%%%

sm = S.ica_res.sm;
tc = S.ica_res.tc;

if ~isempty(bad_components)
    if use_montage
        dat_inv = pinv_plus(meg_dat', S.ica_res.ica_params.num_ics);
        tra = (eye(numel(chan_inds)) - dat_inv*(tc(bad_components,:)'*sm(:,bad_components)'))';
        
        montage             =  [];
        montage.tra         =  tra;
        montage.labelnew    =  D.chanlabels(chan_inds);
        montage.labelorg    =  D.chanlabels(chan_inds);
        
        [~,locs] = ismember(montage.labelnew,D.sensors('MEG').label);
        
        
        montage.chanunitnew =  D.sensors('MEG').chanunit(locs);
        montage.chanunitorg =  D.sensors('MEG').chanunit(locs);
        montage.chantypenew =  lower(D.sensors('MEG').chantype(locs));
        montage.chantypeorg =  lower(D.sensors('MEG').chantype(locs));
        
        
        if ~isempty(badchannels) % CHECK THIS!!
            [~,bad_inds] = intersect(montage.labelorg,D.chanlabels(D.badchannels),'stable');
            montage.tra(bad_inds,:) = 0;
            for bi = bad_inds'
                montage.tra(bi,bi) = 1;
            end
        end
        
        S_montage                =  [];
        S_montage.D              =  fullfile(D.path,D.fname);
        S_montage.montage        =  montage;
        S_montage.keepothers     =  true;
        S_montage.updatehistory  =  1;
        
        Dmontaged = spm_eeg_montage(S_montage);
        
        % rename montaged file
        
        S_copy         = [];
        S_copy.D       = Dmontaged;
        S_copy.outfile = fullfile(D.path, ['A' D.fname]);
        Dclean = spm_eeg_copy(S_copy);

        Dmontaged.delete;
        
    else
        [dir,nam,~]=fileparts(fullfile(D.path,D.fname));
        fname_out=[dir '/A' nam '.dat'];
        meg_dat_clean=meg_dat-(sm(:,bad_components)*tc(bad_components,:));
        
        Dclean=clone(D,fname_out,size(D));
        Dclean(chan_inds,:)=meg_dat_clean;   % changed by DM
        Dclean.save;
    end
else
    Dclean=D;
    msg = sprintf('\n%s\n%s\n%','No bad components have been selected. No de-noising has been applied to ', fullfile(D.path, D.fname));
    fprintf(msg);
end

res = fullfile(Dclean.path, Dclean.fname);

end


