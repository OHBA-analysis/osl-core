function [fname_out,fig_handles,fig_names,fig_titles,S] = osl_africa(S)
% AfRICA - ArteFact Rejection using Independent Component Analysis
% performs ICA denoising of MEG data using either semi-manual or automated 
% identification artefact components.
%
% [fname_out,fig_handles,fig_names,fig_titles,S] = osl_africa(S)
%
% REQUIRED INPUTS:
%
% S.D           - SPM MEG object filename
%
% S.ica_file    - .mat file in which to save ica results
%
%
% OPTIONAL INPUTS 
%
% S.todo        - structure with fields:
%                   .ica    - [0/1] to run ica decomposition
%                   .ident  - [0/1] to run artefact identification
%                   .remove - [0/1] to run artefact removal
%
% S.ident       - structure with fields:
%                   .func - function handle to identification function
%                           e.g. @identify_artefactual_components_manual
%                   .{extra fields depending on .func}                   
%
% S.logfile     - write logfile [0/1], default = 1
%
% S.modality    - modality to use, default = 'MEG'
%
% S.do_plots    - produce diagnostic plots, default = 0
%
% S.used_maxfilter - [0/1] if Maxfilter has been used, default = 0
%
% S.precompute_topos   - pre-compute and save IC spatial map topos after ica is computed for use in ident
%
% Written by Henry Luckhoo and Adam Baker



% Check SPM File Specification:
try
    S.D = char(S.D);
    [pathstr,filestr] = fileparts(S.D);
    S.D = fullfile(pathstr,[filestr '.mat']); % force .mat suffix
    D = spm_eeg_load(S.D);
catch
    if isfield(S,'fname')
        warning('S.fname is deprecated, please use S.D instead')
    end
    error('SPM file specification not recognised or incorrect');
end

if ~isfield(S,'ica_file') || isempty(S.ica_file)
    error('No file specified in S.ica_file for saving AFRICA results')
end

[pathstr,name] = fileparts(S.ica_file);
if isempty(pathstr)
    pathstr = D.path;
end
S.ica_file = fullfile(pathstr,[name '.mat']);

if ~isdir(pathstr) % Create directory if it doesn't exist
    mkdir(pathstr);
end

if isfield(S,'logfile') && S.logfile == 1
    logfile = fullfile(pathstr,[name '_log.txt']);
  if exist(logfile,'file')
    unix(['rm ' logfile]);
  end
  	diary(logfile);
end

if not(isfield(S,'modality'))
    S.modality = 'MEG';
end

if not(isfield(S,'do_plots'))
    S.do_plots = 0;
end

if ~isfield(S,'precompute_topos');
    S.precompute_topos = 0; 
end

if ~isfield(S,'todo');
    S.todo = struct; 
end

if ~isfield(S.todo,'ica');
    S.todo.ica = 1; 
end

if ~isfield(S.todo,'ident');
    S.todo.ident = 1; 
end

if ~isfield(S.todo,'remove');
    S.todo.remove = 1; 
end

if ~isfield(S,'ident')
    S.ident = [];
end

if ~isfield(S.ident,'func')
    S.ident.func = @identify_artefactual_components_manual;
end

if ~isfield(S,'used_maxfilter'); 
    S.used_maxfilter = 0;
end

fig_handles = []; 
fig_names   = [];
fig_titles  = [];


%%%%%%%%%%%%%%%%%%%%%%%%% PERFORM SENSOR SPACE ICA %%%%%%%%%%%%%%%%%%%%%%%%
if S.todo.ica
    S.ica_res = perform_sensorspace_ica(S);
    if isfield(S,'ica_file')
                
        % Precompute topographies
        if S.precompute_topos
            try
                topos = S.ica_res.topos;
            catch
                topos = [];
                sm = S.ica_res.sm;
                modalities = unique(D.chantype(find(strncmpi(S.modality,D.chantype,3)))); %#ok
                for m = 1:numel(modalities)
                    disp(['Precomputing sensor topographies for modality ' modalities{m}]);
                    topos = [topos component_topoplot(D,sm,modalities(m))];
                end
            end
            S.ica_res.topos   = topos;
        end
        
        save(S.ica_file,'S');
        msg = sprintf('\n%s%s\n%','Saving ICA results to ', S.ica_file);
        fprintf(msg);

    else
        msg = sprintf('\n%s\n%','Results not being saved.');
        fprintf(msg);
    end

elseif any(structfun(@istrue,S.todo))
    % Load from file
    if isfield(S,'ica_file') && exist(S.ica_file,'file')==2
        msg = sprintf('\n%s%s\n%','Loading previous ICA results from ', S.ica_file);
        fprintf(msg);
        icafile = load(S.ica_file);
        if isfield(S,'ident');    icafile.S.ident    = S.ident;    end
        if isfield(S,'do_plots'); icafile.S.do_plots = S.do_plots; end
        if isfield(S,'todo');    icafile.S.todo    = S.todo;    end
        S = icafile.S; clear icafile
    end 
end


%%%%%%%%%%%%%%%%%%%%%%%%%% CLASSIFY BAD COMPONENTS %%%%%%%%%%%%%%%%%%%%%%%%
if S.todo.ident
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
        error('ICA decomposition needs to be run. Set S.todo.ica = 1 and try again');
    end
end


%%%%%%%%%%%%%%%%%% REMOVE BAD COMPONENTS FROM THE DATA %%%%%%%%%%%%%%%%%%%%
if S.todo.remove
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

end % MAIN FUNCTION







function ica_res = perform_sensorspace_ica(S)

D = spm_eeg_load(S.D);

if strcmp(S.modality,'EEG')
    chantype = 'EEG';
else
    chantype = {'MEG','MEGANY'}; % need both flags to capture MEG, MEGMAG, MEGGRAD and MEGPLANAR channel types.
end

% Good channels:
chan_inds = indchantype(D,chantype,'GOOD');

% Good timepoints/trials
good_samples = ~all(badsamples(D,':',':',':'));
good_samples = reshape(good_samples,1,D.nsamples*D.ntrials);

% Select data:
icadata = D(chan_inds,:,:);

% Remove trial structure:
icadata = reshape(icadata,size(icadata,1),[]);

% Select good timepoints
icadata = icadata(:,good_samples);

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
    num_ics       = num_ics_default;
    last_eig      = num_ics_default;
    nonlinearity  = 'tanh';
    ica_approach  = 'symm';
    stabilization = 'on';
    max_iter      = 1000;
end

num_ics  = min(num_ics, size(icadata,1)); % added by DM
last_eig = min(last_eig,size(icadata,1)); % added by DM

%%%%%%%%%%%%%%%%%%%%  MINIMUM EIGENVALUE NORMALISATION %%%%%%%%%%%%%%%%%%%%
if strcmp(S.modality,'EEG')  % added by DM
    norm_vec = max(abs(icadata(:)))/1000*ones(size(icadata,1),1);
else
    norm_vec = ones(numel(chan_inds),1);
    if any(strcmp(D.chantype,'MEGMAG')) && any(strcmp(D.chantype,'MEGPLANAR'))
        mag_min_eig = svd(cov(icadata(strcmp(D.chantype(chan_inds),'MEGMAG'),:)')); 
        mag_min_eig = mean(mag_min_eig(mag_cutoff-2:mag_cutoff));
        
        plan_min_eig = svd(cov(icadata(strcmp(D.chantype(chan_inds),'MEGPLANAR'),:)')); 
        plan_min_eig = mean(plan_min_eig(plan_cutoff-2:plan_cutoff));
        
        norm_vec(strcmp(D.chantype(chan_inds),'MEGMAG'))    = mag_min_eig; 
        norm_vec(strcmp(D.chantype(chan_inds),'MEGPLANAR')) = plan_min_eig;
    else
        norm_vec = norm_vec*min(svd(cov(icadata(:,:)')));
    end
    norm_vec = sqrt(norm_vec);
    
end

eigs_preNorm  = svd(cov(icadata'));

% Apply normalisation
icadata = icadata ./ repmat(norm_vec,1,size(icadata,2));

eigs_postNorm = svd(cov(icadata'));

if S.do_plots
    figure; 
    semilogy(eigs_preNorm);
    ho;
    semilogy(eigs_postNorm,'r--');
    title('Raw and normalised eigen spectra'); legend('Raw', 'Normalised');
end

%%%%%%%%%%%%%%%%%% AUTOMATIC DIMENSIONALITY ESTIMATION %%%%%%%%%%%%%%%%%%%%
if 0 == num_ics
    num_ics = spm_pca_order(icadata');
end
if 0 == last_eig
    last_eig = num_ics;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% ICA DECOMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[tc,sm,~] = fastica(icadata,                           ...
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
ica_res.D                       = S.D;
ica_res.tc                      = tc;
ica_res.sm                      = sm .* repmat(norm_vec,1,num_ics);

%%%%%%%%%%%%%%%%%%% ESTIMATE MISSING CHANNELS AND EPOCHS %%%%%%%%%%%%%%%%%% 

sm_full              = zeros(D.nchannels, size(sm,2));
sm_full(chan_inds,:) = ica_res.sm;

bad_timepoints = all(badsamples(D,':',':',':'));
bad_timepoints = reshape(bad_timepoints,1,D.nsamples*D.ntrials);

if ~isempty(indchantype(D,chantype,'BAD'))
    subdata = D(indchantype(D,chantype,'BAD'),:,:);
    subdata = reshape(subdata,size(subdata,1),[]);
    sm_full(indchantype(D,chantype,'BAD'),:) = subdata(:,~bad_timepoints)*pinv(ica_res.tc);
end

subdata = D(chan_inds,:,:);
subdata = reshape(subdata,size(subdata,1),[]);
subdata = subdata(:,bad_timepoints); % why is subdat aempty?
tc_full = zeros(ica_res.ica_params.num_ics,D.ntrials*D.nsamples);
tc_full(:,~bad_timepoints) = ica_res.tc;
tc_full(:,bad_timepoints)  = (subdata'*pinv(sm_full(chan_inds,:)'))';

ica_res.sm = sm_full;
ica_res.tc = tc_full;

end


function res = remove_bad_components(S)

%%%%%%%%%%%%%%%%%%%%%%%% LOAD AND PREPARE MEG DATA %%%%%%%%%%%%%%%%%%%%%%%%

D = spm_eeg_load(S.D);

if strcmp(S.modality,'EEG')
    chantype = 'EEG';
    modality = 'EEG';
    use_montage = 1;
else
    chantype = {'MEG','MEGANY'};
    modality = 'MEG';
    use_montage = 1;
end

% Good channels:
chan_inds = indchantype(D,chantype, 'GOOD');

badchannels    = D.badchannels;
bad_components = unique(S.ica_res.bad_components);
megdata        = D(chan_inds,:,:);
megdata        = reshape(megdata,size(megdata,1),[]);

%%%%%%%%%%%%%%%%%%% REMOVE BAD COMPONENTS USING MONTAGE %%%%%%%%%%%%%%%%%%%

sm = S.ica_res.sm;
tc = S.ica_res.tc;

tra = eye(D.nchannels);

if use_montage
    dat_inv = pinv_plus(megdata', S.ica_res.ica_params.num_ics);
    tra(chan_inds,chan_inds) = (eye(numel(chan_inds)) - dat_inv*(tc(bad_components,:)'*sm(chan_inds,bad_components)'))';
    
    % Set up montage with full channel list:
    montage             =  [];
    montage.tra         =  tra;
    montage.labelnew    =  D.chanlabels;
    montage.labelorg    =  D.chanlabels;
    
    % Exclude channels that are not MEEG:
    xchans = setdiff(1:D.nchannels,indchantype(D,chantype));
    montage.tra(xchans,:)    = [];
    montage.tra(:,xchans)    = [];
    montage.labelorg(xchans) = [];
    montage.labelnew(xchans) = [];
    
    [~,indx] = ismember(montage.labelnew,D.sensors(modality).label);
            
    montage.chanunitnew =  D.sensors(modality).chanunit(indx);
    montage.chanunitorg =  D.sensors(modality).chanunit(indx);
    montage.chantypenew =  D.sensors(modality).chantype(indx);
    montage.chantypeorg =  D.sensors(modality).chantype(indx);

        
    S_montage                =  [];
    S_montage.D              =  fullfile(D.path,D.fname);
    S_montage.montage        =  montage;
    S_montage.keepothers     =  true;
    S_montage.updatehistory  =  1;
    
    Dmontaged = osl_montage(S_montage);
    
    % rename montaged file
    S_copy         = [];
    S_copy.D       = Dmontaged;
    S_copy.outfile = fullfile(D.path, ['A' D.fname]);
    Dclean = spm_eeg_copy(S_copy);
    
    Dmontaged.delete;
    
else
    [dir,nam,~] = fileparts(fullfile(D.path,D.fname));
    fname_out = [dir '/A' nam '.dat'];
    meg_dat_clean = megdata-(sm(chan_inds,bad_components)*tc(bad_components,:));
    
    Dclean=clone(D,fname_out,size(D));
    Dclean(chan_inds,:) = meg_dat_clean;   % changed by DM
    Dclean.save;
end


res = fullfile(Dclean.path, Dclean.fname);

end


