function [fname_out,fig_handles,fig_names,fig_titles,S] = osl_africa(S)
% AfRICA - ArteFact Rejection using Independent Component Analysis
%
% Written by Henry Luckhoo and Adam Baker


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
if S.todo.ica
    S.ica_res = perform_sensorspace_ica(S);
    if isfield(S,'ica_file')
        if ~isdir(fileparts(S.ica_file))
            mkdir(fileparts(S.ica_file));
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

D = spm_eeg_load(S.fname);

if strcmp(S.modality,'EEG')
    chantype = 'EEG';
else
    chantype = 'MEGANY';
end

% Good channels:
chan_inds = indchantype(D,chantype,'GOOD');

% Good samples:
if D.ntrials == 1
    sample_inds = find(~all(badsamples(D,':',':',1)));
else
    sample_inds = 1:D.nsamples;
end

% Good trials:
trial_inds = indtrial(D,D.condlist,'GOOD');

% Select data:
icadata = D(chan_inds,sample_inds,trial_inds);

% Remove trial structure:
icadata = reshape(icadata,size(icadata,1),[]);


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
ica_res.D                       = S.fname;
ica_res.tc                      = tc;
ica_res.sm                      = sm .* repmat(norm_vec,1,num_ics);

%%%%%%%%%%%%%%%%%%% ESTIMATE MISSING CHANNELS AND EPOCHS %%%%%%%%%%%%%%%%%% 

sm_full = zeros(numel(indchantype(D,chantype)),num_ics);
sm_full(chan_inds,:) = ica_res.sm;

bad_timepoints = all(badsamples(D,':',':',':'));
bad_timepoints = reshape(bad_timepoints,1,D.nsamples*D.ntrials);

if ~isempty(indchantype(D,chantype,'BAD'))
    subdata = D(indchantype(D,chantype,'BAD'),:,:);
    subdata = reshape(subdata,size(subdata,1),[]);
    sm_full(indchantype(D,chantype,'BAD'),:) = subdata(:,~bad_timepoints)*pinv(ica_res.tc);
end

subdata = D(indchantype(D,chantype),:,:);
subdata = reshape(subdata,size(subdata,1),[]);
subdata = subdata(:,bad_timepoints);
tc_full = zeros(ica_res.ica_params.num_ics,D.ntrials*D.nsamples);
tc_full(:,~bad_timepoints) = ica_res.tc;
tc_full(:,bad_timepoints)  = (subdata'*pinv(sm_full'))';

ica_res.sm = sm_full;
ica_res.tc = tc_full;

end


function res = remove_bad_components(S)

%%%%%%%%%%%%%%%%%%%%%%%% LOAD AND PREPARE MEG DATA %%%%%%%%%%%%%%%%%%%%%%%%

D = spm_eeg_load(S.fname);

if strcmp(S.modality,'EEG')
    chantype = 'EEG';
    use_montage = 0;
else
    chantype = 'MEGANY';
    use_montage = 1;
end

% Good channels:
chan_inds = indchantype(D,chantype);

badchannels    = D.badchannels;
bad_components = unique(S.ica_res.bad_components);
megdata        = D(chan_inds,:,:);
megdata        = reshape(megdata,size(megdata,1),[]);

%%%%%%%%%%%%%%%%%%% REMOVE BAD COMPONENTS USING MONTAGE %%%%%%%%%%%%%%%%%%%

sm = S.ica_res.sm;
tc = S.ica_res.tc;

if use_montage
    dat_inv = pinv_plus(megdata', S.ica_res.ica_params.num_ics);
    tra = (eye(numel(chan_inds)) - dat_inv*(tc(bad_components,:)'*sm(:,bad_components)'))';
    
    montage             =  [];
    montage.tra         =  tra;
    montage.labelnew    =  D.chanlabels(chan_inds);
    montage.labelorg    =  D.chanlabels(chan_inds);
    
    [~,indx] = ismember(montage.labelnew,D.sensors('MEG').label);
    
    
    montage.chanunitnew =  D.sensors('MEG').chanunit(indx);
    montage.chanunitorg =  D.sensors('MEG').chanunit(indx);
    montage.chantypenew =  lower(D.sensors('MEG').chantype(indx));
    montage.chantypeorg =  lower(D.sensors('MEG').chantype(indx));
    
    
    if ~isempty(badchannels)
        [~,indx] = ismember(D.sensors('MEG').label,montage.labelnew);
        montage.tra(indx,:) = 0;
        for bi = indx'
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
    meg_dat_clean=megdata-(sm(:,bad_components)*tc(bad_components,:));
    
    Dclean=clone(D,fname_out,size(D));
    Dclean(chan_inds,:) = meg_dat_clean;   % changed by DM
    Dclean.save;
end


res = fullfile(Dclean.path, Dclean.fname);

end


