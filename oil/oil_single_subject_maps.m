function oil = oil_single_subject_maps(oil)
% oil_single_subject_maps.m
%
% Syntax: oil = oil_single_subject_maps(oil);
%
% Function to estimate subject specific COPE and T-statistic maps from 
% group ICA decompositions. OIL stages 1-4 must have been run.
%
% Maps are saved in [oil.source_recon.dirname '/' oil.enveloping.name '/' oil.concat_subs.name '/' oil.ica.name '/' oil.ica_first_level.name].
%
% Maps names are stored in oil.ica_first_level.results.cope_files & oil.ica_first_level.results.tstat_files.
%
% Henry Luckhoo (henry.luckhoo@trinity.ox.ac.uk)
%
% Version 1.6
% 281113

% 28-11-2013 Updated to v 1.6 by Giles Colclough
% (giles.colclough@magd.ox.ac.uk)
% incorporating changes suggested by Duncan Astle to prevent errors when
% subject maps are stored as .mat files, not niftis. 

%% Setup

% Setup stats directory to save the results into
save_dir = [oil.source_recon.dirname '/' oil.enveloping.name '/' oil.concat_subs.name '/' oil.ica.name '/' oil.ica_first_level.name];
if ~isdir(save_dir), mkdir(save_dir); end;

Nics     = oil.ica.num_ics;
Nvoxels  = size(nii.quickread([oil.source_recon.dirname '/' oil.enveloping.name '/' oil.enveloping.results.source_space_envelopes_results_fnames{oil.concat_subs.sessions_to_do(1)}], oil.enveloping.gridstep),1);
subj_ind = oil.concat_subs.results.subj_ind;
Nsubs    = length(subj_ind)-1;


%% Spatial Basis Set Regression
% set-up design matrix (identity).
con = mat2cell(eye(Nics), Nics, ones(1,Nics)); 

% get spatial basis and time courses
isSpatialBasisSet = isfield(oil.ica_first_level, 'spatial_basis_set');
isSpatialICA      = strcmp(oil.ica.temp_or_spat,'spatial');

if isSpatialBasisSet || isSpatialICA, % need to calculate pseudo tcs
    if isSpatialBasisSet, 
        
        % load in spatial basis set
        [~,~,scales] = nii.load(oil.ica_first_level.spatial_basis_set);
        
        if scales(1) ~= oil.enveloping.gridstep;
            oil.ica_first_level.spatial_basis_set = nii.resample(...
                oil.ica_first_level.spatial_basis_set, ...
                [oil.ica_first_level.spatial_basis_set, '_ds', ...
                   num2str(oil.enveloping.gridstep) 'mm'], ...
                oil.enveloping.gridstep);
        end
    
        spat_bas = nii.quickread(oil.ica_first_level.spatial_basis_set, ...
                                 oil.enveloping.gridstep);
        
    else %isSpatialICA still
        spat_bas = transpose(oil.ica.results.sICs);
    end
    
    % load concatenated data (nonorm option appears to be removed)
    % Henry's dual-reg approach suggests using weights-normalised data here

        ica_concat_fname = fullfile(oil.source_recon.dirname, ...
                                    oil.enveloping.name, ...
                                    oil.concat_subs.name, ...
                                    oil.concat_subs.results.concat_file);

    try % assume nifti
        concat_data = nii.quickread(ica_concat_fname, ...
                                    oil.enveloping.gridstep);
    catch ME
        % check to see if saves as .mat
        if strcmp(ME.identifier, 'MATLAB:FileIO:InvalidFid'),
            tmp         = load(ica_concat_fname);
            concat_data = tmp.ica_concat;
            clear tmp;
        else
            rethrow(ME);
        end%if
    end%try   
    
    % estimate timecourses by multiple regression
    Nsamples        = size(concat_data,2);
    Nics            = size(spat_bas,2);
    oil.ica.num_ics = Nics;
    
    x       = demean(spat_bas,1);
    R       = chol(x' * x);
    Rinv    = inv(R);
    pinvxtx = Rinv * Rinv';
    pinvx   = pinvxtx * x';
    
    tc_est  = zeros(Nics,Nsamples); % declare memory
    
    for t = 1:Nsamples;
        y = demean(concat_data(:,t));
        [tc_est(:,t)] = glm_fast_for_meg(y, x, pinvxtx, pinvx, con, 0);
    end
    
    oil.ica_first_level.results.pseudo_tcs = normalise(tc_est,2);
    
    % set mixing matrix to new basis set
    oil.ica_first_level.results.mixing_matrix = spat_bas;
    % select time courses to use
    tcs2use = oil.ica_first_level.results.pseudo_tcs ;
    
else% temporal ica
%     spat_bas = oil.ica_first_level.results.mixing_matrix; % - already set
    tcs2use = oil.ica.results.tICs;
end%if spatial ica

%% Second Stage

comp_mean = zeros(Nics, Nsubs);
comp_var  = zeros(Nics, Nsubs);
resid = zeros(Nvoxels,1);
for subnum=1:Nsubs;
    
    % Subject Specific Setup - use nonorm data
    sing_sub_dat=nii.quickread([oil.source_recon.dirname '/' oil.enveloping.name '/' oil.enveloping.results.source_space_envelopes_NoWeightsNorm_results_fnames{oil.concat_subs.sessions_to_do(subnum)}],oil.enveloping.gridstep);
    sing_sub_tc=normalise(tcs2use(:,subj_ind(subnum):subj_ind(subnum+1)-1)',1);  % extract subject specific part of the ICA time courses and normalise to unit std and mean.
        
    % Multiple Regression
    x=sing_sub_tc;
    pinvxtx=pinv(x'*x);
    pinvx=pinv(x);
    copeout=zeros(Nics,Nvoxels);varcopeout=zeros(Nics,Nvoxels);
    for v=1:Nvoxels
        y=demean(sing_sub_dat(v,:)');
        [copeout(:,v), varcopeout(:,v)]=glm_fast_for_meg(y,x,pinvxtx,pinvx,con,0);
        resid(v) = sqrt(sum((y - x*copeout(:,v)).^2));
    end
    
    comp_var(:,subnum)  = var(copeout,[],2);
    comp_mean(:,subnum)=pinv(normalise(copeout',1))*mean(sing_sub_dat,2);
   
    
    % Save COPES and tstats to NII files
    cope_files{subnum} = [save_dir '/ica_copes_' oil.enveloping.results.source_space_envelopes_results_fnames{oil.concat_subs.sessions_to_do(subnum)}];
    nii.quicksave(copeout',cope_files{subnum},oil.enveloping.gridstep);
    tstat_files{subnum} = [save_dir '/ica_tstats_' oil.enveloping.results.source_space_envelopes_results_fnames{oil.concat_subs.sessions_to_do(subnum)}];
    nii.quicksave(copeout'./sqrt(varcopeout'),tstat_files{subnum},oil.enveloping.gridstep);
end

%% 
oil.ica_first_level.results.component_means     = comp_mean;
oil.ica_first_level.results.component_variances = comp_var;
oil.ica_first_level.results.cope_files          = cope_files;
oil.ica_first_level.results.tstat_files         = tstat_files;
end%osl_single_subject_maps