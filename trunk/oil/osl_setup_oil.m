function oil = osl_setup_oil(S)

% oil = osl_setup_oil_for_ica(S)
%
% Creates an OAT (OHBA's easy Analysis Tool) struct with all the appropriate settings which can then be passed to
% oil=osl_run_ica(oil); to use the IC pipeline.
%
% INPUTS:
% 
% General Parameters - stored in oil.
% S.oil_name - the string name of the OAT (e.g. S.oil_name = '/home/myusername/4-8Hz_espm8_ctf_2back';).
% S.paradigm - the typoe of data being used. (S.paradigm = 'task'; OR S.paradigm = 'rest'; 
% S.to_do    - sets the oil.to_do (e.g. S.to_do = [1 1 1 1 0 0]).
%
% Source Reconstruction Parameters - copy these across from an oat. 
%
% Enveloping Parameters - stored in oil.enveloping
% S.hilbert_downsampling - the downsampling window length is seconds (default 1s).
% S.spatial_downsampling - the spatial resampling after enveloping (default 8mm).
% S.spatial_smoothing    - FWHM of Gaussian Kernal for spatial smoothing (default 4mm).
% S.timewindow           - the subset time window for the enveloping (to minimise edge effects).
%
% Concatenation Parameters - stored in oil.concat_subs
% S.ica_subs2use - the subset of subjects passed into the beamformer (default is to use all subjects).
% S.normalise_subjects - global variance normalisation of each subject prior to concatenation (default is 0).
%
% ICA Parameters - stored in oil.ica
% S.num_ics - the number of components sought
% S.last_eig - the dimension of the PCA prior to ICA (default equal to num_ics).
% S.temp_or_spat - select either 'temporal' or 'spatial' ICA. (default is temporal).
% S.icasso_its - number of ICASSO iterations (default is 0).
%
% First Level ICA Statistics - stored in oil.ica_first_level. The design
% matrix and contrasts are only relevant here for task data.
% S.design_matrix_summary - this is a parsimonious description of the design matrix.
% S.contrast - E.g. Two different contrasts, for a design matrix with 3 regressors: S.contrast{1}=[1 0 0]'; S.contrast{2}=[0 1 -1]'; 
%
% Group Level ICA Statistics - stored in oil.ica_group_level
% S.group_design_matrix - usually a unity vector whose length is equal to the number of subjects (to give group mean).
% S.group_contrast - usually [1].
%
% HL 15.08.12

%%%%%%%%%%%%%%%%%%%
% General Settings
try oil.fname=S.oil_name; catch, oil.fname=[S.spm_files{1} '.oil']; end
try oil.to_do=S.to_do; catch, oil.to_do=[1 1 1 1 0 0]; end
try oil.paradigm=S.paradigm; catch, error('User must specify the paradigm type: S.paradigm, It can be either rest or task.'); end

%%%%%%%%%%%%%%%%%%%
% source recon settings
oil.source_recon = S.source_recon;

%%%%%%%%%%%%%%%%%%%
% Enveloping Settings
oil.enveloping=[];
try oil.enveloping.window_length=S.hilbert_downsampling; catch, oil.enveloping.window_length=1;end;    
try oil.enveloping.gridstep=S.spatial_downsampling; catch, oil.enveloping.gridstep=8;end;
try oil.enveloping.ss=S.spatial_smoothing; catch, oil.enveloping.ss=4;end;
try oil.enveloping.timewindow=S.envelope_timewindow; catch, oil.enveloping.timewindow='all';end;

%%%%%%%%%%%%%%%%%%%
% Concatenation Settings
oil.concat_subs=[];
try oil.concat_subs.sessions_to_do=S.ica_subs2use; catch,oil.concat_subs.sessions_to_do = 1:length(oil.source_recon.D_continuous);end;
try oil.concat_subs.normalise_subjects=S.normalise_subjects; catch,oil.concat_subs.normalise_subjects=0;end;

%%%%%%%%%%%%%%%%%%%
% ICA Settings
oil.ica=[];
try oil.ica.temp_or_spat=S.temp_or_spat; catch, oil.ica.temp_or_spat='temporal';end;
try oil.ica.use_gm_mask=S.use_gm_mask; catch, oil.ica.use_gm_mask=0;end;
try oil.ica.num_ics=S.num_ics; catch, oil.ica.num_ics=25;end;
try oil.ica.icasso_its=S.icasso_its; catch, oil.ica.icasso_its=0;end;
try oil.ica.last_eig=S.last_eig; catch, oil.ica.last_eig=oil.ica.num_ics;end;
try oil.ica.other_oils=S.other_oils; catch, end;
try oil.ica.nuisance_oil_names=S.nuisance_oil_names; catch, end;
try oil.ica.normalise_vox=S.normalise_vox; catch, oil.ica.normalise_vox=0; end;

%%%%%%%%%%%%%%%%%%%
% First Level Stats Settings
oil.ica_first_level=[];
% required settings:
if strcmp(oil.paradigm,'task')
    try oil.ica_first_level.design_matrix_summary=S.design_matrix_summary; catch, error('S.design_matrix_summary not specified'); end;
    try oil.ica_first_level.contrast=S.contrast; catch, error('S.contrast not specified'); end;
end
oil.ica_first_level.name='ica_first_level';
try, oil.ica_first_level.time_range=S.ica_first_level_time_range; catch, oil.ica_first_level.time_range=oil.source_recon.time_range; end;
try, oil.ica_first_level.use_robust_glm=S.use_robust_glm; catch, oil.ica_first_level.use_robust_glm=0; end; % do robust GLM (uses bisquare via Matlab robustfit fn)
try, oil.ica_first_level.cope_type=S.cope_type; catch, oil.ica_first_level.cope_type='cope'; end % cope type to input to group GLM (from first level), set to 'coape', 'cope', or 'acope'

%%%%%%%%%%%%%%%%%%%
% Group Level Stats Settings

oil.ica_group_level=[];
try, oil.ica_group_level.group_design_matrix=S.group_design_matrix; catch, oil.ica_group_level.group_design_matrix=ones(1,length(oil.source_recon.D)); end;
try, oil.ica_group_level.group_contrast=S.group_contrast; catch, oil.ica_group_level.group_contrast{1}=[1]; end;
try, oil.ica_group_level.comps2use=S.group_comps2use; catch, oil.ica_group_level.comps2use = 1:oil.ica.num_ics; end;

