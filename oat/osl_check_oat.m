function oat = osl_check_oat(oatin)
% oat = osl_check_oat(oat)
%
% Checks an OAT (OHBA's easy Analysis Tool) struct has all the appropriate
% settings, which can then be passed to osl_run_oat to do an OAT
% analysis. Throws an error if any required inputs are missing, fills other
% settings with default values.
%
% Required inputs:
%
% For oat.source_recon:
%
% oat.source_recon.time_range; Time range (from to) within trial, in secs,
% need to all be the same duration and one for each condition.
% 
% oat.source_recon.conditions; list of conditions to include in the
% analysis, e.g. oat.source_recon.conditions={'Motorbike','Neutral face'};
%
% oat.source_recon.D_continuous; list of continuous time SPM MEEG object
% file paths to run the analysis on (list order should correspond to
% oat.source_recon.mri and oat.source_recon.D_epoched fields if provided)
% AND/OR
% oat.source_recon.D_epoched; list of epoched SPM MEEG object file paths to
% run the analysis on (list order should correspond to oat.source_recon.mri
% and oat.source_recon.D_epoched fields if provided)
%
% e.g. oat.source_recon.D_epoched{1}='subject1';
% oat.source_recon.D_epoched{2}='subject2' 
%
% Optional inputs:
%
% See inside this function (e.g. use "type osl_check_oat") to see the other
% optional settings, or just look at the fields in the output oat!
%
% MWW 2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optional oat settings:
try oat.to_do=oatin.to_do; oatin = rmfield(oatin,'to_do'); catch, oat.to_do=[1 1 1 1]; end;
try oat.do_plots=oatin.do_plots; oatin = rmfield(oatin,'do_plots'); catch, oat.do_plots=0; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% source recon settings
oat.source_recon=[];

% required settings:
try oat.source_recon.time_range = oatin.source_recon.time_range; oatin.source_recon = rmfield(oatin.source_recon,'time_range'); oat.source_recon=ft_checkopt(oat.source_recon,'time_range',{'ascendingdoublebivector'}); catch, error('oatin.source_recon.time_range not specified correctly'); end; % Time range (from to) within trial, in secs, need to all be the same duration and one for each condition
try oat.source_recon.conditions = oatin.source_recon.conditions; oatin.source_recon = rmfield(oatin.source_recon,'conditions'); catch, error('oatin.source_recon.conditions not specified, list of conditions (conditions) to include in the analysis'); end; % e.g. conditions={'Motorbike','Neutral face');
try oat.source_recon.D_continuous=oatin.source_recon.D_continuous; oatin.source_recon = rmfield(oatin.source_recon,'D_continuous'); catch, oat.source_recon.D_continuous=[]; end; % list of continuous SPM MEEG object file paths to run the analysis on (list order should correspond to oat.source_recon.mri and oat.source_recon.D_epoched fields if provided)
% AND/OR:
try oat.source_recon.D_epoched=oatin.source_recon.D_epoched; oatin.source_recon = rmfield(oatin.source_recon,'D_epoched'); catch, oat.source_recon.D_epoched=[]; end; % list of epoched SPM MEEG object file paths - this will be used to obtain the bad channel and bad trial info. (list order should correspond to oat.source_recon.mri and oat.source_recon.D_continuous fields if provided)

% check list of SPM MEEG filenames input
if(~isempty(oat.source_recon.D_epoched)),
    Ds=oat.source_recon.D_epoched;
elseif(~isempty(oat.source_recon.D_continuous)),
    Ds=oat.source_recon.D_continuous;
else
    error('Either oat.source_recon.D_continuous, or oat.source_recon.D_epoched need to be specified');
end;
num_sessions=length(Ds);

% auto detect datatype so that defaults can be set intelligently

chant='unknown';
try
    D=spm_eeg_load(Ds{1});
    chant=D.chantype;
catch, end;

if sum(ismember(unique(chant),'MEGGRAD')), 
    disp('Detected CTF data. Using default CTF settings');  
    datatype='ctf';
elseif sum(ismember(unique(chant),'MEGPLANAR')) && ~sum(ismember(unique(chant),'MEGPCACOMP'))
    disp('Detected Elekta Neuromag 306 data. Using default Elekta Neuromag 306 settings.'); 
    datatype='neuromag';
elseif sum(ismember(unique(chant),'EEG')), 
    disp('Detected EEG data. Using default EEG settings');  
    datatype='eeg';
else
    disp('Unknown datatype. Using default Elekta Neuromag 306 settings');  
    datatype='neuromag';    
end

% optional settings:
try oat.source_recon.freq_range = oatin.source_recon.freq_range; oatin.source_recon = rmfield(oatin.source_recon,'freq_range'); oat.source_recon=ft_checkopt(oat.source_recon,'freq_range',{'ascendingdoublebivector'}); catch, oat.source_recon.freq_range=[]; end; % Frequency range (from to) in Hz
try oat.source_recon.dirname=oatin.source_recon.dirname; oatin.source_recon = rmfield(oatin.source_recon,'dirname'); catch, oat.source_recon.dirname=[Ds{1} '.oat']; end; % directory name which will be created and within which all results associated with this source recon will be stored
try oat.source_recon.sessions_to_do=oatin.source_recon.sessions_to_do; oatin.source_recon = rmfield(oatin.source_recon,'sessions_to_do'); catch, oat.source_recon.sessions_to_do=1:num_sessions; end; % list of sessions indexes to run source_recon on.
try oat.source_recon.epochinfo=oatin.source_recon.epochinfo; oatin.source_recon = rmfield(oatin.source_recon,'epochinfo'); catch, oat.source_recon.epochinfo=[]; end; % epochinfo for passing to spm_eeg_epoch (see help spm_eeg_epoch). if empty, will look inside D_epoched for D_epoched.epochinfo
try oat.source_recon.method=oatin.source_recon.method; oatin.source_recon = rmfield(oatin.source_recon,'method'); oat.source_recon = ft_checkopt(oat.source_recon,'method','char',{'none','beamform','beamform_bilateral','mne_eye','mne_diag_datacov'}); catch,  oat.source_recon.method='beamform'; warning(['oat.source_recon.method either not specified or incorrectly specified, Setting it to: ' oat.source_recon.method]); end; % can be 'beamform' or 'none' or 'beamform_bilateral' (for a sensor space analysis)
try oat.source_recon.bandstop_filter_mains=oatin.source_recon.bandstop_filter_mains; oatin.source_recon = rmfield(oatin.source_recon,'bandstop_filter_mains'); catch, oat.source_recon.bandstop_filter_mains = 0; end % defaults to applying a bandstop filter at 50Hz and higher harmonics
try oat.source_recon.artefact_chanlabel = oatin.source_recon.artefact_chanlabel; oatin.source_recon = rmfield(oatin.source_recon,'artefact_chanlabel'); catch, oat.source_recon.artefact_chanlabel = []; end % the label of the channel in the input data which contains a vector picking out periods with blinks / artefact. If empty, assumed no artefact rejection desired.

try oat.source_recon.modalities=oatin.source_recon.modalities; oatin.source_recon = rmfield(oatin.source_recon,'modalities'); 
    oat.source_recon = ft_checkopt(oat.source_recon,'modalities','cell',{{'MEGMAG'},{'MEGPLANAR'},{'MEGMAG';'MEGPLANAR'}});      
catch, 
    switch datatype
        case 'neuromag'
            oat.source_recon.modalities={'MEGPLANAR', 'MEGMAG'}; 
        case 'ctf' 
            oat.source_recon.modalities={'MEG'}; 
        case 'eeg' 
            oat.source_recon.modalities={'EEG'}; 
        otherwise
            error('datatype not set properly. \n'); % should not be here
    end;
    warning(['oat.source_recon.modalities not set, or not set properly. Will set to default:' oat.source_recon.modalities{:}]);
    
end; % modalities to include

if isfield(oatin.source_recon, 'session_names'),                           %GC 2014
    % This defines the root filename for each subject/session
    oat.source_recon.session_names = oatin.source_recon.session_names; 
    oatin.source_recon             = rmfield(oatin.source_recon, ...
                                             'session_names');
else
    % Define a generic default - the file roots go session1, session2, ...
    oat.source_recon.session_names = arrayfun(@(a) sprintf('session%d', a), ...
                                              (1:num_sessions)', ...
                                              'UniformOutput', false);
end%if session_names defined
try oat.source_recon.normalise_method=oatin.source_recon.normalise_method; oatin.source_recon = rmfield(oatin.source_recon,'normalise_method'); catch, oat.source_recon.normalise_method='mean_eig'; end; % normalisation method to normalise the different channel modalities. Either 'min_eig' (recommended) or 'mean_eig' or 'none' (not recommended)
try oat.source_recon.regpc=oatin.source_recon.regpc; oatin.source_recon = rmfield(oatin.source_recon,'regpc'); catch, oat.source_recon.regpc=0; end; % data covariance regularisation
try oat.source_recon.gridstep=oatin.source_recon.gridstep; oatin.source_recon = rmfield(oatin.source_recon,'gridstep'); catch, oat.source_recon.gridstep=7; end; % space between dipoles in beamformer grid in mm
try oat.source_recon.mask_fname=oatin.source_recon.mask_fname; oatin.source_recon = rmfield(oatin.source_recon,'mask_fname'); catch, end;
try oat.source_recon.forward_meg = oatin.source_recon.forward_meg; oatin.source_recon = rmfield(oatin.source_recon,'forward_meg'); catch, oat.source_recon.forward_meg='Single Shell'; end % MEG head forward model set to 'Single Shell' or 'MEG Local Spheres'
  
% pca_dim is a num_sessions vector with the rank for each session to be use for the pca
% dimensionality reduction. If pca_dim=-1 then spm_pca_order is used to estimate pca_dim. 
% If oat.source_recon.do_reduce_pca=1 then this is used prior to any source recon, 
% if oat.source_recon.do_reduce_pca=0 then this is only used to determine the rank of the beamformer data covariance. 
% Need to make sure that this includes a correction for the number of ICs removed by AFRICA:
switch datatype
    case 'neuromag'
        try oat.source_recon.pca_dim=oatin.source_recon.pca_dim; oatin.source_recon = rmfield(oatin.source_recon,'pca_dim'); catch, oat.source_recon.pca_dim=50; end; 
    case 'ctf' 
        try oat.source_recon.pca_dim=oatin.source_recon.pca_dim; oatin.source_recon = rmfield(oatin.source_recon,'pca_dim'); catch, oat.source_recon.pca_dim=260; end; 
    case 'eeg' 
        try oat.source_recon.pca_dim=oatin.source_recon.pca_dim; oatin.source_recon = rmfield(oatin.source_recon,'pca_dim'); catch, oat.source_recon.pca_dim=-1; end; 
end;
if length(oat.source_recon.pca_dim)>1
    if length(oat.source_recon.pca_dim)~= num_sessions,
        error('length(oat.source_recon.pca_dim)~= num_sessions');
    end;
else
    oat.source_recon.pca_dim=ones(num_sessions,1)*oat.source_recon.pca_dim;
end;

try oat.source_recon.force_pca_dim = oatin.source_recon.force_pca_dim;  oatin.source_recon = rmfield(oatin.source_recon,'force_pca_dim'); 
catch,  
    switch datatype
        case 'neuromag'
            oat.source_recon.force_pca_dim=0; 
        otherwise
            oat.source_recon.force_pca_dim=0; 
    end;
end; % Force the beamformer to use the initial pca_dim rather than attempting to correct it.

try oat.source_recon.mni_coords            = oatin.source_recon.mni_coords;            oatin.source_recon = rmfield(oatin.source_recon,'mni_coords');            catch,                                             end % num_coords x 3
try oat.source_recon.type                  = oatin.source_recon.type;                  oatin.source_recon = rmfield(oatin.source_recon,'type');                  oat.source_recon = ft_checkopt(oat.source_recon,'type','char',{'Scalar','Vector'}); catch, oat.source_recon.type='Scalar';             end % 'Scalar' or 'Vector'
try oat.source_recon.bandstop_filter_mains = oatin.source_recon.bandstop_filter_mains; oatin.source_recon = rmfield(oatin.source_recon,'bandstop_filter_mains'); catch, oat.source_recon.bandstop_filter_mains = 0; end % defaults to applying a bandstop filter at 50Hz and higher harmonics

% HMM optional settings
try oat.source_recon.hmm_num_states = oatin.source_recon.hmm_num_states; oatin.source_recon = rmfield(oatin.source_recon,'hmm_num_states'); catch, oat.source_recon.hmm_num_states = 0; end; % 0 indicates not to do any HMM stuff, -1 indicates to adaptively determine the number of states to use based on oat.source_recon.hmm_av_class_occupancy
try oat.source_recon.hmm_num_starts = oatin.source_recon.hmm_num_starts; oatin.source_recon = rmfield(oatin.source_recon,'hmm_num_starts'); catch, oat.source_recon.hmm_num_starts = 1; end; 
try oat.source_recon.hmm_pca_dim    = oatin.source_recon.hmm_pca_dim;    oatin.source_recon = rmfield(oatin.source_recon,'hmm_pca_dim');    catch, oat.source_recon.hmm_pca_dim    = 40; end; 
try oat.source_recon.hmm_block      = oatin.source_recon.hmm_block;      oatin.source_recon = rmfield(oatin.source_recon,'hmm_block');      catch, end; 
try oat.source_recon.hmm_av_class_occupancy=oatin.source_recon.hmm_av_class_occupancy; oatin.source_recon = rmfield(oatin.source_recon,'hmm_av_class_occupancy');  catch, oat.source_recon.hmm_av_class_occupancy=35; end; % average HMM occupancy to aim for (in secs) when using the oat.source_recon.hmm_num_states=-1 option

% report settings
if ~isfield(oatin.source_recon,'report')
    oatin.source_recon.report=struct;
end;
oat.source_recon.report=struct;
try oat.source_recon.report.do_source_variance_maps=oatin.source_recon.report.do_source_variance_maps; oatin.source_recon.report = rmfield(oatin.source_recon.report,'do_source_variance_maps'); 
catch
    oat.source_recon.report.do_source_variance_maps=1;
end;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FIRST LEVEL (WITHIN SUBJECT) GLM settings

if ~isfield(oatin,'first_level')
    oatin.first_level=struct;
end;
oat.first_level=struct;

% no required settings

% optional settings:

% GLM stuff:
%
% Either set oatin.first_level.design_matrix explicitly (usually for continuous data),
% or set oatin.first_level.design_matrix_summary (usually for continuous
% data):
%
% 1) oatin.first_level.design_matrix_summary is a parsimonious description of the design matrix.
% It contains a vector for each regressor, i.e. oatin.first_level.design_matrix_summary{reg},
% where reg is a regressor no., and the vector contains a value for each
% condition/trigger that will be duplicated over all trials of that
% condition/trigger in the resulting design matrix to create the (num_regressors x num_trials)
% design matrix.
% E.g. for an analysis with two conditions then use the following to create a
% design matrix with two regressors that pick out trials belonging to condition 1 and condition 2:
% oatin.first_level.design_matrix_summary={}; oatin.first_level.design_matrix_summary{1}=[1 0];oatin.first_level.design_matrix_summary{2}=[0 1];');
%
% 2) oatin.first_level.design_matrix_summary can contain a list of full file paths to text files containing each
% subject's design matrix. Each text file is num_trials x num_regressors for a subject,
% where the num_trials (and trial order) assumes that the D.badtrials trials have been removed.
% outputs the design matrix, x, which is num_trials x num_regressorsIn which case you also need to set
% oatin.first_level.trial_rejects, which is a text file containing a list of trial indices ( for
% the loaded in Xsummary design matrix ) to indicate any further trials ( above and beyond the D.badtrials trials )
% that you do not want to include in the analysis (e.g. for behavioural reasons), and which will get set to 0 in the design matrix.

try oat.first_level.doGLM=oatin.first_level.doGLM; oatin.first_level = rmfield(oatin.first_level,'doGLM'); catch, oat.first_level.doGLM=1; end; % flag to determine whether glm is run or not, if not run then the glm_input data is output instead.
try oat.first_level.design_matrix=oatin.first_level.design_matrix; oatin.first_level = rmfield(oatin.first_level,'design_matrix'); catch, end; % design matrix should be num_PES x num_tpts
try oat.first_level.design_matrix_summary=oatin.first_level.design_matrix_summary; oatin.first_level = rmfield(oatin.first_level,'design_matrix_summary'); catch, end;

% contrasts:
try oat.first_level.contrast=oatin.first_level.contrast;
    oatin.first_level = rmfield(oatin.first_level,'contrast');
catch
    if(isfield(oat.first_level,'design_matrix')), oat.first_level.contrast{1}=ones(size(oat.first_level.design_matrix,1),1);
    elseif(isfield(oat.first_level,'design_matrix_summary')), oat.first_level.contrast{1}=ones(length(oat.first_level.design_matrix_summary{1}),1);
    else oat.first_level.contrast{1}=1;
    end;
end;
try oat.first_level.contrast_name=oatin.first_level.contrast_name; oatin.first_level = rmfield(oatin.first_level,'contrast_name'); catch, if(isfield(oat.first_level,'contrast')),for c=1:length(oat.first_level.contrast), oat.first_level.contrast_name{c}=['c' num2str(c)];end;end;end; % Each contrast should be [num_regressors x 1], e.g. contrast{1}=[1 -1]';

% other settings:
try oat.first_level.sessions_to_do=oatin.first_level.sessions_to_do; oatin.first_level = rmfield(oatin.first_level,'sessions_to_do'); catch, oat.first_level.sessions_to_do=oat.source_recon.sessions_to_do; end; % list of subject indexes to run first_level on.
try oat.first_level.time_range=oatin.first_level.time_range; oatin.first_level = rmfield(oatin.first_level,'time_range'); catch, oat.first_level.time_range=oat.source_recon.time_range; end;
try oat.first_level.baseline_timespan=oatin.first_level.baseline_timespan; oatin.first_level = rmfield(oatin.first_level,'baseline_timespan'); catch, if(~isempty(oat.first_level.time_range)), oat.first_level.baseline_timespan=[oat.first_level.time_range(1) 0]; else, oat.first_level.baseline_timespan=[];end; end;
try oat.first_level.mask_fname=oatin.first_level.mask_fname; oatin.first_level = rmfield(oatin.first_level,'mask_fname'); catch, end;
try oat.first_level.do_weights_normalisation=oatin.first_level.do_weights_normalisation; oatin.first_level = rmfield(oatin.first_level,'do_weights_normalisation'); catch, oat.first_level.do_weights_normalisation=1; end; % do beamformer weights normalisation

try oat.first_level.use_robust_glm=oatin.first_level.use_robust_glm; oatin.first_level = rmfield(oatin.first_level,'use_robust_glm'); catch, oat.first_level.use_robust_glm=0; end; % do robust GLM (uses bisquare via Matlab robustfit fn)
try oat.first_level.bc=oatin.first_level.bc; oatin.first_level = rmfield(oatin.first_level,'bc'); catch,
    try
        oat.first_level.bc=ones(size(oat.first_level.contrast,2),1);
    catch
        oat.first_level.bc=[];
    end;
end;% (num_contrasts x 1) binary vector indicating whether to do baseline correction after fitting the GLM on the COPE timecourse at each voxel (contrast-wise)
try oat.first_level.sensor_space_combine_planars = oatin.first_level.sensor_space_combine_planars; oatin.first_level = rmfield(oatin.first_level,'sensor_space_combine_planars'); % can be 'dont_combine' or 'combine_cartesian' (e.g. if neuromag).
catch, 
    if strcmp(datatype,'neuromag'),
        oat.first_level.sensor_space_combine_planars = 'combine_cartesian'; 
    else
        oat.first_level.sensor_space_combine_planars = 'dont_combine'; 
    end;
    
end % how to treat the planar gradiometers if running a sensorspace oat

try oat.first_level.bc_trialwise=oatin.first_level.bc_trialwise; oatin.first_level = rmfield(oatin.first_level,'bc_trialwise'); catch, oat.first_level.bc_trialwise=0; end;% binary number indicating whether to do baseline correction for each trial separately at each voxel on the data BEFORE fitting the GLM. If doing TF analysis this is applied to the power time courses at each freq.
try oat.first_level.trial_rejects=oatin.first_level.trial_rejects; oatin.first_level = rmfield(oatin.first_level,'trial_rejects'); catch, end; % used to indicate further rejected trials if setting design matrix (via oat.first_level.design_matrix_summary) using a file
try oat.first_level.name=oatin.first_level.name; oatin.first_level = rmfield(oatin.first_level,'name'); catch, oat.first_level.name='first_level'; end;
try oat.first_level.time_moving_av_win_size = oatin.first_level.time_moving_av_win_size; oatin.first_level = rmfield(oatin.first_level,'time_moving_av_win_size'); catch, oat.first_level.time_moving_av_win_size = []; end; % temporal moving average window size in seconds. 
try oat.first_level.post_glm_time_moving_av_win_size = oatin.first_level.post_glm_time_moving_av_win_size; oatin.first_level = rmfield(oatin.first_level,'post_glm_time_moving_av_win_size'); catch, oat.first_level.post_glm_time_moving_av_win_size = []; end; % temporal moving average window size in seconds. 
try oat.first_level.mni_coords = oatin.first_level.mni_coords;  oatin.first_level = rmfield(oatin.first_level,'mni_coords'); catch, end; % num_coords x 3

% HMM settings
try oat.first_level.hmm_do_glm_statewise = oatin.first_level.hmm_do_glm_statewise; oatin.first_level = rmfield(oatin.first_level,'hmm_do_glm_statewise'); catch, oat.first_level.hmm_do_glm_statewise = 0; end;

% TF settings
% general
try oat.first_level.tf_method=oatin.first_level.tf_method;  oatin.first_level = rmfield(oatin.first_level,'tf_method'); catch, oat.first_level.tf_method='none'; end;  % % available TF methods, are: 'none', 'hilbert', 'morlet', 'hanning'.  
try oat.first_level.tf_freq_range=oatin.first_level.tf_freq_range; oatin.first_level = rmfield(oatin.first_level,'tf_freq_range'); catch, disp('Using oat.source_recon.freq_range for oat.first_level.tf_freq_range'); oat.first_level.tf_freq_range=oat.source_recon.freq_range; end; % only used if doing frequency analysis, i.e. if tf_method is NOT 'none'
try oat.first_level.tf_num_freqs=oatin.first_level.tf_num_freqs; oatin.first_level = rmfield(oatin.first_level,'tf_num_freqs'); catch, oat.first_level.tf_num_freqs=1; end; % this corresponds to the num bins used for the TF transform on the source recon freq range
try oat.first_level.time_average=oatin.first_level.time_average; oatin.first_level = rmfield(oatin.first_level,'time_average');catch, oat.first_level.time_average=0;end; % all averaging is done after any TF transform but prior to calculcating GLM stats
%try oat.first_level.freq_average=oatin.first_level.freq_average; oatin.first_level = rmfield(oatin.first_level,'freq_average'); catch, oat.first_level.freq_average=0;end; % all averaging is done after any TF transform but prior to calculcating GLM stats
try oat.first_level.tf_multitaper_ncycles = oatin.first_level.tf_multitaper_ncycles; oatin.first_level = rmfield(oatin.first_level,'tf_multitaper_ncycles'); catch, oat.first_level.tf_multitaper_ncycles = []; end
if strcmp(oat.first_level.tf_method,'none') && oat.first_level.tf_num_freqs>1,
    warning('oat.first_level.tf_method is none, forcing oat.first_level.tf_num_freqs=1');
    oat.first_level.tf_num_freqs=1;
end;

% downsampling settings
% post_movingaverage_downsample_factor is used in source space after recon weights, abs, and moving average are applied 
% post_tf_downsample_factor is used after TF transform in sensor space on complex values
try oat.first_level.post_movingaverage_downsample_factor = oatin.first_level.post_movingaverage_downsample_factor; oatin.first_level = rmfield(oatin.first_level,'post_movingaverage_downsample_factor'); catch, oat.first_level.post_movingaverage_downsample_factor = []; end % divisive factor

if isfield(oatin.first_level,'post_tf_downsample_factor') && isfield(oatin.first_level,'tf_time_step')
    error('Please set either post_tf_downsample_factor or tf_time_step but not both: they are different ways of controlling the same variable')
elseif isfield(oatin.first_level,'post_tf_downsample_factor')
    oat.first_level.post_tf_downsample_factor = oatin.first_level.post_tf_downsample_factor;
    oatin.first_level = rmfield(oatin.first_level,'post_tf_downsample_factor');

elseif isfield(oatin.first_level,'tf_time_step')
    D=spm_eeg_load(Ds{1});
    fs = D.fsample;
    oat.first_level.post_tf_downsample_factor = oatin.first_level.tf_time_step*fs;
    oatin.first_level = rmfield(oatin.first_level,'tf_time_step');
else
    oat.first_level.post_tf_downsample_factor = [];
    if strcmp(oat.first_level.tf_method,'hanning');
        error('oat.first_level.tf_time_step or oat.first_level.post_tf_downsample_factor not set!!  This is required for Hanning.');
    end
end

try oat.first_level.oat.first_level.sessions_to_dotf_morlet_factor = oatin.first_level.morlet_factor; oatin.first_level = rmfield(oatin.first_level,'tf_morlet_factor'); catch, oat.first_level.tf_morlet_factor=6; end;

% hilbert
try oat.first_level.tf_hilbert_freq_res=oatin.first_level.tf_hilbert_freq_res;  oatin.first_level = rmfield(oatin.first_level,'tf_hilbert_freq_res'); catch, oat.first_level.tf_hilbert_freq_res=2; end;  % Freq res (width of freq bins) to use for Hilbert transform in Hz
try oat.first_level.tf_hilbert_do_bandpass_for_single_freq=oatin.first_level.tf_hilbert_do_bandpass_for_single_freq;  oatin.first_level = rmfield(oatin.first_level,'tf_hilbert_do_bandpass_for_single_freq'); catch, oat.first_level.tf_hilbert_do_bandpass_for_single_freq=0; end;  % do not redo band pass filtering (under the assumption it has been done as part of, or before, the source recon. This helps to avoid edge effects with epoched data in particular)
try oat.first_level.tf_hilbert_freq_ranges=oatin.first_level.tf_hilbert_freq_ranges; oatin.first_level = rmfield(oatin.first_level,'tf_hilbert_freq_ranges'); catch, oat.first_level.tf_hilbert_freq_ranges=[]; end; % /this is a matrix num_freq_bands x 2 

% morlet
try oat.first_level.tf_morlet_factor=oatin.first_level.tf_morlet_factor; oatin.first_level = rmfield(oatin.first_level,'tf_morlet_factor'); catch, oat.first_level.tf_morlet_factor=6; end;

% hanning
if isfield(oatin.first_level,'tf_hanning_ncycles')
    oat.first_level.tf_hanning_ncycles  = oatin.first_level.tf_hanning_ncycles;
    oatin.first_level = rmfield(oatin.first_level,'tf_hanning_ncycles');
else
    if strcmp(oat.first_level.tf_method,'hanning');
        error('oat.first_level.tf_hanning_ncycles not set!!');
    end
    oat.first_level.tf_hanning_ncycles = [];
end % the time step to move the tf window over

% cope type
try
    oat.first_level.cope_type=oatin.first_level.cope_type; % cope type to input to group GLM (from first level), set to 'coape', 'cope', or 'acope'
    oatin.first_level = rmfield(oatin.first_level,'cope_type');
catch,
    if(strcmp(oat.first_level.tf_method,'none')),
        oat.first_level.cope_type='acope';
    else,
        oat.first_level.cope_type='cope';
    end;
end;

% flags to output trialwise data
try oat.first_level.save_trialwise_data = oatin.first_level.save_trialwise_data; oatin.first_level = rmfield(oatin.first_level,'save_trialwise_data'); catch; oat.first_level.save_trialwise_data = 0; end
try oat.first_level.trialwise_directory = oatin.first_level.trialwise_directory; oatin.first_level = rmfield(oatin.first_level,'trialwise_directory'); catch; oat.first_level.trialwise_directory = []; end

% parcellation settings:
if ~isfield(oatin.first_level,'parcellation')
    oatin.first_level.parcellation=struct;
end;
oat.first_level.parcellation=struct; % args to pass to osl_apply_parcellation

try oat.first_level.parcellation.do = oatin.first_level.parcellation.do; oatin.first_level.parcellation = rmfield(oatin.first_level.parcellation,'do'); catch; oat.first_level.parcellation.do = 0; end
try oat.first_level.parcellation.parcellation = oatin.first_level.parcellation.parcellation; oatin.first_level.parcellation = rmfield(oatin.first_level.parcellation,'parcellation'); catch; oat.first_level.parcellation.do = 0; end
try oat.first_level.parcellation.orthogonalisation = oatin.first_level.parcellation.orthogonalisation; oatin.first_level.parcellation = rmfield(oatin.first_level.parcellation,'orthogonalisation'); catch; oat.first_level.parcellation.orthogonalisation = 'symmetric'; end
try oat.first_level.parcellation.method = oatin.first_level.parcellation.method; oatin.first_level.parcellation = rmfield(oatin.first_level.parcellation,'method'); catch; oat.first_level.parcellation.method = 'spatialBasis'; end
try oat.first_level.parcellation.normalise_voxeldata = oatin.first_level.parcellation.normalise_voxeldata; oatin.first_level.parcellation = rmfield(oatin.first_level.parcellation,'normalise_voxeldata'); catch; oat.first_level.parcellation.normalise_voxeldata = 0; end

%% establish if we are dealing with epoched or continuous data (this will be
% used to call the correct osl_run_first_level* fn
oat.first_level.is_epoched=(~isempty(oat.source_recon.D_epoched) | ~isempty(oat.source_recon.epochinfo));
try oatin.first_level=rmfield(oatin.first_level,'is_epoched'); catch, end;
if(~oat.first_level.is_epoched)
    disp('No oat.source_recon.D_epoched set. OAT will do a continuous data time-wise GLM');
else
    disp('oat.source_recon.D_epoched set. OAT will do an epoched data trial-wise GLM');
end;

try oat.first_level.do_glm_demean = oatin.first_level.do_glm_demean; oatin.first_level = rmfield(oatin.first_level,'do_glm_demean'); catch, oat.first_level.do_glm_demean = ~oat.first_level.is_epoched; end % demean GLM (data and design matrix)

% report settings
if ~isfield(oatin.first_level,'report')
    oatin.first_level.report=struct;
end;
oat.first_level.report=struct;

try oat.first_level.report.first_level_cons_to_do=oatin.first_level.report.first_level_cons_to_do; oatin.first_level.report = rmfield(oatin.first_level.report,'first_level_cons_to_do'); 
catch, 
     if(isfield(oat.first_level,'contrast')),
         oat.first_level.report.first_level_cons_to_do=1:length(oat.first_level.contrast); 
     end;
end; % Contrasts to be plotted in diagnostic report (first one listed will be the one that determines the max voxel and time point.
try oat.first_level.report.modality_to_do=oatin.first_level.report.modality_to_do; oatin.first_level.report = rmfield(oatin.first_level.report,'modality_to_do'); catch, oat.first_level.report.modality_to_do=oat.source_recon.modalities{1}; end; % Indicates which modality to do diagnostic report - only relevant for a sensor space analysis
try oat.first_level.report.time_range=oatin.first_level.report.time_range; oatin.first_level.report = rmfield(oatin.first_level.report,'time_range'); catch, oat.first_level.report.time_range=oat.first_level.time_range; end; % Indicates which modality to do diagnostic report - only relevant for a sensor space analysis
try oat.first_level.report.freq_range=oatin.first_level.report.freq_range; oatin.first_level.report = rmfield(oatin.first_level.report,'freq_range'); catch, oat.first_level.report.freq_range=oat.first_level.tf_freq_range; end; % Indicates which modality to do diagnostic report - only relevant for a sensor space analysis
try oat.first_level.report.max_method=oatin.first_level.report.max_method; oatin.first_level.report = rmfield(oatin.first_level.report,'max_method'); catch, oat.first_level.report.max_method='max_abs_tstat'; end; % method to use to find location, time and freq of interest for report: can be 'max_abs_tstat' or 'max_tstat'

try oat.first_level.warp_fields=oatin.first_level.warp_fields; oatin.first_level = rmfield(oatin.first_level,'warp_fields'); catch, oat.first_level.warp_fields=[]; end; % warp fields

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subjects stats setting

if ~isfield(oatin,'subject_level')
    oatin.subject_level=struct;
end;

oat.subject_level=struct;

% no required settings

% optional settings:
try oat.subject_level.session_index_list=oatin.subject_level.session_index_list; oatin.subject_level = rmfield(oatin.subject_level,'session_index_list');catch, for ii=1:num_sessions,oat.subject_level.session_index_list{ii}=ii; end; end; % This is a vector of cells of length num_subjects, with each cell being a vector listing which sessions belong to each subject. Needs to be of the correct size to include all subjects even if oat.subject_level.subjects_to_do is being used to select a subset of subjects.
num_subjects=length(oat.subject_level.session_index_list);
try oat.subject_level.subjects_to_do=oatin.subject_level.subjects_to_do; oatin.subject_level = rmfield(oatin.subject_level,'subjects_to_do'); catch, oat.subject_level.subjects_to_do=oat.first_level.sessions_to_do; end; % list of subject indexes to run on.
try oat.subject_level.name=oatin.subject_level.name; oatin.subject_level = rmfield(oatin.subject_level,'name'); catch, oat.subject_level.name='sub_level'; end;

try oat.subject_level.compute_laterality=oatin.subject_level.compute_laterality;
    oatin.subject_level = rmfield(oatin.subject_level,'compute_laterality');
catch
    oat.subject_level.compute_laterality = [];
end

try oat.subject_level.merge_contraipsi=oatin.subject_level.merge_contraipsi;
    oatin.subject_level = rmfield(oatin.subject_level,'merge_contraipsi');
catch
    oat.subject_level.merge_contraipsi = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% group stats setting

if ~isfield(oatin,'group_level')
    oatin.group_level=struct;
end;

oat.group_level=struct;

% no required settings

% optional settings:
try oat.group_level.subjects_to_do=oatin.group_level.subjects_to_do; oatin.group_level = rmfield(oatin.group_level,'subjects_to_do'); catch, oat.group_level.subjects_to_do=oat.subject_level.subjects_to_do; end; % list of subject indexes to run on.
num_subjects=length(oat.group_level.subjects_to_do);
try oat.group_level.group_design_matrix=oatin.group_level.group_design_matrix; oatin.group_level = rmfield(oatin.group_level,'group_design_matrix');catch, oat.group_level.group_design_matrix=ones(1,num_subjects); end; % Num_regressors x Num_subjects (where num_subjects needs to be length(subjects_to_do))
try oat.group_level.group_contrast=oatin.group_level.group_contrast; oatin.group_level = rmfield(oatin.group_level,'group_contrast');catch, oat.group_level.group_contrast{1}=[1]; end; % Each contrast should be [num_regressors x 1], e.g. group_contrast{1}=[1 -1]';
try oat.group_level.group_contrast_name=oatin.group_level.group_contrast_name; oatin.group_level = rmfield(oatin.group_level,'group_contrast_name');catch, if(isfield(oat.group_level,'group_contrast')),for c=1:length(oat.group_level.group_contrast), oat.group_level.group_contrast_name{c}=['gc' num2str(c)];end;end; end; % Each contrast can be given a name
try oat.group_level.glm_method=oatin.group_level.glm_method; oatin.group_level = rmfield(oatin.group_level,'glm_method');catch, oat.group_level.glm_method='ols'; end; % ols or fixed_effects
try oat.group_level.freq_range=oatin.group_level.freq_range; oatin.group_level = rmfield(oatin.group_level,'freq_range');catch, oat.group_level.freq_range=[]; end;
try oat.group_level.time_range=oatin.group_level.time_range; oatin.group_level = rmfield(oatin.group_level,'time_range');catch, oat.group_level.time_range=[]; end;
try oat.group_level.name=oatin.group_level.name; oatin.group_level = rmfield(oatin.group_level,'name'); catch, oat.group_level.name='group_level'; end;
try oat.group_level.spatial_smooth_fwhm=oatin.group_level.spatial_smooth_fwhm; oatin.group_level = rmfield(oatin.group_level,'spatial_smooth_fwhm'); catch,  oat.group_level.spatial_smooth_fwhm=0; end; % Gaussian spatial smoothing to apply, via std dev in mm (set to 0 for none). Note that this is applied prior to applying the group level mask
try oat.group_level.time_smooth_std=oatin.group_level.time_smooth_std; oatin.group_level = rmfield(oatin.group_level,'time_smooth_std'); catch, oat.group_level.time_smooth_std=0; end;% Gaussian temporal smoothing to apply, via std dev in secs (set to 0 for none). Note that this is applied prior to applying the group level time window
try oat.group_level.space_average=oatin.group_level.space_average; oatin.group_level = rmfield(oatin.group_level,'space_average'); catch, oat.group_level.space_average=0; end; % all averaging is done prior to calculcating GLM stats (i.e. on the GLM "data")
try oat.group_level.freq_average=oatin.group_level.freq_average; oatin.group_level = rmfield(oatin.group_level,'freq_average'); catch, oat.group_level.freq_average=0; end; % all averaging is done prior to calculcating GLM stats (i.e. on the GLM "data")
try oat.group_level.time_average=oatin.group_level.time_average; oatin.group_level = rmfield(oatin.group_level,'time_average'); catch, oat.group_level.time_average=0; end; % all averaging is done prior to calculcating GLM stats (i.e. on the GLM "data")
try oat.group_level.use_tstat=oatin.group_level.use_tstat; oatin.group_level = rmfield(oatin.group_level,'use_tstat'); catch, oat.group_level.use_tstat=0; end; % if set to 1 uses tstats as input to group GLM (from first level), if set to 0 uses cope
try oat.group_level.use_robust_glm=oatin.group_level.use_robust_glm; oatin.group_level = rmfield(oatin.group_level,'use_robust_glm'); catch, oat.group_level.use_robust_glm=0; end; % do robust GLM (uses bisquare via Matlab robustfit fn)
try oat.group_level.mask_fname=oatin.group_level.mask_fname; oatin.group_level = rmfield(oatin.group_level,'mask_fname'); catch, end;
try oat.group_level.mni_coords=oatin.group_level.mni_coords; oatin.group_level = rmfield(oatin.group_level,'mni_coords'); catch, end;
try oat.group_level.group_varcope_time_smooth_std=oatin.group_level.group_varcope_time_smooth_std; oatin.group_level = rmfield(oatin.group_level,'group_varcope_time_smooth_std'); catch, oat.group_level.group_varcope_time_smooth_std=0.04; end;
try oat.group_level.group_varcope_spatial_smooth_fwhm=oatin.group_level.group_varcope_spatial_smooth_fwhm; oatin.group_level = rmfield(oatin.group_level,'group_varcope_spatial_smooth_fwhm'); catch, oat.group_level.group_varcope_spatial_smooth_fwhm=12; end;
try oat.group_level.store_lower_level_copes=oatin.group_level.store_lower_level_copes; oatin.group_level = rmfield(oatin.group_level,'store_lower_level_copes'); catch, oat.group_level.store_lower_level_copes=0; end;
try oat.group_level.first_level_contrasts_to_do=oatin.group_level.first_level_contrasts_to_do; oatin.group_level = rmfield(oatin.group_level,'first_level_contrasts_to_do'); catch, oat.group_level.first_level_contrasts_to_do=1:length(oat.first_level.contrast); end;

% report settings
if ~isfield(oatin.group_level,'report')
    oatin.group_level.report=struct;
end;
oat.group_level.report=struct;
try oat.group_level.report.group_level_cons_to_do=oatin.group_level.report.group_level_cons_to_do; oatin.group_level.report = rmfield(oatin.group_level.report,'group_level_cons_to_do'); 
catch, 
     if(isfield(oat.group_level,'group_contrast')),
         oat.group_level.report.group_level_cons_to_do=1:length(oat.group_level.group_contrast); 
     end;
end; % Contrasts to be plotted in diagnostic report (first one listed will be the one that determines the max voxel and time point.
try oat.group_level.report.first_level_cons_to_do=oatin.group_level.report.first_level_cons_to_do; oatin.group_level.report = rmfield(oatin.group_level.report,'first_level_cons_to_do'); 
catch, 
     if(isfield(oat.first_level,'contrast')),
         oat.group_level.report.first_level_cons_to_do=1:length(oat.first_level.contrast); 
     end;
end; % Contrasts to be plotted in diagnostic report (first one listed will be the one that determines the max voxel and time point.
try oat.group_level.report.modality_to_do=oatin.group_level.report.modality_to_do; oatin.group_level.report = rmfield(oatin.group_level.report,'modality_to_do'); catch, oat.group_level.report.modality_to_do=oat.source_recon.modalities{1}; end; % Indicates which modality to do diagnostic report - only relevant for a sensor space analysis
try oat.group_level.report.time_range=oatin.group_level.report.time_range; oatin.group_level.report = rmfield(oatin.group_level.report,'time_range'); catch, end; % Indicates which modality to do diagnostic report - only relevant for a sensor space analysis
try oat.group_level.report.freq_range=oatin.group_level.report.freq_range; oatin.group_level.report = rmfield(oatin.group_level.report,'freq_range'); catch, end; % Indicates which modality to do diagnostic report - only relevant for a sensor space analysis
try oat.group_level.report.max_method=oatin.group_level.report.max_method; oatin.group_level.report = rmfield(oatin.group_level.report,'max_method'); catch, oat.group_level.report.max_method='max_abs_tstat'; end; % method to use to find location, time and freq of interest for report: can be 'max_abs_tstat' or 'max_tstat'
try oat.group_level.report.show_lower_level_copes=oatin.group_level.report.show_lower_level_copes; oatin.group_level.report = rmfield(oatin.group_level.report,'show_lower_level_copes'); catch, oat.group_level.report.show_lower_level_copes=1; end; % flag to indicate whether lower level copes (the group GLM input) time courses should be plotted
if(oat.group_level.report.show_lower_level_copes), oat.group_level.store_lower_level_copes=1;end;

try oat.group_level.report.show_lower_level_cope_maps=oatin.group_level.report.show_lower_level_cope_maps; oatin.group_level.report = rmfield(oatin.group_level.report,'show_lower_level_cope_maps'); catch, oat.group_level.report.show_lower_level_cope_maps=0; end; % flag to indicate whether lower level copes (the group GLM input) spatial maps should be plotted


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% copy any results
try oat.source_recon.results_fnames=oatin.source_recon.results_fnames;
    oatin.source_recon = rmfield(oatin.source_recon,'results_fnames'); catch, end;
try oat.first_level.results_fnames=oatin.first_level.results_fnames;
    oatin.first_level = rmfield(oatin.first_level,'results_fnames'); catch, end;
try oat.subject_level.results_fnames=oatin.subject_level.results_fnames;
    oatin.subject_level = rmfield(oatin.subject_level,'results_fnames'); catch, end;
try oat.group_level.results_fnames=oatin.group_level.results_fnames;
    oatin.group_level = rmfield(oatin.group_level,'results_fnames'); catch, end;
try oat.results=oatin.results; oatin = rmfield(oatin,'results'); catch, end;
try oat.source_recon.results=oatin.source_recon.results;
    oatin.source_recon = rmfield(oatin.source_recon,'results'); catch, end;
try oat.first_level.results=oatin.first_level.results;
    oatin.first_level = rmfield(oatin.first_level,'results'); catch, end;
try oat.subject_level.results=oatin.subject_level.results;
    oatin.subject_level = rmfield(oatin.subject_level,'results'); catch, end;
try oat.group_level.results=oatin.group_level.results;
    oatin.group_level = rmfield(oatin.group_level,'results'); catch, end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check people haven't set any weird fields

weirdfields = fieldnames(oatin.source_recon.report);
if ~isempty(weirdfields)
    disp('The following oat.source_recon.report settings were not recognized by osl_check_oat');
    
    for iprint = 1:numel(weirdfields)
        disp([' ' weirdfields{iprint} ' '])
    end
    error('Invalid oat.source_recon.report settings');
end % if ~isempty(weirdfields)

weirdfields = fieldnames(oatin.first_level.report);
if ~isempty(weirdfields)
    disp('The following oat.first_level.report settings were not recognized by osl_check_oat');
    
    for iprint = 1:numel(weirdfields)
        disp([' ' weirdfields{iprint} ' '])
    end
    error('Invalid oat.first_level.report settings');
end % if ~isempty(weirdfields)

weirdfields = fieldnames(oatin.group_level.report);
if ~isempty(weirdfields)
    disp('The following oat.group_level.report settings were not recognized by osl_check_oat');
    
    for iprint = 1:numel(weirdfields)
        disp([' ' weirdfields{iprint} ' '])
    end
    error('Invalid oat.group_level.report settings');
end % if ~isempty(weirdfields)

oatin.source_recon = rmfield(oatin.source_recon,'report');
oatin.first_level = rmfield(oatin.first_level,'report');
oatin.group_level = rmfield(oatin.group_level,'report');

weirdfields = fieldnames(oatin.first_level.parcellation);
if ~isempty(weirdfields)
    disp('The following oat.first_level.parcellation settings were not recognized by osl_check_oat');
    
    for iprint = 1:numel(weirdfields)
        disp([' ' weirdfields{iprint} ' '])
    end
    error('Invalid oat.first_level.parcellation settings');
end % if ~isempty(weirdfields)
oatin.first_level = rmfield(oatin.first_level,'parcellation');

weirdfields = fieldnames(oatin.source_recon);
if ~isempty(weirdfields)
    disp('The following oat.source_recon settings were not recognized by osl_check_oat');
    
    for iprint = 1:numel(weirdfields)
        disp([' ' weirdfields{iprint} ' '])
    end
    error('Invalid oat.source_recon settings');
end % if ~isempty(weirdfields)

weirdfields = fieldnames(oatin.first_level);
if ~isempty(weirdfields)
    disp('The following oat.first_level settings were not recognized by osl_check_oat');
    
    for iprint = 1:numel(weirdfields)
        disp([' ' weirdfields{iprint} ' '])
    end
    error('Invalid oat.first_level settings');
end % if ~isempty(weirdfields)

weirdfields = fieldnames(oatin.subject_level);
if ~isempty(weirdfields)
    disp('The following oat.subject_level settings were not recognized by osl_check_oat');
    
    for iprint = 1:numel(weirdfields)
        disp([' ' weirdfields{iprint} ' '])
    end
    error('Invalid oat.subject_level settings');
end % if ~isempty(weirdfields)

weirdfields = fieldnames(oatin.group_level);
if ~isempty(weirdfields)
    disp('The following oat.group_level settings were not recognized by osl_check_oat');
    
    for iprint = 1:numel(weirdfields)
        disp([' ' weirdfields{iprint} ' '])
    end
    error('Invalid oat.group_level settings');
end % if ~isempty(weirdfields)

try oatin = rmfield(oatin,'osl2_version');catch, end;
try, oatin = rmfield(oatin,'date');catch, end;
oatin = rmfield(oatin,'source_recon');
oatin = rmfield(oatin,'first_level');
oatin = rmfield(oatin,'subject_level');
oatin = rmfield(oatin,'group_level');
try oatin = rmfield(oatin,'fname'); catch, end;

weirdfields = fieldnames(oatin);
if ~isempty(weirdfields)
    disp('The following oat settings were not recognized by osl_check_oat');
    
    for iprint = 1:numel(weirdfields)
        disp([' ' weirdfields{iprint} ' '])
    end
    error('Invalid osl_check_oat settings');
end % if ~isempty(weirdfields)

%% add osl version

oat.osl2_version=osl2_version;
