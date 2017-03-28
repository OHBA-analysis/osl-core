%% PRACTICAL: OAT CONTINUOUS BEAMFORMER CONTRAST
%
% In this practical we will estimate neural activity in the brain's source
% space using a beamformer algorithm, task dependant differences in this source
% activity is then quantified using a GLM. This will go through the following
% steps:
% 
%   1) Set-up an OAT analysis: source_recon and first_level
%   2) Run source space GLM fitting and contrasts
%   3) Check coregistration
%   4) Save and view t-stat volumes
%   5) GLM analysis in a ROI
%   6) Time-frequency analysis in an ROI
%   7) Whole brain time-frequency analysis
%
% This practical we will work with a single subject's data from a blocked
% finger tapping experiment (CTF data courtesy of Matt Brookes), and
% perform a time-frequency analysis in the beta band using a time-wise GLM
% in source space.   
 
% Data can be downloaded from: 
% www.fmrib.ox.ac.uk/~woolrich/ft_subject1_data.tar.gz

%% SETUP THE MATLAB PATHS
%
% Sets the Matlab paths to include OSL. Change these paths so that they 
% correspond to the setup on your computer. You will also need to ensure 
% that fieldtrip and spm are not in your matlab path (as they are included
% within OSL).%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

OSLDIR = getenv('OSLDIR');
    
% This cell sets the Matlab paths to include OSL. Change the osldir path so
% that it corresponds to the setup on your computer before running the cell. 


%% INITIALISE GLOBAL SETTINGS FOR THIS ANALYSIS
%
% This cell sets the directory that OAT will work in. Change the workingdir
% variable to correspond to the correct directory on your computer before
% running the cell.

% directory where the data is:
%datadir=['/Users/abaker/Data/']; % directory where the data is
%workingdir=[datadir '/ctf_fingertap_subject1_data_osl2']; % this is the directory where the SPM files will be stored in
workingdir=[practical_dir '/homedir/vols_data/osl_testdata/osl2_testdata_dir/ctf_fingertap_subject1_data_osl2']; % directory where the data is

cmd = ['mkdir ' workingdir]; if ~exist(workingdir, 'dir'), unix(cmd); end % make dir to put the results in

cd(workingdir);
clear spm_files;

%% SET UP THE LIST OF SUBJECTS FOR THE ANALYSIS
%
% Specify a list of the fif files, structual files (not applicable for this
% practical) and SPM files (which will be created). It is important to make
% sure that the order of these lists is consistent across sessions. 
% Note that here we only have 1 subject, but more generally there would be
% more than one, e.g.:
%
% fif_files{1}=[testdir '/fifs/sub1_face_sss.fif']; 
% fif_files{2}=[testdir '/fifs/sub2_face_sss.fif']; 
% etc...
%
% spm_files{1} = [workingdir '/sub1_face_sss.mat'];
% spm_files{2} = [workingdir '/sub2_face_sss.mat'];
% etc...
% set up a list of SPM MEEG object file names (we only have one here)
% set up a list of SPM MEEG object file names (we only have one here)
% set up a list of SPM MEEG object file names (we only have one here)
%
% As we will be working in source space for this practical we will also define
% a list of structural MRI files, one for each participant. These will be used
% during coregistration and source localisation to define the source space for
% each participant
spm_files={[workingdir '/dsubject1.mat']};
         
%structural_files = {[workingdir '/subject1_struct.nii']};      
structural_files = {[workingdir '/subject1_struct_canon.nii']};      

pos_files = {[workingdir '/subject1_Motorcon.pos']};      


cleanup_files=0; % flag to indicate that you want to clean up files that are no longer needed


%% Sort out fiducials:
for f = 1:length(spm_files)
    
    spm_file = prefix(spm_files{f},'');
    D = spm_eeg_load(spm_file);
    
    fID = fopen(pos_files{f});
    fid_data = textscan(fID, '%s %f %f %f');
    fclose(fID);
    
    fid_new = [];
    
    % Fiducials:
    fid_inds = [find(strcmpi(fid_data{1},'nasion'))
        find(strcmpi(fid_data{1},'left'))
        find(strcmpi(fid_data{1},'right'))];
    
    fid_new.fid.pnt = [fid_data{2}(fid_inds) fid_data{3}(fid_inds) fid_data{4}(fid_inds)] * 10;
    fid_new.fid.label = {'nas';'lpa';'rpa'};
    
    % Headshape:
    hs_inds = setdiff(2:length(fid_data{1}),fid_inds);
    fid_new.pnt = [fid_data{2}(hs_inds) fid_data{3}(hs_inds) fid_data{4}(hs_inds)] * 10;   
    fid_new.unit = 'mm';
    
    % Labels:
    fid_labels = fid_data{1};
        
    D = fiducials(D,fid_new);
    D.save;
    
end


%% Run coregistration & forward model
for f = 1:length(spm_files)
    
    S = [];
    S.D             = prefix(spm_files{f},'');
    S.mri           = structural_files{f};
    S.fid.label     = struct('nasion','nas','lpa','lpa','rpa','rpa');
    S.useheadshape  = 1;
    S.use_rhino     = 1;
    S.useCTFhack    = 0;
    S.forward_meg   = 'Single Shell';
    osl_headmodel(S);
    
end

%% Check & correct coregistration:
for f = 1:length(spm_files)
    D = spm_eeg_load(prefix(spm_files{f},''));
    %rhino_manual(D)
    rhino_display(D)
    %waitfor(gcf)
end

%% LOAD AND VIEW RAW DATA
%
% The data we are using in this analysis in contained in an spm object, this
% can be loaded using spm_eeg_load. Take a look at the contents by typing "D"
% into the command line in matlab. Note that this is continuous data with
% 232000 timepoints sampled at 250Hz
%
% Next look at the data in oslview, Note that there is only one sensor type (as
% this is CTF data) and that there are no obvious artefacts to mark.

subnum = 1;             
D = spm_eeg_load(spm_files{subnum});


% check data in oslview
oslview(D);
       
%% SET UP WHOLE BRAIN BEAMFORMER IN OAT
%
% In previous practicals, the oat.source_recon parameters defined some simple
% data preparation to be done prior to sensor space analysis. In this practical
% we will define options for a whole brain source reconstruction  which
% requires some additional options.
%
% The conditions, freq range and time range are the same as we have used in
% previous analyses. The following parameters are used for the source analysis:
%   'method': This defines the source reconstruction method to be used. other options include 'beamform_bilateral' and 'mne_datacov'
%   'normalise_method': This defines how the data will be normalised
%   'gridstep': This is the distance (in mm) between points to be reconstructed, the spatial resolution of the analysis
%   'forward_meg': This is the forward model to be used
%   'mri': The path to the structural volume
%   'dirname': The output path
%
% Several other options exist to perform the source localisation based on
% covariance estimates from within HMM states defined from the sensors (see
% Woolrich et al 2013 NeuroImage, Baker et al 2014 elife) We will not be
% performing state-by-state analyses in this practical so do_hmm is set to "0"


oat=[];
oat.source_recon.D_continuous=spm_files;
oat.source_recon.conditions={'Undefined'};
oat.source_recon.freq_range=[4 30]; % frequency range in Hz
oat.source_recon.time_range=[300,32*30];

% beamformer parameters.
oat.source_recon.method='beamform';
oat.source_recon.gridstep=8; % in mm, using a lower resolution here than you would normally, for computational speed
do_hmm=0;
if(do_hmm)
    oat.source_recon.hmm_num_states=13;    
    oat.source_recon.hmm_num_starts=1;
    oat.first_level.hmm_do_glm_statewise=0;
end;

oat.source_recon.dirname=[workingdir '/subj1_results_beta_rhino_hmm' num2str(do_hmm) '_' oat.source_recon.method];

oat.source_recon.normalise_method='none';

oat = osl_check_oat(oat);

if(0)
    % load in previously run HMM to avoid having to rerun it
    oat=osl_load_oat(oat);
    recon=oat_load_results(oat,oat.source_recon.results_fnames{1});
    
    oat.source_recon.hmm_block={};
    oat.source_recon.hmm_block{1}=recon.block;
end;

%% CHECK OAT SOURCE_RECON AND RUN BEAMFORMER
%
% osl_check_oat will have added some extra parameters to your OAT structure,
% take a look by printing the oat and oat.source_recon structures
%
% Next, we want to run the beamformer. We define oat.to_do such that only the
% source_recon stage will be completed and define an output directory name
% before calling osl_run_oat.

oat
oat.source_recon

oat.to_do=[1 0 0 0];
oat.source_recon.pca_dim=250;

oat = osl_run_oat(oat);

%% ESTABLISH REGRESSOR FOR CONTINUOUS TIME GLM.
%
% This section generates the predictors for our GLM analysis based on the order
% of the finger tapping blocks. There were four conditions:
% Left hand tap
% Right hand tap
% Rest
% Both hands tap
% Rest at start
%
% The order of these blocks in recorded in block_order and their length in seconds
% in block_length.
%
% This should be setup to correspond to the time window
% over which the source recon is carried out.

D=spm_eeg_load(spm_files{1});

% This should be setup to correspond to the same time window
% as the full time window for D
x=zeros(length(D.time),5);

block_length=30; %s
block_order=[5 5 5 5 5 5 5 5 5 5 4 3 2 1 2 3 1 4 3 4 1 3 2 1 4 4 2 1 3 3 4 1 4 3 1 2 1 2 3 4 3 4 1 2 3 4 1 2];

% [Left, Right, Rest, Both, Rest_at_start]
% [  1     2      4     8     16 ]
% figure;plot(D.time,squeeze(D(1,:,:)))
% emacs ~/vols_data/From_Nottingham_with_Love/JRH_MotorCon_20100429_01_FORMARK.ds/MarkerFile.mrk 

tres=1/(D.fsample);
tim=1;
for tt=1:length(block_order),    
    x(tim:tim+block_length/tres-1,block_order(tt))=1;
    tim=tim+block_length/tres;
end;

figure;plot(D.time,x);

%% SET UP FIRST LEVEL GLM ANALYSIS CONTRASTS
%
% This cell defines the contrast vectors for the first level analysis. The OAT
% structure created by the source_recon we ran earlier is loaded in prior to
% defining the contrasts.
%
% Each contrast is a vector containing a weight per condition defining how
% the condition parameter estimates are to be compared. Each vector will
% produce a different t-map across the sensors. Contrasts 1 and 2 describe
% positive correlations between each sensor's activity and the presence
% of a motorbike or face stimulus respectively. Contrast 3 tests whether
% each sensors activity is larger for faces than motorbikes.
%
% Several other parameters are set afterwards defining the time-frequency
% methods to be used in the first-level analysis, downsampling factor and
% the name of this first-level analysis.
% Note, to look for changes in power specifically in the beta band 
% (since the motor system is known to produce beta band changes) 
% we have set "oat.first_level.tf_freq_range=[13 30];"


oat.source_recon.dirname=[workingdir '/subj1_results_beta_rhino_hmm' num2str(do_hmm) '_' oat.source_recon.method];
oat=osl_load_oat(oat.source_recon.dirname, oat.first_level.name,'sub_level','group_level');

oat.first_level.name=['first_level'];

% GLM stuff:
oat.first_level.design_matrix=x';
oat.first_level.contrast{1}=[-1 0 1 0 0]'; % rest-left
oat.first_level.contrast{2}=[0 -1 1 0 0]'; % rest-right
oat.first_level.contrast{3}=[0  0 1 -1 0]'; % rest-both
oat.first_level.contrast_name{1}='rest-left';
oat.first_level.contrast_name{2}='rest-right';
oat.first_level.contrast_name{3}='rest-both';
oat.first_level.bc=[0 0 0];
oat.first_level.report.first_level_cons_to_do=1:3;
oat.first_level.cope_type='cope';
oat.first_level.tf_hilbert_freq_res=diff(oat.first_level.tf_freq_range);
oat.first_level.tf_method='hilbert';
oat.first_level.post_tf_downsample_factor=1; 
oat.first_level.tf_hilbert_freq_ranges=[13 30];
oat.first_level.time_moving_av_win_size=[];
oat.first_level.tf_hilbert_do_bandpass_for_single_freq=0;
oat.first_level.do_weights_normalisation=1;
oat.first_level.hmm_do_glm_statewise=0; 
oat.first_level.do_glm_demean=0;

oat.first_level.parcellation.do=1;
if oat.first_level.parcellation.do
    tilde='/Users/woolrich/';
    addpath([tilde 'Dropbox/vols_scripts/MEG-ROI-nets']);

    parc_file=[tilde '/homedir/vols_data/hmm_investigations/parcellations/fmri_d100_parcellation_with_PCC_reduced_2mm'];
    parcellationfile = [parc_file '_ss5mm_ds8mm'];
    
    parcellationfile = [tilde '/homedir/vols_data/hmm_investigations/parcellations/aal2mni_cortical_4d_8mm'];
    parcellationfile = [tilde '/homedir/vols_data/hmm_investigations/parcellations/AAL_4d_8mm'];
    
    oat.first_level.parcellation.parcellation=parcellationfile;
    oat.first_level.parcellation.orthogonalisation = 'symmetric';
    oat.first_level.parcellation.method            = 'spatialBasis';
    oat.first_level.parcellation.normalise_voxeldata = 0;
end
oat.first_level.name=[oat.first_level.name '_parc' num2str(oat.first_level.parcellation.do)];

%% RUN FIRST-LEVEL CONTINUOUS TIME OAT ANALYSIS
%
% Run OAT analysis.

oat.to_do=[1 1 0 0];

oat = osl_check_oat(oat);
oat = osl_run_oat(oat);

%% SAVE ADDITIONAL NIFTI VOLUMES AND VIEW IN FSLVIEW
%
% Additional Nifti volumes can be saved out by calling oat_save_nii_stats and
% viewed by calling fslview. In this case we are saving out first level
% contrasts 1,2 and 3 in all frequencies. fslview is then called to view the
% Nifti volume in the first frequency

S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.first_level_contrasts=[3]; % list of first level contrasts to output
S2.resamp_gridstep=oat.source_recon.gridstep;

clear statsdir;
for ff=1:size(oat.first_level.tf_hilbert_freq_ranges,1),
    S2.freq_bin=ff;
    S2.stats_dir=[oat.source_recon.dirname '/' oat.first_level.name 'new_f' num2str(mean(oat.first_level.tf_hilbert_freq_ranges(ff,:))) '_stats_dir'];

    [statsdir{ff},times]=oat_save_nii_stats(S2);    
end;

% make sure you view results using fslview
contrast_num=3;
runcmd(['fslview ' statsdir{1} '/cope' num2str(contrast_num) '_2mm ' statsdir{1} '/tstat' num2str(contrast_num) '_2mm &']);

%% CREATE AN ROI MASK
% 
% Often we will want to view the GLM results from a specific location, perhaps
% defined anatomically or with a localiser scan.
%
% This section creates a mask by thresholding the tstat volume from contrast 3
% in the previous OAT analysis. This is done with a call to fslmaths 

% First calculate an ROI mask 
oat.source_recon.dirname=[workingdir '/subj1_results_beta_hmm' num2str(do_hmm) '_' oat.source_recon.method];
oat.first_level.name=['first_level'];
oat=osl_load_oat(oat.source_recon.dirname, oat.first_level.name,'sub_level','group_level');
statsdir=[oat.source_recon.dirname '/' oat.first_level.name 'new_f' num2str(mean(oat.first_level.tf_hilbert_freq_ranges(ff,:))) '_stats_dir'];

% use FSL maths to threshold to create mask
con=3;
thresh=50;
runcmd(['fslmaths ' statsdir '/tstat' num2str(con) '_2mm -thr ' num2str(thresh) ' ' statsdir '/tstat' num2str(con) '_2mm_mask']);

% view the mask in fslview:
runcmd(['fslview ' statsdir '/tstat' num2str(con) '_2mm_mask &']);

%% VIEW DATA WITHIN ROI MASK
%
% The first level analysis can be adapted to output the time-series that the
% GLM would be run on without actually fitting the model. This is often useful
% to see the raw data that we are going to model. This section generates the
% average raw data within the mask and plots it against the regressor we
% created earlier.

% Look at Hilbert Envelope timecourse averaged over ROI
oat.source_recon.dirname=[workingdir '/subj1_results_beta_hmm' num2str(do_hmm) '_' oat.source_recon.method];
oat.first_level.name=['first_level'];
oat=osl_load_oat(oat.source_recon.dirname, oat.first_level.name,'sub_level','group_level');

statsdir=[oat.source_recon.dirname '/' oat.first_level.name 'new_f' num2str(mean(oat.first_level.tf_hilbert_freq_ranges(ff,:))) '_stats_dir'];
oat.first_level.mask_fname=[statsdir '/tstat' num2str(con) '_2mm_mask'];
oat.first_level.design_matrix=x';
oat.first_level.contrast{1}=[-1 0 1 0 0]'; % rest-left
oat.first_level.contrast{2}=[0 -1 1 0 0]'; % rest-right
oat.first_level.contrast{3}=[0  0 1 -1 0]'; % rest-both
oat.first_level.bc=[0 0 0];
%oat.first_level.diagnostic_cons_to_do=1:3;
oat.first_level.tf_hilbert_freq_res=diff(oat.first_level.tf_freq_range);

oat.first_level.tf_method='hilbert';
oat.first_level.post_movingaverage_downsample_factor=1; 
oat.first_level.post_tf_downsample_factor=2; 
oat.first_level.tf_hilbert_freq_ranges=[13 30];
oat.first_level.time_moving_av_win_size=1;
oat.first_level.tf_hilbert_do_bandpass_for_single_freq=0;
oat.first_level.do_glm_demean=0;
oat.first_level.name='roi';
oat.first_level.doGLM=0; % this turns the GLM off and just outputs the reconstructed data (that would be input into the GLM)

oat.to_do=[0 1 0 0];
oat = osl_run_oat(oat);

stats=oat_load_results(oat,oat.first_level.results_fnames{1});
figure;plot(stats.glm_input_times,normalise(squeeze(mean(stats.glm_input_data,1))));

% compare to design matrix:
time_ind = ismembc2(stats.glm_input_times, D.time); %use undocumented built-in function
time_ind(time_ind==0) = [];

ho;plot(D.time(time_ind),x(time_ind,1),'r','LineWidth',2);
plot(D.time(time_ind),x(time_ind,2),'g','LineWidth',2);
plot(D.time(time_ind),x(time_ind,4),'k','LineWidth',2);
legend('beta','left','right','both');
plot4paper('time(secs)','power');

%% COMPUTE TIME-FREQUENCY REPRESENTATION WITHIN ROI - 1
%
% Another use for outputting the raw data is to compute the time-frequency
% information within an ROI. In order to do this here, we first will need to
% redo the beamformer with a wider frequency band

% Load in previous oat
oat.source_recon.dirname=[workingdir '/subj1_results_beta_hmm' num2str(do_hmm) '_' oat.source_recon.method];
oat.first_level.name=['first_level'];
oat=osl_load_oat(oat.source_recon.dirname, oat.first_level.name,'sub_level','group_level');

con=3;
oat.source_recon.mask_fname=[statsdir '/tstat' num2str(con) '_2mm_mask'];

% need new OAT dirname for new source recon
oat.source_recon.dirname=[workingdir '/subj1_results_wideband.oat'];

oat.source_recon.freq_range=[4 48];
oat.to_do=[1 0 0 0];
oat = osl_run_oat(oat);

%% COMPUTE TIME-FREQUENCY REPRESENTATION WITHIN ROI - 2
%
% Once the source reconstruction has been repeated with a wider band we can use
% the first_level to generate the raw data that we would use in a GLM analysis.
% This can subsequently be plotted against the regressors we defined earlier.
% In this case we extract the raw source estimates in several frequency bands.

oat.first_level.mask_fname=[statsdir '/tstat' num2str(con) '_2mm_mask'];
oat.first_level.tf_method='morlet';
oat.first_level.tf_morlet_factor=6;
if(isfield(oat.first_level,'tf_hilbert_freq_ranges') ),oat.first_level=rmfield(oat.first_level,'tf_hilbert_freq_ranges');end;
oat.first_level.tf_num_freqs=15;
oat.first_level.tf_freq_range=[4 48];
oat.first_level.tf_hilbert_freq_res=5;
oat.first_level.doGLM=0; % does not fit GLM and will output timeseries that would have been input into GLM
oat.first_level.post_tf_downsample_factor=2; 
oat.first_level.do_glm_demean=0;
oat.first_level.time_moving_av_win_size=10;

oat.to_do=[0 1 0 0];
oat = osl_run_oat(oat);

%% COMPUTE TIME-FREQUENCY REPRESENTATION WITHIN ROI - 3
%
% Here we will load in the results of our new OAT analysis and plot the
% power across time and frequency against the regressors from the GLM fit.
%
% Notice the large increase in power between 10 and 25Hz during the rest
% condition.

stats=oat_load_results(oat,oat.first_level.results_fnames{1});

figure;imagesc(stats.glm_input_times,stats.glm_input_frequencies,squeeze(mean(stats.glm_input_data,1))');colorbar;
axis xy;

% compare to design matrix:
ho;
x3=x(:,3);
x3(x3==0)=nan;
plot(D.time(time_ind),10*x3(time_ind),'k','LineWidth',5);ho;
legend('rest');
plot4paper('beta power','time(secs)');

