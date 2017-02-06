% This practical will work with a single subject's data from an emotional
% faces experiment (Elekta Neuromag data). You can get the data from:
% https://sites.google.com/site/ohbaosl/practicals/practical-data/emotional-face-processing-elekta-neuromag-data
%
% Work your way through the script cell by cell using the supplied dataset.
% As well as following the instructions below, make sure that you read all
% of the comments (indicated by %), as these explain what each step is
% doing. Note that you can run a cell (marked by %%) using the ?Cell? drop
% down menu on the Matlab GUI.

%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

global OSLDIR;

% set this to where you have downloaded OSL and the practical data:
practical_dir='/home/mwoolrich/Desktop';
osldir=[practical_dir '/osl2'];

practical_dir='/Users/woolrich';
osldir=[practical_dir '/Dropbox/osl2'];

addpath(osldir);
osl_startup(osldir);

%%%%%%%%%%%%%%%%%%
%% INITIALISE GLOBAL SETTINGS FOR THIS ANALYSIS

% directory where the data is:
workingdir=[practical_dir '/homedir/vols_data/osl_testdata/osl2_testdata_dir/faces_subject1_data_osl2']; % directory where the data is

%workingdir = '/Users/andrew/Software/Matlab/osl_test_data/face_data_sub1_osl2_new';

cmd = ['mkdir ' workingdir]; if ~exist(workingdir, 'dir'), unix(cmd); end % make dir to put the results in

clear spm_files_continuous spm_files_epoched;

% set up a list of SPM MEEG object file names (we only have one here)
spm_files_continuous{1}=[workingdir '/fdspm_meg19.mat'];
spm_files_epoched{1}=[workingdir '/Sedspm_meg1.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this section we will do a wholebrain beamformer, followed by a trial-wise
% GLM that will correspond to a comparison of the ERFs for the different
% conditions.

%%%%%%%%%%%%%%%%%%%
%% SETUP THE OAT:
oat=[];
%oat.source_recon.D_continuous=spm_files_continuous;
oat.source_recon.D_epoched=spm_files_epoched;
oat.source_recon.conditions={'Motorbike','Neutral face','Happy face','Fearful face'};
oat.source_recon.freq_range=[1 40]; % frequency range in Hz
oat.source_recon.time_range=[-0.2 0.4]; % time range in secs
oat.source_recon.method='beamform'; % mne_eye, mne_datacov mne_diag_datacov, mne_other_hmm_states_datacov
%oat.source_recon.method='beamform_bilateral'; % mne_eye, mne_datacov mne_diag_datacov
oat.source_recon.normalise_method='mean_eig';

oat.source_recon.gridstep=8; % in mm, using a lower resolution here than you would normally, for computational speed
oat.source_recon.dirname=[spm_files_continuous{1} '_erf_wideband_' oat.source_recon.method]; % directory the oat and results will be stored in

%oat.source_recon.forward_meg='MEG Local Spheres';
oat.source_recon.forward_meg='Single Shell';
oat.source_recon.modalities{1}={'MEGPLANAR', 'MEGMAG'};
oat.source_recon.report.do_source_variance_maps=0;

% design_matrix_summary is a parsimonious description of the design matrix.
% It contains values design_matrix_summary{reg,cond}, where reg is a regressor no. and cond
% is a condition no. This will be used (by expanding the conditions over
% trials) to create the (num_regressors x num_trials) design matrix:
design_matrix_summary={};
design_matrix_summary{1}=[1 0 0 0];design_matrix_summary{2}=[0 1 0 0];design_matrix_summary{3}=[0 0 1 0];design_matrix_summary{4}=[0 0 0 1];
oat.first_level.design_matrix_summary=design_matrix_summary;

% contrasts to be calculated:
oat.first_level.contrast={};
oat.first_level.contrast{1}=[3 0 0 0]'; % motorbikes
oat.first_level.contrast{2}=[0 1 1 1]'; % faces
oat.first_level.contrast{3}=[-3 1 1 1]'; % faces-motorbikes
oat.first_level.contrast_name={};
oat.first_level.contrast_name{1}='motorbikes';
oat.first_level.contrast_name{2}='faces';
oat.first_level.contrast_name{3}='faces-motorbikes';
oat.first_level.report.first_level_cons_to_do=[2 1 3];
oat.first_level.time_range=[-0.1 0.3];
oat.first_level.post_tf_downsample_factor=1;
oat.first_level.name=['wholebrain_first_level'];

oat.first_level.parcellation.do=1;
if oat.first_level.parcellation.do
    tilde='/Users/woolrich/';
    addpath([tilde 'Dropbox/vols_scripts/MEG-ROI-nets']);

    parc_file=[tilde '/Dropbox/vols_scripts/parcellations/fmri_d100_parcellation_with_PCC_reduced_2mm'];
    parcellationfile = [parc_file '_ss5mm_ds8mm'];
    
    %parcellationfile = [tilde '/homedir/vols_data/hmm_investigations/parcellations/aal2mni_cortical_4d_8mm'];

    oat.first_level.parcellation.parcellation=parcellationfile;
    oat.first_level.parcellation.orthogonalisation = 'symmetric';
    oat.first_level.parcellation.method            = 'spatialBasis';
    oat.first_level.parcellation.normalise_voxeldata = 0;
    oat.first_level.name=[oat.first_level.name '_parc' num2str(oat.first_level.parcellation.do)];
end

oat = osl_check_oat(oat);

%%%%%%%%%%%%%%%%%%%
%% RUN THE OAT:

oat.to_do=[0 1 0 0];
oat = osl_run_oat(oat);

% report = oat_source_recon_report(oat);
% report = oat_first_level_stats_report(oat,oat.first_level.results_fnames{1});disp(report.html_fname);

%%%%%%%%%%%%%%%%%%%
%% OUTPUT SUBJECT'S NIFTII FILES
% Having run the GLM on our source space data, we would like to inspect the
% results for our single subject.
% We can do this by saving the contrast of parameter estimates (COPEs) and
% t-statistics for each of our contrasts to NIFTI images.

S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.first_level_contrasts=[3]; % list of first level contrasts to output
%S2.resamp_gridstep=oat.source_recon.gridstep;
S2.resamp_gridstep=oat.source_recon.gridstep;
[statsdir,times,count]=oat_save_nii_stats(S2);

%%%%%%%%%%%%%%%%%%%
%% VIEW NIFTII RESULTS IN FSLVIEW
% We can now view the nifti images containing our GLM results in FSL, here
% we are running fslview from the matlab command line, but you do not need
% to - you can run it from the UNIX command line instead.

mni_brain=[osldir '/std_masks/MNI152_T1_' num2str(S2.resamp_gridstep) 'mm_brain'];

% INSPECT THE RESULTS OF A CONTRAST IN FSLVIEW
contrast=3;
runcmd(['fslview ' mni_brain ' ' [statsdir '/cope' num2str(contrast) '_' num2str(S2.resamp_gridstep) 'mm'] ' ' [statsdir '/tstat' num2str(contrast) '_' num2str(S2.resamp_gridstep) 'mm']  ' ' [statsdir '/tstat' num2str(contrast) '_mip_' num2str(S2.resamp_gridstep) 'mm'] ' &']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Investigating location of interest using an MNI coordinate
% In this section we will interrogate the wholebrain OAT (run above) using
% a specified MNI coordinate.

mni_coord=[27,-64,-18]; % set this to the MNI coord of interest

S2=[];
S2.vox_coord=mni_coord;
S2.stats=oat.first_level.results_fnames{1};
S2.oat=oat;
S2.first_level_cons_to_do=oat.first_level.report.first_level_cons_to_do; % plots all of these contrasts

[vox_ind_used] = oat_plot_vox_stats(S2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Investigating regions of interest using an MNI mask
% In this section we will interrogate the wholebrain OAT (run above) using
% an ROI mask

% Spatially average the results over an ROI
S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.mask_fname=[OSLDIR '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex'];
[stats,times,mni_coords_used]=oat_output_roi_stats(S2);

% plot
S2=[];
S2.stats=stats;
S2.oat=oat;
S2.first_level_cons_to_do=oat.first_level.report.first_level_cons_to_do; % plots all of these contrasts
[vox_ind_used] = oat_plot_vox_stats(S2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROI Time-frequency analysis
% This section will re-rerun the first level analysis using a ROI in
% the temporal occiptial fusiform cortex, and do a time-frequency
% trial-wise GLM first-level analysis

% We are going to use the wholebrain OAT (which was run above), to make use of the settings
% and source_recon results already in there

% Give the first level analysis a new name to avoid copying over previous
% first-level analyses:
oat.first_level.name='roi_tf_first_level';
oat.first_level.mask_fname=[OSLDIR '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex'];

%oat.source_recon.mask_fname=[OSLDIR '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex'];
%oat.source_recon=rmfield(oat.source_recon,'mask_fname');
%oat.first_level=rmfield(oat.first_level,'mask_fname');
%oat.source_recon.mni_coords=[26   -54   -16];

oat.first_level.tf_freq_range=[4 48]; % frequency range in Hz
oat.first_level.time_range=[-0.1 0.3];
oat.first_level.tf_num_freqs=12;
oat.first_level.tf_method='hilbert';
oat.first_level.tf_hilbert_freq_res=8;
oat.first_level.post_tf_downsample_factor=2;

oat = osl_check_oat(oat);

oat.to_do=[1 1 0 0];
oat = osl_run_oat(oat);

% Spatially average the results over an ROI
S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.mask_fname=[OSLDIR '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex'];
[stats,times,mni_coords_used]=oat_output_roi_stats(S2);

%%%%%%%%
%% Plot the time-frequency image for faces-motorbikes
S2=[];
S2.stats=stats;
S2.oat=oat;
S2.first_level_cons_to_do=3; % plots all of these contrasts
[vox_ind_used] = oat_plot_vox_stats(S2);

%%%%%%%%
%% Plot the time course of a single freq bin
freqbin=nearest(stats.frequencies,8); % find bin for 8Hz

S2=[];
S2.stats=stats;
S2.oat=oat;
S2.freq_inds=freqbin;
S2.first_level_cons_to_do=3; % plots all of these contrasts
[vox_ind_used] = oat_plot_vox_stats(S2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Broadband wholebrain TF amplitude analysis
% This section will re-rerun the first level analysis over the whole brain,
% and do a time-frequency trial-wise GLM first-level analysis across a
% wide 4-40Hz band

% LOAD PREVIOUSLY RUN OAT
% We are going to use the wholebrain OAT (which was run above for the ERF
% analysis), to make use of the settings already in there
oat.source_recon.dirname=[spm_files_continuous{1} '_erf_wideband_' oat.source_recon.method]; % directory the oat and results will be stored in
oat.first_level.name=['wholebrain_first_level'];
oat=osl_load_oat(oat);
try, oat.first_level=rmfield(oat.first_level,'freq_average'); catch, end;
res=oat_load_results(oat,oat.source_recon.results_fnames{1});

% RUN THE OAT:
oat.to_do=[1 1 0 0];
oat.source_recon.gridstep=14; % in mm, using a lower resolution here than you would normally, for computational speed
oat.source_recon.freq_range=[4 40]; % broadband

oat.first_level.tf_num_freqs=1;
oat.first_level.tf_method='hilbert';
oat.first_level.tf_freq_range=oat.source_recon.freq_range;
oat.first_level.post_tf_downsample_factor=2;
oat.first_level.name=['wholebrain_tf_first_level'];

oat = osl_run_oat(oat);

%% OUTPUT SUBJECT'S NIFTII FILES
% output and view niis
S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.first_level_contrasts=[3]; % list of first level contrasts to output
S2.freq_bin=1;
S2.stats_dir=[oat.source_recon.dirname '/' oat.first_level.name '_f' num2str(S2.freq_bin) '_stats_dir']; % directory to put niis in
S2.resamp_gridstep=oat.source_recon.gridstep;
[statsdir,times]=oat_save_nii_stats(S2);

% view results using fslview
mni_brain=[osldir '/std_masks/MNI152_T1_' num2str(S2.resamp_gridstep) 'mm_brain'];
contrast_num=4;
runcmd(['fslview ' mni_brain ' ' statsdir '/tstat' num2str(contrast_num) '_' num2str(S2.resamp_gridstep) 'mm &']);


%% Compute a laterality contrast for the faces condition

oat.subject_level.compute_laterality = [2];
%oat.subject_level.merge_contraipsi = [1 2]; % this will flip 1 and not 2.

oat.to_do = [ 0 0 1 0];

oat = osl_check_oat( oat );
oat = osl_run_oat( oat );

res=oat_load_results(oat,oat.subject_level.results_fnames{1});

S2=[];
S2.oat=oat;
S2.stats_fname=oat.subject_level.results_fnames{1};
S2.first_level_contrasts=[1:4]; % list of first level contrasts to output
S2.freq_bin=1;
S2.stats_dir=[oat.source_recon.dirname '/' oat.first_level.name '_f' num2str(S2.freq_bin) '_stats_dir']; % directory to put niis in
S2.resamp_gridstep=oat.source_recon.gridstep;
[statsdir,times]=oat_save_nii_stats(S2);

fslview([ workingdir '/fdspm_meg19.mat_erf_wideband_beamform.oat/wholebrain_first_level_f1_stats_dir/tstat4_6mm.nii.gz']);
