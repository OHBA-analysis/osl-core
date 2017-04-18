%% BEAMFORMER ANALYSIS WITH OAT
% This practical will work with a single subject's data from an emotional
% faces experiment (Elekta Neuromag data). 
%
% Work your way through the script cell by cell using the supplied dataset.
% As well as following the instructions below, make sure that you read all
% of the comments and understand each step as you go.

%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

% setenv('OSLDIR','/Users/andrew/Software/Matlab/osl2.1/')
% addpath(genpath(getenv('OSLDIR')));
% osl_startup(osldir);

%% INITIALISE GLOBAL SETTINGS FOR THIS ANALYSIS
% This cell sets the directory that OAT will work in. Change the workingdir variable to correspond to the correct directory on your computer before running the cell.

% directory where the data is:
datadir = fullfile(osldir,'example_data','faces_singlesubject','spm_files');
%datadir = fullfile(osldir,'example_data','faces_subject1_data_osl2');

% directory to put the analysis in
workingdir = fullfile(osldir,'example_data','faces_subject1_data_osl2');

%% SET UP THE SUBJECTS FOR THE ANALYSIS
%
% Specify a list of the fif files, structural files (not applicable for this practical) and SPM files (which will be created). It is important to make sure that the order of these lists is consistent across sessions. Note that here we only have 1 subject, but more generally there would be more than one. For example:
%
% fif_files{1}=[testdir '/fifs/sub1_face_sss.fif'];
% fif_files{2}=[testdir '/fifs/sub2_face_sss.fif'];
% etc...
% spm_files{1} = [workingdir '/sub1_face_sss.mat'];
% spm_files{2} = [workingdir '/sub2_face_sss.mat'];
% etc...

spm_files_continuous{1}=[datadir '/Aface_meg1.mat'];
spm_files_epoched{1}=[datadir '/eAface_meg1.mat'];


%% SETUP THE OAT:
%
% The oat.source_recon options define the parameters for the data
% filtering, windowing and beamforming. We define the D files, conditions
% and time-frequency options in the same way as the sensorspace OAT. THe 
% In this section we will do a wholebrain beamformer, followed by a trial-wise
% GLM that will correspond to a comparison of the ERFs for the different
% conditions. The source-space projection is defined by a new set of
% options. 
%
% The critical options are:
% * method - Which algorithm to use, We will use 'beamform'
% * normalise_method - How to normalise the magnetometers and gradiometers
% * gridstep - This sets the resolution of the grid through the brain. We
% are using 8mm which is lower than usual but faster to compute.
% * forward_meg - This specifies the forward model used.
% * modalidities - Defines which types of sensor to use.
% * do_source_variance_maps - If set to 1, this outputs an optional sanity check

% These options are the same as the sensorspace OAT and define input-data,
% conditions, filtering and windowing.
oat=[];
oat.source_recon.D_continuous=spm_files_continuous;
oat.source_recon.D_epoched=spm_files_epoched;
oat.source_recon.conditions={'Motorbike','Neutral face','Happy face','Fearful face'};
oat.source_recon.freq_range=[3 20]; % frequency range in Hz
oat.source_recon.time_range=[-0.2 0.4]; % time range in secs

% These options specify the source reconstruction that will take place.
oat.source_recon.method='beamform';
oat.source_recon.normalise_method='mean_eig';
oat.source_recon.gridstep=8; % in mm
oat.source_recon.forward_meg='Single Shell';
oat.source_recon.modalities{1}={'MEGPLANAR', 'MEGMAG'};
oat.source_recon.report.do_source_variance_maps=1;

oat.source_recon.dirname=[workingdir '/beamformer_erf']; % directory the oat and results will be stored in


%% SPECIFIY FIRST LEVEL OPTIONS
%
% These options are the same as the sensorspace ERF tutorial.
%
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

oat = osl_check_oat(oat);

%% RUN THE OAT:
%
% We only need to run the source_recon and first_level stages in this
% tutorial, the subject_level and group_level options can be set to 0.
%
% The OAT will produce and close a number of figures as it processes, we
% will discuss what they mean once it has finished (this takes a couple of
% minutes)

oat.to_do=[1 1 0 0];
oat = osl_run_oat(oat);


%% VIEW OAT  REPORT
%
% Once finished, the OAT will print a link to a html report. Click to open
% it in the matlab html browser.
%
% *Source Recon*
%
% Firstly, click to open the "Session_1_reoprt". This contains the plots
% relevent to the source_recon. In our case this is the sensor
% normalisation. Note that the eigenspectrum of the sensordata 
% ('Pre-normalised log eigenspectrum') drops off sharply at
% around 64, this indicates the rank of the sensor data which limited by
% maxfilter de-noising.
%
% The Pre-normalised variances shows the sensor-by-sensor variance across
% time. Channels 200-300 are the Magnetometers and have much higher
% variance than the Gradiometers. We need to normalise the data to remove
% this disparity or the information in the Magnetometers will dominate the
% reconstruction.
% 
% In the normalised_eigs plot below, we can see that the 'mean_eig' sensor
% normalisation has both removed the shelf in the eigenspectrum and brought
% the sensor variances in line with each other.
%
% *First Level*
%
% Click back to the first page and on to the 'First level (epoched) report'
% and then the 'session1 report'. This contains a results summary from the
% voxel-wise GLM
%
% The first figure is the design matrix from the GLM, this is the same as
% the sensorspace OAT.
%
% Next we see a set of time-series from the voxel with the peak response to
% the Faces contrast. The COPE is on the left and the t-stat on the right.
% Note the large peak around 150ms after stimulus onset.
%
% Finally the source maps for each contrast are shown for the peak time-point in the Faces
% contrast. This time-point is around 136ms and shows peaks in lateral and
% medial occipital cortex. These are right hemisphere lateralised for the
% Faces and Faces-Motorbikes condition.
%
% <<osl_example_beamfomer_oat_05.png>>

%% OUTPUT SUBJECT'S NIFTII FILES
%
% The html report gives a brief summary of the results but we would like to
% go into more detail.
%
% We can do this by saving the contrast of parameter estimates (COPEs) and
% t-statistics for each of our contrasts to NIFTI images.

S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.first_level_contrasts=[3,1,2]; 
S2.resamp_gridstep=oat.source_recon.gridstep;
[statsdir,times,count]=oat_save_nii_stats(S2);

%% OPEN NIFTII RESULTS IN FSLVIEW
% We can now view the nifti images containing our GLM results in FSL, here
% we are running fslview from the matlab command line, but you do not need
% to - you can run it from the UNIX command line instead.

mni_brain=[osldir '/std_masks/MNI152_T1_' num2str(S2.resamp_gridstep) 'mm_brain.nii.gz'];

% Inspect the results of an OAT contrast in FSLVIEW
contrast=3;
runcmd(['fslview ' mni_brain ' ' [statsdir '/tstat' num2str(contrast) '_' num2str(S2.resamp_gridstep) 'mm'] ' &']);

%% VIEW RESULTS IN FLSVIEW
%
% The previous command should have opened FSLVIEW with the t-stats for the
% faces contrast. We need to tweat the viewing settings to see the results
% well.
%
% Firstly, make sure the 'tstat2_8mm' image in selected by clicking on it
% once in the bottom window. Next set the 'Min' and 'Max' to 15 and 25
% respectively in the two boxes in the middel at the top of the screen.
%
% We are currently looking at the first time-point in the source_recon
% window which corresponds to -100ms, in which not much is happening.
%
% Change the 'Volume' in the bottom left to 49, this corresponds to around
% 100ms after stimulus onset. You should see a response in early visual
% cortex at the back of the brain. You can cross-check the times referred
% to in the Volumes with the 'times' variable returned by
% 'oat_save_nii_stats' above (Note that FSLVIEW indexes Volumes from 0 and
% Matlab indexes times from 1)
%
% Now change the 'Volume' to 59, corresponding to around 140ms after
% stimulus onset. The peak of activation jumps to the right hemiphere
% visual cortex.
%
% You can futher investigate the results by clicking on the
% 'Tools->Timeseries' option in the top menu. This will bring up and new
% window showing the t-value across time for the voxel under the green
% curser. You can navigate across space by clicking on a part of the brain
% and time by clicking at a time-point on the timeseries.
%
% Try changing to code in the previous cell to bring up the tstats for
% contrast 3 'Face-Motorbikes'. Following the same visualisation
% instructions you should be able to see that the early time-window (Volume
% 49) does not contain any large responses whereas the later response
% (Volume 59) still occurs. This indicates that there are no large
% differences between Faces and Motorbikes in the the early response,
% whereas the later response does differ.

%% INVESTIGATING LOCATION OF INTEREST USING AN MNI COORDINATE
%
% We may want to see the results across all contrasts for a single ROI to
% gain another perspective on our results.
%
% Here we will interrogate the wholebrain OAT (run above) using
% a specified MNI coordinate. This will bring up a the COPE and tstat
% estimates across time for a voxel in visual cortex. Note the prominant
% response around 100ms
%
% Try changing the code to run at 32,-64,-18
%
% This corresponds to a point in Right Hemisphere Fusiform Cortex. Note
% that the 100ms response does not appear here, rather the later Face
% specific response is more dominant.

mni_coord=[4,-82,-8]; % Visual Cortex Voxel

S2=[];
S2.vox_coord=mni_coord;
S2.stats=oat.first_level.results_fnames{1};
S2.oat=oat;
S2.first_level_cons_to_do=oat.first_level.report.first_level_cons_to_do; % plots all of these contrasts

[vox_ind_used] = oat_plot_vox_stats(S2);

%% INVESTIGATING REGIONS OF INTEREST USING AN MNI MASK
%
% We can also interrogate the wholebrain OAT (run above) using
% an ROI mask rather than a single voxel.
%
% In this case we will look at the Right Hemisphere Fusiform Cortex.

% Apply a mask and spatially average the results over an ROI
S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.mask_fname=[osldir '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex_8mm.nii.gz'];
[stats,times,mni_coords_used]=oat_output_roi_stats(S2);

% Plot the COPEs and t-stats within the ROI
S2=[];
S2.stats=stats;
S2.oat=oat;
S2.first_level_cons_to_do=oat.first_level.report.first_level_cons_to_do; % plots all of these contrasts
[vox_ind_used] = oat_plot_vox_stats(S2);


%% ROI TIME-FREQUENCY ANALYSIS
%
% We can extend the OAT to run in time-frequency space. Here we will run a
% time-frequency OAT restricted to the RH Fusiform ROI used above.
%
% Most of the OAT settings do not need to be changed. We do need to define
% the mask in the source recon stage and add the time-frequency
% decomposition parameters in the first level.


% Give the first level analysis a new name to avoid copying over previous
% first-level analyses:
oat.first_level.name='roi_tf_first_level';

% Re-use the mask we used above
oat.source_recon.mask_fname=[osldir '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex'];

% Add first level source recon options
oat.first_level.tf_freq_range=[4 48]; % frequency range in Hz
oat.first_level.time_range=[-0.1 0.3];
oat.first_level.tf_num_freqs=12;
oat.first_level.tf_method='hilbert';
oat.first_level.tf_hilbert_freq_res=8;
oat.first_level.post_tf_downsample_factor=2;

% Check and run the source_recon and first_level OAT stges
oat = osl_check_oat(oat);

oat.to_do=[1 1 0 0];
oat = osl_run_oat(oat);

% Spatially average the results over an ROI
S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.mask_fname=[osldir '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex_8mm.nii.gz'];
[stats,times,mni_coords_used]=oat_output_roi_stats(S2);


%% SOURCE ROI TIME-FREQUENCY PLOT
%
% As with the voxel example able, we can interrogate the OAT to extract the
% time-frequency plots from specific contrasts.

S2=[];
S2.stats=stats;
S2.oat=oat;
S2.first_level_cons_to_do=3; % plots all of these contrasts
[vox_ind_used] = oat_plot_vox_stats(S2);


%% SOURCE ROI SINGLE-FREQUENCY POWER PLOT
%
% We may also isolate the results from a single frequency by defining the
% freq_inds parameter

freqbin=nearest(stats.frequencies,8); % find bin for 8Hz

S2=[];
S2.stats=stats;
S2.oat=oat;
S2.freq_inds=freqbin;
S2.first_level_cons_to_do=3; % plots all of these contrasts
[vox_ind_used] = oat_plot_vox_stats(S2);

%% BROADBAND WHOLEBRAIN TF ANALYSIS
%
% Next we will extend the OAT to run a time-frequecy analsis across the
% whole brain.
%
% We start by loading in our first OAT analysis, as with the second analysis
% only a few parameters need to be changed. Firstly we reduce the gridstep
% in the source recon (just to ensure it runs within the session) and
% secondly we add the first_level time-frequency decomposition options.
%
% These options will do a time-frequency trial-wise GLM first-level analysis across a
% wide 4-40Hz band

% LOAD PREVIOUSLY RUN OAT
% We are going to use the wholebrain OAT (which was run above for the ERF
% analysis), to make use of the settings already in there
oat.source_recon.dirname=[workingdir '/beamformer_erf']; % Make sure this matches the dirname used above
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
%
% Again, we can output the nifti files from this participant to explore the
% results in more detail

% output and view niis
S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.first_level_contrasts=[2]; % list of first level contrasts to output
S2.freq_bin=1;
S2.stats_dir=[oat.source_recon.dirname '/' oat.first_level.name '_f' num2str(S2.freq_bin) '_stats_dir']; % directory to put niis in
S2.resamp_gridstep=oat.source_recon.gridstep;
[statsdir,times]=oat_save_nii_stats(S2);

% view results using fslview
mni_brain=[osldir '/std_masks/MNI152_T1_' num2str(S2.resamp_gridstep) 'mm_brain.nii.gz'];
contrast_num=2;
runcmd(['fslview ' mni_brain ' ' statsdir '/tstat' num2str(contrast_num) '_' num2str(S2.resamp_gridstep) 'mm &']);
