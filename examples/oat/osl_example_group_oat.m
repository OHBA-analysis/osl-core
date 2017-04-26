%% OAT 4 - Group Analysis
%
% This practical will work with a group-level OAT analysis run across 40
% participants in the faces-motorbikes task we looked at in single subject
% practical sessions.
% 
% Work your way through the script cell by cell using the supplied dataset.
% As well as following the instructions below, make sure that you read all
% of the comments (indicated by %), as these explain what each step is
% doing. Note that you can run a cell (marked by %%) using the ?Cell? drop
% down menu on the Matlab GUI.    
%

%% LOAD PRE-RUN GROUP OAT
%
% A group analysis needs the source_recon, first_level and subject_level
% analyses to be completed for each participant and the group_level to be
% run across all the lower level results.
%
% This takes quite a while to run for a 40 subject dataset, so to save time
% for this tutorial we have provided an completed OAT analysis which we can load in and
% explore.
%
% This section loads the  OAT analysis for which the first 3 stages (source recon,
% first-level GLM, subject-level averaging) have already been run.
% Note that the 1st level contrasts that have been run are:
% S2.contrast{1}=[3 0 0 0]'; % motorbikes
% S2.contrast{2}=[0 1 1 1]'; % faces
% S2.contrast{3}=[-3 1 1 1]'; % faces-motorbikes
%

% the directory containing the group OAT
workingdir = fullfile(osldir,'example_data','faces_group');

% Load in the previous OAT
oat = osl_load_oat([osldir 'example_data/faces_group/sourcespace'], 'wholebrain_first_level','sub_level','group_level'); 

%% EXAMINE THE LOWER LEVEL ANALYSES
%
% Take a look at the options for the source-recon and first-level analyses.
% These are very similar to the parameters we chose for the single subject
% ERF beamformer tutorial
%
% Take a moment to familiarise yourself with the options. In particular,
% look at the first level GLM settings and which contrasts have been run.

disp(oat.source_recon);
disp(oat.first_level);

%% SUBJECT LEVEL OPTIONS
%
% The subject level allows for fixed-effects combination of different
% datasets recorded from a single participant. To clarify some terminology:
%
% * session - a session at the |source_recon| or |first_level| is a single
% dataset acquired at a single point in time
% * subject - a single participant who has contributed one or more sessions
% to the analysis.
%
% If we have a design in which one subject has contributed several
% sessions, we can combine the COPEs from each session to get a subject
% level estimate. This subject level estimate is then passed into the
% group_level analysis.
%
% There are some extra options for computing laterality contrasts at the
% subject level. This is useful if you was to compare a response or effect
% in the left or right hemisphere.
%
% The key parameter in the subject level is
% |oat.subject_level.session_index_list|. This is a cell array with a cell
% containing the indices of all the first_level sessions contributed by
% that subject.
%
% In our case, each subject contributed a single session so this indexing
% is a simple one-to-one matching between the |first_level| and
% |subject_level|

disp(oat.subject_level);
disp(oat.subject_level.session_index_list);

%% SETUP GROUP-LEVEL OPTIONS
%
% This section defines the parameters for the group_level in the OAT analysis.
% Note that we do not actually run this analyses here, but these options
% define the analysis as it was run
%
% The group design matrix is defined as a single vector of ones, this
% calculates a group average across all participants. Finally the first and
% group level contrasts to run for the report are set up. IN this case we are
% running first level contrast 3 (faces > motorbikes) and group level contrast
% 1 (grand mean).

oat.group_level.name='group_level';
oat.group_level.subjects_to_do=[1:42];

% Spatial and temporal averaging options
oat.group_level.time_range=[-0.1 0.3];
oat.group_level.space_average=0;
oat.group_level.time_average=0;
oat.group_level.time_smooth_std=0; % secs
oat.group_level.use_tstat=0;

% Spatial and temporal smoothing options
oat.group_level.spatial_smooth_fwhm=0; % mm
oat.group_level.group_varcope_time_smooth_std=0;
oat.group_level.group_varcope_spatial_smooth_fwhm=100; % smooths the variance of the group copes. It is recommended to do this.

% Set up design matrix and contrasts
oat.group_level.group_design_matrix='';
oat.group_level=rmfield(oat.group_level,'group_design_matrix');
oat.group_level.group_design_matrix=ones(1,length(oat.group_level.subjects_to_do));
oat.group_level.group_contrast=[];   
oat.group_level.group_contrast{1}=[1]';  
oat.group_level.group_contrast_name={};
oat.group_level.group_contrast_name{1}='mean';

% Define which contrasts to perform for the report
oat.group_level.first_level_contrasts_to_do=[1:3]; % list of first level contrasts to run the group analysis on
oat.group_level.report.first_level_cons_to_do=[2,1,3]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.group_level_cons_to_do=[1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.show_lower_level_copes=0;
oat.group_level.report.show_lower_level_cope_maps=0;

%% CHECK OAT
%
% The OAT structure that we have created should be passed to osl_run_oat to
% perform the pipeline as defined. However the OAT structure should also
% contain oat.to_do which is a list of binary values indicating which stages of
% the OAT to run.
%
% As the OAT analysis we loaded at the start of this practical has already
% completed we do not need to run the oat here.

oat.to_do=[0 0 0 1]; % run group-level stage only

oat = osl_check_oat(oat);


%% CREATE A STATS REPORT
%
% We will  create a stats report (without rerunning OAT), where we look 
% for the maximum stat for a different first level contrast 
% (specified by S.first_level_con), within the time range S.time_range:

oat.results.plotsdir =fullfile(osldir,'example_data','faces_group','sourcespace.oat','plots');

oat.group_level.report.time_range=[0.07 0.13];
oat.group_level.report.first_level_cons_to_do=[2,1,3]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.group_level_cons_to_do=[1]; % we want to isolate the group average response
oat.group_level.report.show_lower_level_cope_maps=0; % 
report = oat_group_level_stats_report(oat,oat.group_level.results_fnames);

%% VIEW GROUP LEVEL STATS REPORT
%
% Click the link to open the stats report after the previous cell has
% completed.
%
% The report shows a range of plots which are useful for checking the
% analysis has run properly and summarising the key results.
%
% On first page we see:
%
% *Group-level design matrix* - This is the GLM design matrix from our
% analysis
%
% *Average lower level STDCOPEs* - These show the variance of the COPE
% estimate averaged over time for each subject for each contrast. If any subject has an
% extremely high or low variance it might indicate that their are an
% outlier and we should double check their lower level results. Perhaps
% they have very few good trials in one condition or their coregistration
% has gone wrong.
%
% *Stats for GC1* - These are the group average COPEs and t-stats for the lower level
% contrasts.
%
% *Lower level COPES* - These show the COPES for each individual subject
% for each lower-level contrast. The thick line indicates the group mean. 
% Again it is good to check these to identify any outliers.
%
% * Lower-level STDCOPES* - These are the lower-level COPEs (as we saw
% above) from the maximal time-point in the faces contrast.
%
% There are no apparent outliers in our diagnostic plots so we can click
% the links at the top of the page to explore the results in more detail.
% Click the link to open the results for the faces contrast.
%
% Here we see the group level COPE and t-stats across the brain for the
% faces contrast. The time-point is chosen by the maximum response in the
% faces condition. There is a clear occipital response around 100ms after
% stimulus onset. Explore the results from the other contrasts before
% moving on to the next cell.
%
% <<osl_example_group_oat_maxt_smap_c2_gc1.png>>


%% OUTPUT GROUP'S NIFTII FILES
%
% As with the single subject analyses we can output nifti files containig
% the group-level results across the whole brain and the whole experimental
% epoch.
%
% This section saves out a nifti volume containing the results from first_level
% contrast 3 (faces > motorbikes) across the whole group

S2=[];
S2.oat=oat;
S2.stats_fname=oat.group_level.results_fnames;
S2.first_level_contrasts=[3,2,1]; % list of first level contrasts to output
S2.group_level_contrasts=[1]; % we want to isolate the group average response
S2.resamp_gridstep=8;
[statsdir,times]=oat_save_nii_stats(S2);

%% VIEW NIFTII RESULTS IN FSLVIEW
%
% Run the code below to open FSLVIEW with the group-level results for the
% faces-motorbikes contrast along with a structural image and a mask of the
% right hemisphere fusiform cortex.
%
% Once FSLVIEW is open, explore the volume to find the green Fusiform Mask.
% This is brain region which is closely associated with the processing of
% visual information and is particularly specialised for processing faces.
%
% Make sure the |tstat3_gc1_8mm| volume is highlighted by clicking on it in
% the white box at the bottom of the FSLVIEW window.
%
% Change the Volume in the bottom left to 52, this corresponds to 100ms
% after stimulus onset. We can see a very strong response in primary visual
% cortex in the medial part of the occipital lobes.
%
% Now change the Volume to 63, this corresponds to ~150ms after stimulus
% onset. The primary visual response has finised and now the region
% highlighted by the Green fusiform mask has a much stronger response.
%

% display the times associated with each volume in the nifti output.
disp(times)

% Find the paths to the relevant results
tstat = fullfile( oat.source_recon.dirname,'wholebrain_first_level_sub_level_group_level_dir','tstat3_gc1_8mm.nii.gz' );
fusiform = fullfile(osldir,'example_data','faces_groupdata2','structurals','Right_Temporal_Occipital_Fusiform_Cortex_8mm.nii.gz');

% Display the results in FSLVIEW
fslview( {fusiform; tstat}, [0 5;10 15], {'Green';'Red-Yellow'} );

%% INVESTIGATING LOCATIONS OF INTEREST USING AN MNI COORDINATE
%
% We can isolate the results from a specific voxel using
% |oat_plot_vox_stats|. Here we will look at the group results for all
% three lower level contrasts in the visual cortex.
%
% Try changing the co-ordinate to 32,-64,-18. What differences do you
% notice

mni_coord=[4,-82,-8];

S2=[];
S2.oat=oat;
S2.stats=oat.group_level.results_fnames;
S2.vox_coord=mni_coord;
S2.first_level_cons_to_do = [2,1,3];
S2.group_level_cons_to_do = [1];
oat_plot_vox_stats(S2);

%% INVESTIGATING REGIONS OF INTEREST USING AN MNI MASK
%
% In this section we will interrogate the wholebrain OAT (run above) using 
% an ROI mask. This will provide the results averaged across all voxels in
% the defined mask.
%
% We will use the right hemisphere Fusiform Cortex mask which we used in
% FSLVIEW in a previous section.

% load OAT analysis for which the first 4 stages have already been run
oatdir=[osldir 'example_data/faces_groupdata2/sourcespace.oat'];
oat = osl_load_oat(oatdir,'wholebrain_first_level','sub_level','group_level'); 

% Spatially average the results over an ROI
S2=[];
S2.oat=oat;
S2.stats_fname=oat.group_level.results_fnames;
S2.mask_fname=[osldir '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex_8mm.nii.gz'];
[stats,times,mni_coords_used]=oat_output_roi_stats(S2);

% plot
S2=[];
S2.stats=stats;
S2.oat=oat;
S2.first_level_cons_to_do = [2,1,3];
S2.group_level_cons_to_do = [1];
[vox_ind_used] = oat_plot_vox_stats(S2);
