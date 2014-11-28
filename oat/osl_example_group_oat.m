%% PRACTICAL: GROUP ANALYSIS
%
% This practical will work with a single subject's data from an emotional
% faces experiment (Elekta Neuromag data). You can get the data from:
% https://sites.google.com/site/ohbaosl/practicals/practical-data/emotional-face-processing-elekta-neuromag-data
% 
% Work your way through the script cell by cell using the supplied dataset.
% As well as following the instructions below, make sure that you read all
% of the comments (indicated by %), as these explain what each step is
% doing. Note that you can run a cell (marked by %%) using the ?Cell? drop
% down menu on the Matlab GUI.    
%
% get data from: www.fmrib.ox.ac.uk/~woolrich/faces_group_data.tar.gz

%% SETUP THE MATLAB PATHS
%
% Sets the Matlab paths to include OSL. Change these paths so that they 
% correspond to the setup on your computer. You will also need to ensure 
% that fieldtrip and spm are not in your matlab path (as they are included
% within OSL).%% SETUP THE MATLAB PATHS

global OSLDIR;
    
% This cell sets the Matlab paths to include OSL. Change the osldir path so
% that it corresponds to the setup on your computer before running the cell. 
osldir = '/home/mwoolrich/homedir/matlab/osl1.5.0_beta';

addpath(osldir);
osl_startup(osldir);

%% INITIALISE GLOBAL SETTINGS FOR THIS ANALYSIS
%
% This cell sets the directory that OAT will work in. Change the workingdir
% variable to correspond to the correct directory on your computer before
% running the cell.
workingdir = '/home/mwoolrich/Desktop/faces_group_data_osl1.5.0';

%% Load pre-run OAT analysis
%
% Load OAT analysis for which the first 3 stages (source recon,
% first-level GLM, subject-level averaging) have already been run.
% Note that the 1st level contrasts that have been run are:
% S2.contrast{1}=[3 0 0 0]'; % motorbikes
% S2.contrast{2}=[0 1 1 1]'; % faces
% S2.contrast{3}=[-3 1 1 1]'; % faces-motorbikes

oatdir=[workingdir '/beamform.oat'];
oat = osl_load_oat(oatdir,'first_level_none','sub_level','group_level'); 

%% Setup group-level stage
%
% This section defines the parameters for the group_level in the OAT analysis.
% Note that we are not performing any spatial or temporal averaging in this
% analysis.
%
% The group design matrix is defined as a single vector of ones, this
% calculates a group average across all participants. Finally the first and
% group level contrasts to run for the report are set up. IN this case we are
% running first level contrast 3 (faces > motorbikes) and group level contrast
% 1 (grand mean).

oat.group_level.time_range=[-0.1 0.3];
oat.group_level.space_average=0;
oat.group_level.time_average=0;
oat.group_level.time_smooth_std=0; % secs
oat.group_level.spatial_smooth_fwhm=0; % mm
oat.group_level.mask_fname='';oat.group_level=rmfield(oat.group_level,'mask_fname');
oat.group_level.use_tstat=0;
oat.group_level.group_varcope_time_smooth_std=0;
oat.group_level.group_varcope_spatial_smooth_fwhm=100; % smooths the variance of the group copes. It is recommended to do this.
oat.group_level.name='group_level';
oat.group_level.subjects_to_do=[3:14];


% Set up design matrix and contrasts
oat.group_level.group_design_matrix='';oat.group_level=rmfield(oat.group_level,'group_design_matrix');
oat.group_level.group_design_matrix=ones(1,length(oat.group_level.subjects_to_do));
oat.group_level.group_contrast=[];   
oat.group_level.group_contrast{1}=[1]';  
oat.group_level.group_contrast_name={};
oat.group_level.group_contrast_name{1}='mean';


% Define which contrasts to perform for the report
oat.group_level.first_level_contrasts_to_do=[1:3]; % list of first level contrasts to run the group analysis on
oat.group_level.report.first_level_cons_to_do=[3]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.group_level_cons_to_do=[1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.show_lower_level_copes=1;
oat.group_level.report.show_lower_level_cope_maps=0;

%% Check OAT
%
% The OAT structure that we have created should be passed to osl_run_oat to
% perform the pipeline as defined. However the OAT structure should also
% contain oat.to_do which is a list of binary values indicating which stages of
% the OAT to run.
%
% As the OAT analysis we loaded at the start of this practical contained a
% completed source_recon, first_level and subject_level analysis, we only need
% to run the group level here

oat.to_do=[0 0 0 1]; % run group-level stage only

oat = osl_check_oat(oat);

%% RUN OAT

oat = osl_run_oat(oat);

%% Create an alternative stats report
%
% The diagnostic report created above is based on showing the results at the
% voxel and time of the maximum statistic for the first-level contrast
% specified first in the setting
% oat.group_level.diagnostic_first_level_cons_to_do
% (i.e. contrast 1), and at the group level contrast specified first in 
% oat.group_level.diagnostic_group_cons_to_do (i.e. contrast 1)
%
% We will not create a new report (without rerunning OAT), where we look 
% for the maximum stat for a different first level contrast 
% (specified by S.first_level_con), within the time range S.time_range:

oat.group_level.report.time_range=[0.13 0.15];
oat.group_level.report.first_level_cons_to_do=[3]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.show_lower_level_cope_maps=1;
report = oat_group_level_stats_report(oat,oat.group_level.results_fnames);

%% OUTPUT GROUP'S NIFTII FILES
%
% This section saves out a nifti volume containing the results from first_level
% contrast 3 (faces > motorbikes) across the whole group

S2=[];
S2.oat=oat;
S2.stats_fname=oat.group_level.results_fnames;
S2.first_level_contrasts=[3]; % list of first level contrasts to output
S2.resamp_gridstep=8;
[statsdir,times]=oat_save_nii_stats(S2);

%% VIEW NIFTII RESULTS IN FSLVIEW
%
% View the results we just saved out

con=3; % view the 3rd contrast (face vs motorbikes)
runcmd(['fslview ' OSLDIR '/std_masks/MNI152_T1_' num2str(S2.resamp_gridstep) 'mm_brain.nii.gz ' OSLDIR '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex_' num2str(S2.resamp_gridstep) 'mm.nii.gz ' statsdir '/cope' num2str(con) '_gc1_' num2str(S2.resamp_gridstep) 'mm.nii.gz ' statsdir '/tstat' num2str(con) '_gc1_' num2str(S2.resamp_gridstep) 'mm.nii.gz ' statsdir '/tstat' num2str(con) '_gc1_mip_' num2str(S2.resamp_gridstep) 'mm.nii.gz &']);

%% Investigating locations of interest using an MNI coordinate
%
% In this section we will interrogate the wholebrain OAT (run above) using 
% an MNI coordinate of interest

% load OAT analysis for which the first 4 stages have already been run
oatdir=[workingdir '/beamform.oat'];
oat = osl_load_oat(oatdir,'first_level_none','sub_level','group_level'); 

mni_coord=[20,-68,-6]; %in RFFA

S2=[];
S2.oat=oat;
S2.stats=oat.group_level.results_fnames;
S2.vox_coord=mni_coord;
oat_plot_vox_stats(S2);

%% Investigating regions of interest using an MNI mask
%
% In this section we will interrogate the wholebrain OAT (run above) using 
% an ROI mask

% load OAT analysis for which the first 4 stages have already been run
oatdir=[workingdir '/beamform.oat'];
oat = osl_load_oat(oatdir,'first_level_none','sub_level','group_level'); 

% Spatially average the results over an ROI
S2=[];
S2.oat=oat;
S2.stats_fname=oat.group_level.results_fnames;
S2.mask_fname=[OSLDIR '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex'];
[stats,times,mni_coords_used]=oat_output_roi_stats(S2);

% plot
S2=[];
S2.stats=stats;
S2.oat=oat;
[vox_ind_used] = oat_plot_vox_stats(S2);

%% Using ROI at the input to the group level
%
% Here we rerun the group level stage of OAT, but restricted to only those
% voxels within a specified mask (oat.group_level.mask_fname). 
% These voxels are then spatially averaged 
% over (due to the oat.group_level.space_average=1 setting). Note that
% the averaging occurs before the group GLM is fit).

%% Setup the ROI OAT
%
% load OAT analysis for which the first 3 stages have already been run
oatdir=[workingdir '/beamform.oat'];
oat = osl_load_oat(oatdir,'first_level_none','sub_level','group_level'); 

% adjust OAT group level settings
oat.group_level.time_range=[-0.15 0.3];
oat.group_level.space_average=1;
oat.group_level.time_average=0;
oat.group_level.time_smooth_std=0; % secs
%oat.group_level.group_varcope_time_smooth_std=0.2;
%oat.group_level.group_varcope_spatial_smooth_fwhm=100;
oat.group_level.spatial_smooth_fwhm=0; % mm
oat.group_level.mask_fname=[OSLDIR '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex'];
oat.group_level.name='roi_average_group_level';
oat.group_level.use_tstat=0;
oat.group_level.first_level_contrasts_to_do=1:3;
oat.group_level.store_lower_level_copes=1;

oat = osl_check_oat(oat);

% only run group level stage of oat:
oat.to_do=[0 0 0 1];

%% Run the ROI OAT and do the plots
oat = osl_run_oat(oat);

% Spatially average the results over the ROI
S2=[];
S2.oat=oat;
S2.stats_fname=oat.group_level.results_fnames;
S2.mask_fname=oat.group_level.mask_fname;
[gstats,times,mni_coords_used]=oat_output_roi_stats(S2);

% plot
S2=[];
S2.stats=gstats;
S2.oat=oat;
[vox_ind_used] = oat_plot_vox_stats(S2);

% individual subject plots:
baseline_correct=1;
con=3;
tmp=((gstats.lower_level_copes{con}));
figure;plot(times,permute(tmp(1,:,:),[2,3,1]));ho;
plot(times,squeeze(median(tmp(1,:,:),2)),'LineWidth',3);
plot4paper('time (secs)','1st level cope'); 
title(['First level cope ' num2str(con)]);

%% Group stats on a single volume within a time window
%
% time window is specifed by oat.group_level.time_range
% desire to do averaging over time within that window is specified by:
% oat.group_level.time_average=1;

% load in OAT with first 3 stages run:
oatdir=[workingdir '/beamform.oat'];
oat = osl_load_oat(oatdir,'first_level_none','sub_level','group_level'); 

% setup and run group level stage of OAT
oat.group_level.time_range=[0.140 0.150];
oat.group_level.time_average=1;

oat.group_level.first_level_contrasts_to_do=[1:3]; % list of first level contrasts to run the group analysis on
oat.group_level.report.first_level_cons_to_do=[3]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.group_level_cons_to_do=[1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat

oat.group_level.name='time_average_group_level';
oat.group_level.use_tstat=0;
oat.group_level.group_varcope_spatial_smooth_fwhm=100;

%oat.group_level=rmfield(oat.group_level,'diagnostic_cons_to_do');

oat = osl_check_oat(oat);

oat.to_do=[0 0 0 1];

oat = osl_run_oat(oat);

% View the single volume results
S2=[];
S2.oat=oat;
S2.stats_fname=oat.group_level.results_fnames;
S2.first_level_contrasts=[3]; % list of first level contrasts to output

[statsdir,times]=oat_save_nii_stats(S2);

con=3;gcon=1;


runcmd(['fslview ' OSLDIR '/std_masks/MNI152_T1_' num2str(oat.source_recon.gridstep) 'mm_brain.nii.gz '  statsdir '/tstat' num2str(con) '_gc' num2str(gcon) '_' num2str(oat.source_recon.gridstep) 'mm.nii.gz ' statsdir '/cope' num2str(con) '_gc' num2str(gcon) '_' num2str(oat.source_recon.gridstep) 'mm.nii.gz &']);

%% Run 3D permutation stats on a single volume
%
% We will now run 3D permutation stats on the time window average volume

% load in OAT from last cell
oatdir=[workingdir '/beamform.oat'];
oat = osl_load_oat(oatdir,'first_level_none','sub_level','time_average_group_level'); 

S=[];
S.oat=oat;
%S.timepoint=[nearest(gstats.times,0.144)];
S.cluster_stats_thresh=6;
S.cluster_stats_nperms=1000; % normally recommend doing 5000 perms
S.first_level_copes_to_do=[3];
S.group_level_copes_to_do=[1];
S.group_varcope_spatial_smooth_fwhm=S.oat.group_level.group_varcope_spatial_smooth_fwhm;
S.write_cluster_script=0;
S.time_range=[0.140 0.150];
S.time_average=1;
S.fsl_version_4p1=0;
[ gstats ] = oat_cluster_permutation_testing( S );

% View permutation stats
con=S.first_level_copes_to_do(1);
runcmd(['fslview ' OSLDIR '/std_masks/MNI152_T1_' num2str(gstats.gridstep) 'mm_brain.nii.gz ' gstats.dir '/tstat' num2str(con) '_gc1_' num2str(gstats.gridstep) 'mm.nii.gz ' gstats.dir '/clustere_tstat' num2str(con) '_gc1_' num2str(gstats.gridstep) 'mm.nii.gz ' gstats.dir '/clustere_corrp_tstat' num2str(con) '_gc1_' num2str(gstats.gridstep) 'mm.nii.gz &']);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROI Time-Freq analysis
%
% We will now run an ROI power analysis over multiple frequency bands
% The oat.group_level.space_average=1 indicates that we want to do spatial
% averaging over the ROI at the group level (the averaging occurs before
% the group GLM is fit).

% load in an OAT for which the first 3 stages have been run to do an ROI
% analysis using the Right_Temporal_Occipital_Fusiform mask
oatdir=[workingdir '/beamform_roi.oat'];
oat = osl_load_oat(oatdir,'first_level_hilbert','sub_level','group_level'); 

% setup Group level of OAT and run
oat.group_level.time_range=oat.first_level.time_range;
oat.group_level.space_average=1;
oat.group_level.time_average=0;
oat.group_level.time_smooth_std=0; % secs
oat.group_level.spatial_smooth_fwhm=0; % mm
oat.group_level.use_tstat=1;
oat.group_level.name='group_level';
oat.group_level.store_lower_level_copes=1;
oat.group_level.subjects_to_do=[3:14];
oat.group_level.time_range=[0 0.35];

oat=osl_check_oat(oat);

oat.to_do=[0 0 0 1];

oat=osl_run_oat(oat);

% Plot the group results

% Spatially average the results over an ROI
S2=[];
S2.oat=oat;
S2.stats_fname=oat.group_level.results_fnames;
S2.mask_fname=[OSLDIR '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex'];
[stats,times,mni_coords_used]=oat_output_roi_stats(S2);

% plot
S2=[];
S2.stats=stats;
S2.oat=oat;
[vox_ind_used] = oat_plot_vox_stats(S2);

%% Do 2D cluster permutation stats 
%
% Perform sign-flip permutation tests on the data from a single point in 
% the time-frequency OAT. This is done on first level contrast 2 (face
% grand mean). This tells us which points in the time frequency analysis
% are statistically significant with p-values corrected for multiple 
% comparisons and cluster size.

% load GLM result
stats=oat_load_results(oat,oat.group_level.results_fnames);

contrast=2; % first level contrast
cluster_forming_threshold=3.5; % threshold used on t-stats to form clusters
num_perms=1000; % num permutations to do - normally recommend doing 5000

% call osl_clustertf to do permutation testing
% this returns corrp, the corrected (over multiple comparisons) P-Value for
% each cluster in a 2D time-freq map
[corrp tstats]= osl_clustertf(permute(stats.lower_level_copes{contrast},[2 4 3 1]),cluster_forming_threshold,num_perms,26);

figure;
subplot(1,2,1);imagesc(stats.times, stats.frequencies, squeeze(tstats));axis xy;
ylabel('frequency (Hz)'); xlabel('time (s)'); colorbar; title(['T-stats' num2str(contrast)]);
subplot(1,2,2);imagesc(stats.times, stats.frequencies, squeeze(corrp));axis xy;
ylabel('frequency (Hz)'); xlabel('time (s)'); colorbar; title(['Cluster extent P-values for cope' num2str(contrast)]);

%%
















































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optional extras: IGNORE if doing the workshop practical
if(0),

   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Permutation stats in 4D
    
    oat = osl_load_oat(oatdir,'wnorm0','sub_level','group_level');    

    S=[];
    S.oat=oat;
    %S.timepoint=[nearest(gstats.times,0.144)];
    
    S.cluster_stats_thresh=4;
    S.cluster_stats_nperms=100;
    S.first_level_copes_to_do=[3];
    S.group_level_copes_to_do=[1];
    S.group_varcope_spatial_smooth_fwhm=0;%S.oat.group_level.group_varcope_spatial_smooth_fwhm;
    S.write_cluster_script=0;    
    S.time_range=[0.130 0.155];
    S.fsl_version_4p1=0;
    S.matlab_exe_name='/Applications/MATLAB_R2012a.app/bin/matlab';
    S.time_average=0;
    [ gstats ] = oat_cluster_permutation_testing( S );

    %%
    
    if(0),
        % if running using script - you will need this extra stuff
        %cluster_stats= % replace ???? with text outputted by oat_cluster_permutation_testing call above
        gstats=oat_load_results(oat,oat.group_level.results_fnames);
        gstats.clusterstats{1,1}=cluster_stats;
        oat_save_results(oat,gstats);
    end;
    
    % VIEW NIFTII RESULTS IN FSLVIEW
    S2=[];
    S2.oat=oat;
    S2.stats_fname=oat.group_level.results_fnames;
    S2.first_level_contrasts=[3]; % list of first level contrasts to output
    S2.resamp_gridstep=8;
    [statsdir,times,clust_times]=oat_save_nii_stats(S2);

    con=3;
    runcmd(['fslview ' OSLDIR '/std_masks/MNI152_T1_' num2str(S2.resamp_gridstep) 'mm_brain.nii.gz ' statsdir '/tstat' num2str(con) '_gc1_' num2str(S2.resamp_gridstep) 'mm.nii.gz &']);
    runcmd(['fslview ' OSLDIR '/std_masks/MNI152_T1_' num2str(S2.resamp_gridstep) 'mm_brain.nii.gz ' statsdir '/tstat' num2str(con) '_gc1_clust4d_corrp_' num2str(S2.resamp_gridstep) 'mm.nii.gz ' statsdir '/tstat' num2str(con) '_gc1_clust4d_tstats_' num2str(S2.resamp_gridstep) 'mm.nii.gz ' statsdir '/tstat' num2str(con) '_gc1_clust4d_' num2str(S2.resamp_gridstep) 'mm.nii.gz &']);
  %%

  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run 3D permutation stats on a single volume
    % We will now run 3D permutation stats on the time window average volume

    % load in OAT from last cell
    %oat = osl_load_oat(oatdir,'first_level','sub_level','time_average_group_level'); 
    oat = osl_load_oat(oatdir,'first_level_none_mask','sub_level','group_level'); 

    % adjust OAT group level settings
oat.group_level.time_range=[-0.15 0.3];
oat.group_level.space_average=0;
oat.group_level.time_average=0;
oat.group_level.time_smooth_std=0; % secs
%oat.group_level.group_varcope_time_smooth_std=0.2;
%oat.group_level.group_varcope_spatial_smooth_fwhm=100;
oat.group_level.spatial_smooth_fwhm=0; % mm
oat.group_level.mask_fname=[OSLDIR '/std_masks/Right_Temporal_Occipital_Fusiform_Cortex'];
oat.group_level.name='roi_average_group_level';
oat.group_level.use_tstat=0;
oat.group_level.subjects_to_do=[1:14];
oat.group_level.first_level_contrasts_to_do=1:3;
oat.group_level.store_lower_level_copes=1;

oat = osl_check_oat(oat);

% only run group level stage of oat:
oat.to_do=[0 0 0 1];

%% Run the ROI OAT and do the plots
oat = osl_run_oat(oat);

% Spatially average the results over the ROI
S2=[];
S2.oat=oat;
S2.stats_fname=oat.group_level.results_fnames;
S2.mask_fname=oat.group_level.mask_fname;
[gstats,times,mni_coords_used]=oat_output_roi_stats(S2);

% plot
S2=[];
S2.stats=gstats;
S2.oat=oat;
[vox_ind_used] = oat_plot_vox_stats(S2);

% individual subject plots:
baseline_correct=1;
con=3;
tmp=((gstats.lower_level_copes{con}));
figure;plot(times,permute(tmp(1,:,:),[2,3,1]));ho;
plot(times,squeeze(median(tmp(1,:,:),2)),'LineWidth',3);
plot4paper('time (secs)','1st level cope'); 
title(['First level cope ' num2str(con)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run 3D permutation stats on a single volume
% We will now run 3D permutation stats on the time window average volume

% load in OAT from last cell
%oat = osl_load_oat(oatdir,'first_level','sub_level','time_average_group_level'); 
    oat = osl_load_oat(oatdir,'first_level_none_mask','sub_level','roi_average_group_level'); 

S=[];
S.oat=oat;
%S.timepoint=[nearest(gstats.times,0.144)];
S.cluster_stats_thresh=4.5;
S.cluster_stats_nperms=1000; % normally recommend doing 5000 perms
S.first_level_copes_to_do=[3];
S.group_level_copes_to_do=[1];
S.group_varcope_spatial_smooth_fwhm=S.oat.group_level.group_varcope_spatial_smooth_fwhm;
S.write_cluster_script=0;
S.time_range=[0.140 0.150];
S.time_average=1;
S.fsl_version_4p1=0;
[ gstats ] = oat_cluster_permutation_testing( S );

% View permutation stats
con=S.first_level_copes_to_do(1);
runcmd(['fslview ' OSLDIR '/std_masks/MNI152_T1_' num2str(gstats.gridstep) 'mm_brain.nii.gz ' gstats.dir '/tstat' num2str(con) '_gc1_' num2str(gstats.gridstep) 'mm.nii.gz ' gstats.dir '/clustere_tstat' num2str(con) '_gc1_' num2str(gstats.gridstep) 'mm.nii.gz ' gstats.dir '/clustere_corrp_tstat' num2str(con) '_gc1_' num2str(gstats.gridstep) 'mm.nii.gz &']);

end;


