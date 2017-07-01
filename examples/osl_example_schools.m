% Work Experience programme
% This is a 1 hour overview of MEG analysis

% Load the ERF results
oat = osl_load_oat('sensorspace_erf.oat')

%% PLOT ERF
D
D = spm_eeg_load(oat.source_recon.D_epoched{1})
D = spm_eeg_load('/Users/romesh/oxford_postdoc/toolboxes/osl/example_data/faces_singlesubject/spm_files/eAface_meg1.mat')
plot(D.time,mean(D(50,:,:),3)); % Now play with different averaging etc.

%% PLOT COLOUR TIMECOURSES
first_level_results=oat.first_level.results_fnames{1}
first_level_results=oat_load_results(oat,first_level_results); 
first_level_cons_to_do=oat.first_level.report.first_level_cons_to_do; % plots all of these contrasts 
S2=oat.first_level.report; 
S2.stats=first_level_results; 
S2.modality=oat.first_level.report.modality_to_do; 
[vox_ind_max time_ind_max freq_ind_max stats max_stat] = oat_find_max_stats( S2 ); 
S2=[]; 
S2.stats=stats; 
S2.chanlabel=stats.chanlabels{vox_ind_max}; % can change to this setting to plot the time courses at the voxel with the max t-stat 
S2.first_level_cons_to_do=first_level_cons_to_do; % plots all of these contrasts 
[fig_handle fig_name fig_title] = oat_plot_vox_stats(S2); 

%% FT COPE or t-state timecourse and topoplots, ask interactive question?
S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.modality='MEGPLANAR'; % can also set this to 'MEGMAG'
S2.first_level_contrast=[3]; % view faces-motorbikes contrast
S2.view_cope=1; % set to 0 to see the t-stat

% calculate t-stat using contrast of absolute value of parameter estimates
[cfg, dats, fig_handle]=oat_stats_multiplotER(S2);


% TIME FREQUENCY RESULTS - MAYBE SKIP
clear all
oat = osl_load_oat('sensorspace_tf.oat')
S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.modality='MEGPLANAR';
S2.first_level_contrast=[1];
S2.cfg.colorbar='yes';
S2.cfg.zlim = [0 35];
S2.view_cope=0;
[cfg, data]=osl_stats_multiplotTFR(S2);
title([oat.first_level.contrast_name{S2.first_level_contrast}]);

% Overview of sensor space analysis and stats testing there
%% OAT 1 - Sensorspace ERF Analysis
datadir = fullfile(osldir,'example_data','faces_singlesubject','spm_files');
workingdir = fullfile(osldir,'example_data','faces_singlesubject');
clear spm_files_continuous spm_files_epoched

spm_files_continuous{1}=[datadir '/Aface_meg1.mat'];
spm_files_epoched{1}=[datadir '/eAface_meg1.mat'];

oat=[];
oat.source_recon.D_epoched=spm_files_epoched; % this is passed in so that the bad trials and bad channels can be read out
oat.source_recon.D_continuous=spm_files_continuous;
oat.source_recon.conditions={'Motorbike','Neutral face','Happy face','Fearful face'};
oat.source_recon.freq_range=[4 100]; % frequency range in Hz
oat.source_recon.time_range=[-0.2 0.4];
oat.source_recon.method='none';
oat.source_recon.normalise_method='none';
oat.source_recon.dirname = [workingdir '/sensorspace_erf'];

Xsummary={};
Xsummary{1}=[1 0 0 0];Xsummary{2}=[0 1 0 0];Xsummary{3}=[0 0 1 0];Xsummary{4}=[0 0 0 1];
oat.first_level.design_matrix_summary=Xsummary;
oat.first_level.contrast={};
oat.first_level.contrast{1}=[3 0 0 0]'; % motorbikes
oat.first_level.contrast{2}=[0 1 1 1]'; % faces
oat.first_level.contrast{3}=[-3 1 1 1]'; % faces-motorbikes
oat.first_level.contrast_name{1}='motorbikes';
oat.first_level.contrast_name{2}='faces';
oat.first_level.contrast_name{3}='faces-motorbikes';

oat.first_level.cope_type='cope';
oat.first_level.report.first_level_cons_to_do=[2 1 3];
oat.first_level.bc=[0 0 0];

oat = osl_check_oat(oat);

oat
oat.source_recon
oat.first_level

oat.to_do=[1 1 0 0];
oat = osl_run_oat(oat);

disp(oat.results.report.html_fname); % show path to web page report


oat.first_level.report.modality_to_do='MEGPLANAR';
report = oat_first_level_stats_report(oat,oat.first_level.results_fnames{1});
open(oat.results.report.html_fname);
disp('oat.results:');
disp(oat.results);
stats1=oat_load_results(oat,oat.first_level.results_fnames{1});
figure;imagesc(stats1.x);title('GLM design matrix');xlabel('regressor no.');ylabel('trial no.');

S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.modality='MEGPLANAR'; % can also set this to 'MEGMAG'
S2.first_level_contrast=[3]; % view faces-motorbikes contrast
S2.view_cope=1; % set to 0 to see the t-stat
[cfg, dats, fig_handle]=oat_stats_multiplotER(S2);

% Overview of beamforming
% Set up data paths

Loaded from file  /Users/romesh/oxford_postdoc/toolboxes/osl/example_data/coreg_example/dsubject1_reduced.mat


datadir = fullfile(osldir,'example_data','faces_singlesubject','spm_files');
spm_files_continuous=[datadir '/Aface_meg1.mat'];
spm_files_epoched=[datadir '/eAface_meg1.mat'];
S = [];
S.D = spm_files_continuous;
S.mri = fullfile(osldir,'example_data','faces_singlesubject','structurals','struct.nii');
S.useheadshape = 1;
S.use_rhino = 0;
S.forward_meg = 'Single Shell';
S.fid.label.nasion = 'Nasion';
S.fid.label.lpa = 'LPA';
S.fid.label.rpa = 'RPA';
D=osl_headmodel(S);
S = [];
S.D = spm_files_epoched;
S.mri = fullfile(osldir,'example_data','faces_singlesubject','structurals','struct.nii');
S.useheadshape = 1;
S.use_rhino = 0;
S.forward_meg = 'Single Shell';
S.fid.label.nasion = 'Nasion';
S.fid.label.lpa = 'LPA';
S.fid.label.rpa = 'RPA';
D=osl_headmodel(S);

% Overview of coregistration
D = spm_eeg_load('rhino_example')
rhino_display(D);

% Overview of source space results and parcellation



% Some network plots