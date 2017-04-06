%% Sensorspace ERF Analysis with OAT
%
% In this practical we will work with a single subject's data from an
% emotional faces task (data courtesy of Susie Murphy) and perform an ERF
% analysis in sensor space.
%
% You will need the following files from the example_data directory:
%
% * Asss_fif_spm12_meg25.mat - an SPM MEEG object that has continuous data that has already been SSS Maxfiltered and downsampled to 250 Hz.
% * eAsss_fif_spm12_meg25.mat - an SPM MEEG object that has the same data epoched into the different task conditions.

%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

setenv('OSLDIR','/Users/andrew/Software/Matlab/osl2.1/')
addpath(genpath(getenv('OSLDIR')));
osl_startup(osldir);

%% INITIALISE GLOBAL SETTINGS FOR THIS ANALYSIS
% This cell sets the directory that OAT will work in. Change the workingdir variable to correspond to the correct directory on your computer before running the cell.

% directory where the data is:
workingdir='/Users/andrew/Projects/OSL2_testdir/meg_workshop/drugface/data/'; % directory where the data is

cmd = ['mkdir ' workingdir]; if ~exist(workingdir, 'dir'), unix(cmd); end % make dir to put the results in

clear spm_files_continuous spm_files_epoched;

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

spm_files_continuous{1}=[workingdir '/Asss_fif_spm12_meg25.mat'];
spm_files_epoched{1}=[workingdir '/eAsss_fif_spm12_meg25.mat'];

%% SETUP SENSOR SPACE SOURCE RECON
% This stage sets up the source reconstruction stage of an OAT analysis. The source_recon stage is always run even for a sensorspace analysis, though in these cases it simply prepares the data for subsequent analysis.
% In this example we define our input files (D_continuous and D_epoched) and conditions before setting a time frequency window from -200ms before stimulus onset to +400ms and 4 to 100Hz. The source recon method is set to 'none' as we are performing a sensorspace analysis
% The oat.source_recon.dirname is where all the analysis will be stored. This includes all the intermediate steps, diagnostic plots and final results. Make sure this directory does not contain any other analyses that might be overwritten.

oat=[];
oat.source_recon.D_epoched=spm_files_epoched; % this is passed in so that the bad trials and bad channels can be read out
oat.source_recon.D_continuous=spm_files_continuous;
oat.source_recon.conditions={'Motorbike','Neutral face','Happy face','Fearful face'};
oat.source_recon.freq_range=[4 100]; % frequency range in Hz
oat.source_recon.time_range=[-0.2 0.4];
oat.source_recon.method='none';
oat.source_recon.normalise_method='none';

% Set this to something specific
oat.source_recon.dirname = 'sensorspace_erf';

%% SETUP THE FIRST LEVEL GLM
% This cell defines the GLM parameters for the first level analysis. Critically this includes the design matrix (in Xsummary) and contrast matrix
% Xsummary is a parsimonious description of the design matrix. It contains values Xsummary{reg,cond}, where reg is a regressor index number and cond is a condition index number. This will be used (by expanding the conditions over trials) to croat_settingse the (num_regressors x num_trials) design matrix:
% Each contrast is a vector containing a weight per condition defining how the condition parameter estimates are to be compared. Each vector will produce a different t-map across the sensors. Contrasts 1 and 2 describe positive correlations between each sensors activity and the presence of a motorbike or face stimulus respectively. Contrast 3 tests whether each sensors activity is larger for faces than motorbikes.

Xsummary={};
Xsummary{1}=[1 0 0 0];Xsummary{2}=[0 1 0 0];Xsummary{3}=[0 0 1 0];Xsummary{4}=[0 0 0 1];
oat.first_level.design_matrix_summary=Xsummary;

% contrasts to be calculated:
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

%% VIEW OAT SETTINGS
% As well as using the settings we specified in the previous cell, calling osl_check_oat has filled in a bunch of other settings as well.

oat
oat.source_recon
oat.first_level

%% CHECK OVER OAT SETTINGS
% The osl_check_oat.m function should be used to setup the settings for these three stages of the pipeline.  On the Matlab command line type
%
% help osl_check_oat
%
% to see what the mandatory fields are. Note that you MUST specify:
%
% * oat.source_recon.time_range
% * oat.source_recon.conditions
% * oat.source_recon.D_continuous
% * oat.source_recon.D_epoched
%
% Descriptions of what these correspond to are also displayed when you type >> help osl_check_oat
%
% We are only interested in the first two stages of an OAT analysis, as we are not doing any group analysis here:
%
% * oat.source_recon
% * oat.first_level
%
% Each of these contains the settings for the relevant stages of the pipeline. You can read the Pipeline Stagesù section of the OSL manual to get more description about these two stages.
%
% Take a look at oat.source_recon.dirname. This will be the name of the directory (full path) where the OAT will be stored, and is given a .oatextension. [Note that each OAT directory is associated with a specific source recon; a new source recon will overwrite an old directory if the same oat.source_recon.dirname is used, and the old source recon results will be lost. Hence, you should ensure that you change oat.source_recon.dirname for a new source recon analysis, if you want to avoid overwriting an old one.]

%% RUN OAT
%
% The OAT structure you have created, oat, should be passed to the function osl_run_oat.m, to run the pipeline.
% However, the OAT structure should contain a structure (oat.to_do), which is a list of binary values indicating which part of the pipeline is to be run.
% E.g. oat.to_do=[1 1 0 0]; will run just the source recon and first level stages, whereas oat.to_do=[1 1 1 1]; will run all four (source-recon, first level, subject-level and group level).
% Here we are not doing any group analysis, so we need:
%
% oat.to_do=[1 1 0 0];

oat.to_do=[1 1 0 0];
oat = osl_run_oat(oat);

%%

oat.first_level.report.modality_to_do='MEGPLANAR';
report = oat_first_level_stats_report(oat,oat.first_level.results_fnames{1});

%% RESULTS
%
% The results are stored in the oat structure and the can be loaded back into matlab using oat_load_results. This is useful for checking over results of the GLM or as the basis for further analyses.

disp('oat.results:');
disp(oat.results);

% load first-level GLM result
stats1=oat_load_results(oat,oat.first_level.results_fnames{1});

%% VIEW THE GLM DESIGN MATRIX
%
% This is the nconditions by ntrials design matrix which is used to fit the GLM
% (NOTE that column 1 is motorbikes, columns 2-4 are faces)

figure;imagesc(stats1.x);title('GLM design matrix');xlabel('regressor no.');ylabel('trial no.');

%% VISUALISE USING FIELDTRIP
% note that this produces an interactive figure, with which you can:
% * draw around a set of sensors
% * click in the drawn box to produce a plot of the time series
% * on the time series plot you can draw a time window
% * and click in the window to create a topoplot averaged over that time
% window (which is itself interactive....!)

S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.modality='MEGPLANAR'; % can also set this to 'MEGMAG'
S2.first_level_contrast=[3]; % view faces-motorbikes contrast
S2.view_cope=1; % set to 0 to see the t-stat

% calculate t-stat using contrast of absolute value of parameter estimates
[cfg, dats, fig_handle]=oat_stats_multiplotER(S2);
