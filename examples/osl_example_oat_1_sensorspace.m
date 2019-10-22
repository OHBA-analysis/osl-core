%% OAT 1 - Sensorspace ERF Analysis
%
% In this practical we will work with a single subject's data from an
% emotional faces task (data courtesy of Susie Murphy) and perform an ERF
% analysis in sensor space.
%
% We will go through the following steps:
%
% # Prepare data for OAT analysis
% # Bandpass filter data and split into epochs
% # Compute a first level GLM analysis with OAT
% # Visualise results with FieldTrip
%
% Please read each cell in turn before copying it's contents either directly  
% into the MatLab console or your own blank script. By the end of this session 
% you should have created your own template analysis script which can be %
% applied to further analysis.
%
% We will work with a single subject's data from an emotional faces task and perform an ERF analysis in sensor space. 
% This will use the following files, which should be in the example_data
% directory:
%
% * Asss_fif_spm12_meg25.mat - an SPM MEEG object that has continuous data that has already been SSS Maxfiltered and downsampled to 250 Hz.
% * eAsss_fif_spm12_meg25.mat - an SPM MEEG object that has the same data epoched into the different task conditions.

%% A QUICK SUMMARY OF OAT
%
% OAT splits the analysis into 4 distinct pipeline stages. These are:
%
% # source reconstruction. Note that this stage is still used in a sensor space analysis in order to setup the data
% # first-level GLM trial-wise analysis
% # subject-level fixed effects averaging over multiple sessions
% # group-level GLM subject-wise analysis
%
% <<oat_stage_summary.png>>
%
% Since we are only analyzing a single session for a single subject here, we will only use the first two stages.

%% INITIALISE GLOBAL SETTINGS FOR THIS ANALYSIS
% This cell sets the directory that OAT will work in. Change the workingdir variable to correspond to the correct directory on your computer before running the cell.

% directory where the data is:
datadir = fullfile(osldir,'example_data','faces_singlesubject','spm_files');

% directory to put the analysis in
workingdir = fullfile(osldir,'example_data','faces_singlesubject');

%% SET UP THE SUBJECTS FOR THE ANALYSIS
%
% Specify a list of the fif files, structural files (not applicable for this practical) and SPM files (which will be created). It is important to make sure that the order of these lists is consistent across sessions. Note that here we only have 1 subject, but more generally there would be more than one. For example:
%
% |fif_files{1}=[testdir '/fifs/sub1_face_sss.fif'];|
%
% |fif_files{2}=[testdir '/fifs/sub2_face_sss.fif'];|
%
% etc...
%
% |spm_files{1} = [workingdir '/sub1_face_sss.mat'];|
%
% |spm_files{2} = [workingdir '/sub2_face_sss.mat'];|
%
% etc...

% clear old spm files
clear spm_files_continuous spm_files_epoched

spm_files_continuous{1} = fullfile(datadir,'Aface_meg1.mat');
spm_files_epoched{1}    = fullfile(datadir,'eAface_meg1.mat');

%% SETUP SENSOR SPACE SOURCE RECON
% This stage sets up the source reconstruction stage of an OAT analysis. The source_recon stage is always run even for a sensorspace analysis, though in these cases it simply prepares the data for subsequent analysis.
% In this example we define our input files (|D_continuous| and
% |D_epoched|) and conditions before setting a time frequency window from -200ms before stimulus onset to 400ms after and from 4Hz to 100Hz. The source recon method is set to 'none' as we are performing a sensorspace analysis.
% The |oat.source_recon.dirname| is where all the analysis will be stored. This includes all the intermediate steps, diagnostic plots and final results.

oat=[];
oat.source_recon.D_epoched=spm_files_epoched; % this is passed in so that the bad trials and bad channels can be read out
oat.source_recon.D_continuous=spm_files_continuous;
oat.source_recon.conditions={'Motorbike','Neutral face','Happy face','Fearful face'};
oat.source_recon.freq_range=[4 100]; % frequency range in Hz
oat.source_recon.time_range=[-0.2 0.4];
oat.source_recon.method='none';
oat.source_recon.normalise_method='none';

% Set this to something specific
oat.source_recon.dirname = fullfile(workingdir,'sensorspace_erf');

%% SETUP THE FIRST LEVEL GLM
% This cell defines the GLM parameters for the first level analysis.
% Critically this includes the design matrix (in Xsummary) and contrast matrix.
% Xsummary is a parsimonious description of the design matrix. It contains
% values |Xsummary{reg,cond}|, where reg is a regressor index number and cond is a condition index number. This will be used (by expanding the conditions over trials) to create the (num_regressors x num_trials) design matrix.
% 
% Each contrast is a vector containing a weight per condition defining how the condition 
% parameter estimates are to be compared. Each vector will produce a different t-map across the sensors. 
% Contrasts 1 and 2 describe positive correlations between each sensors activity and the 
% presence of a motorbike or face stimulus respectively. Contrast 3 tests whether each 
% sensors activity is larger for faces than motorbikes.

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
% The |osl_check_oat.m| function should be used to setup the settings for these three stages of the pipeline.  
%
% to see what the mandatory fields are. Note that you MUST specify:
%
% * |oat.source_recon.time_range|
% * |oat.source_recon.conditions|
% * |oat.source_recon.D_continuous|
% * |oat.source_recon.D_epoched|
%
% Descriptions of what these correspond to are also displayed when you type
% |help osl_check_oat| in to the matlab console.
%
% We are only interested in the first two stages of an OAT analysis, as we are not doing any group analysis here:
%
% * |oat.source_recon|
% * |oat.first_level|
%
% Each of these contains the settings for the relevant stages of the pipeline. You can read the Pipeline Stages