% In this practical we will work with a single subject's data from an
% emotional faces task (data courtesy of Susie Murphy) and perform an ERF
% analysis in sensor space. This dataset can be downloaded from:  
% 
% www.fmrib.ox.ac.uk/~woolrich/faces_subject1_data.tar.gz
% 
% Note that this contains the spm file:
% spm8_meg1.mat
% that is an SPM MEEG object that has continuous data that has already been
% SSS Maxfiltered and downsampled to 250 Hz. 
% 
% and
% espm8_meg1.mat
%
% which is an SPM MEEG object that has the same data epoched into the
% different task conditions.

%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

OSLDIR = getenv('OSLDIR');
    
% set this to where you have downloaded OSL and the practical data:
practical_dir='/home/mwoolrich/Desktop'; 
osldir=[practical_dir '/osl2.0'];    

practical_dir='/Users/woolrich';
osldir=[practical_dir '/homedir/matlab/osl2.0'];    

addpath(osldir);
osl_startup(osldir);

%%%%%%%%%%%%%%%%%%
%% INITIALISE GLOBAL SETTINGS FOR THIS ANALYSIS

% directory where the data is:
workingdir=[practical_dir '/homedir/matlab/osl2_testdata_dir/faces_subject1_data_osl2']; % directory where the data is

cmd = ['mkdir ' workingdir]; if ~exist(workingdir, 'dir'), unix(cmd); end % make dir to put the results in

clear spm_files_continuous spm_files_epoched;

% set up a list of SPM MEEG object file names (we only have one here)
spm_files_continuous{1}=[workingdir '/dspm_meg1.mat'];
spm_files_epoched{1}=[workingdir '/Sedspm_meg1.mat'];

%%%%%%%%%%%%%%%%%%%
%% SETUP SENSOR SPACE AND FIRST LEVEL GLM USING OAT

oat=[];
oat.source_recon.D_continuous=spm_files_continuous;
oat.source_recon.conditions={'Motorbike','Neutral face','Happy face','Fearful face'};
oat.source_recon.D_epoched=spm_files_epoched; % this is passed in so that the bad trials and bad channels can be read out
oat.source_recon.freq_range=[4 100]; % frequency range in Hz
oat.source_recon.time_range=[-0.2 0.4];
oat.source_recon.method='none';
oat.source_recon.normalise_method='none';

% Xsummary is a parsimonious description of the design matrix.
% It contains values Xsummary{reg,cond}, where reg is a regressor no. and cond
% is a condition no. This will be used (by expanding the conditions over
% trials) to croat_settingse the (num_regressors x num_trials) design matrix:
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

%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOK AT OAT SETTINGS
% As well as using the settings we specified in the previous cell, calling osl_check_oat  has filled in a bunch of other settings as well.
oat
oat.source_recon
oat.first_level

%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN OAT
oat.to_do=[1 1 0 0];
oat = osl_run_oat(oat);

%%

oat.first_level.report.modality_to_do='MEGPLANAR';
report = oat_first_level_stats_report(oat,oat.first_level.results_fnames{1});

%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESULTS

disp('oat.results:');
disp(oat.results);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% View the GLM design matrix 
% (NOTE that column 1 is motorbikes, columns 2-4 are faces)

% load first-level GLM result
stats1=oat_load_results(oat,oat.first_level.results_fnames{1});

figure;imagesc(stats1.x);title('GLM design matrix');xlabel('regressor no.');ylabel('trial no.');

%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualise using Fieldtrip
% note that this produces an interactive figure, with which you can:
% - draw around a set of sensors
% - click in the drawn box to produce a plot of the time series
% - on the time series plot you can draw a time window
% - and click in the window to create a topoplot averaged over that time
% window (which is itself interactive....!)

S2=[];
S2.oat=oat;
S2.stats_fname=oat.first_level.results_fnames{1};
S2.modality='MEGMAG'; % can also set this to 'MEGPLANAR'
S2.first_level_contrast=[2];
S2.view_cope=1;

% calculate t-stat using contrast of absolute value of parameter estimates
[cfg, dats, fig_handle]=oat_stats_multiplotER(S2);
