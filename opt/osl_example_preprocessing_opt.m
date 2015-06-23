%%%%%%%%%%%%%%%%%%
% This is a TEMPLATE script for running the OHBA recommended preprocessing
% pipeline on Elekta-Neuromag data (a very similar pipeline will work on
% CTF data as well) using OPT. It works through the following steps:
%
% You'll need to do alter (at the very least) the settings in:
% datadir, fif_files, spm_files, structurals_files

%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

global OSLDIR;
    
osldir = '/Users/andrew/Software/Matlab/osl2.0';

addpath(osldir);
osl_startup(osldir);

%%%%%%%%%%%%%%%%%%
%% INITIALISE GLOBAL SETTINGS FOR THIS ANALYSIS

% directory where the data is:
datadir = '/Users/andrew/Projects/OSL_test/osl2_tutorials/button_press_data_osl2';
 
% Set up the list of subjects and their structural scans for the analysis 
clear raw_fif_files input_files spm_files structural_files;

% Specify a list of the existing raw fif files for subjects for input into
% Maxfilter.
% Note that here we only have 1 subject, but more generally there would be
% more than one, e.g.:
% raw_files{1}=[testdir '/fifs/sub1_face_sss']; 
% raw_files{2}=[testdir '/fifs/sub2_face_sss']; 
% etc...
% OR
% Specify a list of the input files to be converted into SPM
% Note that here we only have 1 subject, but more generally there would be
% more than one, e.g.:
% fif_files{1}=[testdir '/fifs/sub1_face_sss']; 
% fif_files{2}=[testdir '/fifs/sub2_face_sss']; 
% etc...

raw_fif_files{1}=[datadir '/fifs/loc_S02.fif']; 
input_files{1}=[datadir '/fifs/loc_S02_sss1.fif']; 
spm_files{1}=[datadir '/spm_files/loc_S02']; 

% Setup a list of existing structural files, in the same order as spm_files and fif_files:
% Note that here we only have 1 subject, but more generally there would be
% more than one.
structural_files{1}=[datadir '/structs/anat.nii']; % leave empty if no .nii structural file available

%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP OPT

opt=[];

%%%%
% required inputs

%opt.raw_fif_files=raw_fif_files;
opt.input_files=input_files;
%opt.spm_files=spm_files;

opt.datatype='neuromag';

%%%%
% optional inputs

opt.dirname=[datadir '/practical_badseg_africa.opt']; % directory opt settings and results will be stored

% maxfilter settings
opt.maxfilter.do=0; % here we are going to skip the double maxfilter call as this has been run already for us
 
% africa settings
opt.africa.do=0;
opt.africa.ident.artefact_chans={'ECG','EOG'}; % artefact channels
opt.africa.ident.mains_kurt_thresh=0.5;

% highpass filter settings
opt.highpass.do=0;

% bad segments settings
opt.bad_segments.do=0;

% Epoching settings
%
% Here the epochs are set to be from -1s to +2s relative to triggers
% in the MEG data.
opt.epoch.do=1;
opt.epoch.time_range = [-1 2]; % epoch end in secs   
opt.epoch.trialdef(1).conditionlabel = 'StimLRespL';
opt.epoch.trialdef(1).eventtype = 'STI101_down';
opt.epoch.trialdef(1).eventvalue = 11;
opt.epoch.trialdef(2).conditionlabel = 'StimLRespR';
opt.epoch.trialdef(2).eventtype = 'STI101_down';
opt.epoch.trialdef(2).eventvalue = 16;
opt.epoch.trialdef(3).conditionlabel = 'StimRRespL';
opt.epoch.trialdef(3).eventtype = 'STI101_down';
opt.epoch.trialdef(3).eventvalue = 21;
opt.epoch.trialdef(4).conditionlabel = 'StimRRespR';
opt.epoch.trialdef(4).eventtype = 'STI101_down';
opt.epoch.trialdef(4).eventvalue = 26;
opt.epoch.trialdef(5).conditionlabel = 'RespLRespL'; %L but
opt.epoch.trialdef(5).eventtype = 'STI101_down';
opt.epoch.trialdef(5).eventvalue = 13;
opt.epoch.trialdef(6).conditionlabel = 'RespLRespR';
opt.epoch.trialdef(6).eventtype = 'STI101_down';
opt.epoch.trialdef(6).eventvalue = 19;
opt.epoch.trialdef(7).conditionlabel = 'RespRRespL'; % L but press
opt.epoch.trialdef(7).eventtype = 'STI101_down';
opt.epoch.trialdef(7).eventvalue = 23;
opt.epoch.trialdef(8).conditionlabel = 'RespRRespR';
opt.epoch.trialdef(8).eventtype = 'STI101_down';
opt.epoch.trialdef(8).eventvalue = 29;

% coreg settings
opt.coreg.do=0; % not doing coreg here - although normally you would if you want to do subsequent analyses in source space 

% outliers settings
opt.outliers.do=1;

%%%%
% setup settings

opt=osl_check_opt(opt);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOK AT OPT SETTINGS

disp('opt settings:');
disp(opt);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOK AT OPT SUB-SETTINGS

disp('opt.maxfilter settings:');
disp(opt.maxfilter);

disp('opt.downsample settings:');
disp(opt.downsample);

disp('opt.africa settings:');
disp(opt.africa);

disp('opt.highpass settings:');
disp(opt.highpass);

disp('opt.epoch settings:');
disp(opt.epoch);

disp('opt.outliers settings:');
disp(opt.outliers);

disp('opt.coreg settings:');
disp(opt.coreg);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN OPT

opt=osl_run_opt(opt);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESULTS

disp('opt.results:');
disp(opt.results);
