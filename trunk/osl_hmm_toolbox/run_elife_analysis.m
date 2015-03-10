% Runs group HMM analysis of 4-30 Hz band-limited source-space amplitude 
% envelopes based on the data and methods described in:
%
% "Fast transient networks in spontaneous human brain activity", eLife 2014
% Baker, A.P., Brookes, M.J., Rezek, I.A., Smith, S.M., Behrens, T., 
% Probert Smith, P.J., Woolrich, M.W. 
% DOI:10.7554/eLife.01867

% Beamformed data (4-30 Hz)

osldir = '/Users/abaker/Dropbox/Code/OSL2/osl2.0/';    
addpath(osldir);
osl_startup(osldir);

%%
spmdir      = '/Users/abaker/Data/notts_3state/spm12/';
spm_files = strcat([spmdir 'subj'],num2str(sort(repmat(1:10,1,3))','%0.2d'),'_',repmat({'ec','eo','mov'},1,10)');

sessions2use = setdiff(1:3:30,[4:6,7:9,10:12]);
%sessions2use = 1; %for test

spm_files = spm_files(sessions2use);

BFfiles = prefix(spm_files,'fA');
                  
hmmdir  = '/Users/abaker/Scratch/HMMtesting/HMMtest/';

% BFfiles = '/Users/abaker/Scratch/HMMtesting/epoched_testdata/concatfsession1_spm_meeg.mat';
% hmmdir  = '/Users/abaker/Scratch/HMMtesting/HMMtest_epoched/';



mkdir(hmmdir);



%% Run HMM with the default settings:
todo.envelope = 0;
todo.concat   = 0;
todo.infer    = 1;
todo.output   = 1;
    
options.envelope.windowsize = 0.1;
options.concat.log          = 0;
options.concat.norm_subs    = 1;
options.concat.pcadim       = 40;
options.hmm.nstates         = 8;
options.hmm.nreps           = 1;
options.output.method       = 'pcorr';

[HMMresults,statemaps] = osl_hmm_groupinference(BFfiles,hmmdir,todo,options);


%% View state maps
fslview(statemaps);


%% Plot temporal statistics
load(HMMresults)
stats = osl_hmm_stats(hmm,'do_plots');


%% Plot statepath
figure
osl_hmm_plotstatepath(hmm);


%%%%%%%%%%%%%%%%%%%%%%%%%
%% example parcellation call

todo.prepare  = 0;
todo.concat   = 1;
todo.infer    = 1;
todo.output   = 1;
   
options=[];
options.prepare.windowsize  = 0.1;
options.prepare.envelope    = 1;
options.prepare.log         = 0;

use_parcels=1;
if use_parcels
    %options.prepare.parcellation.file = fullfile([tilde '/parcellations/fmri_d100_parcellation_with_PCC_reduced_2mm_voted_2mm_alloc_8mm.nii.gz']);
    options.prepare.parcellation.file = '/Users/abaker/Data/fMRI_parcellations/giles_parcellation/fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz';
    options.prepare.parcellation.method = 'spatialBasis';       
    
    %options.prepare.parcellation.file = [tilde '/parcellations/fmri_d100_parcellation_with_PCC_reduced_2mm_ds8mm'];
    %options.prepare.parcellation.method = 'spatialBasis';
    
    %options.prepare.parcellation.protocol = 'none';
    options.prepare.parcellation.protocol = 'symmetric';
    
    options.concat.pcadim       = -1;
else
    options.concat.pcadim       = 40;

end;
options.concat.whiten       = 1;
options.concat.filename     = ['concat_pcdim' num2str(options.concat.pcadim)];
options.concat.normalisation = 'global';

options.hmm.nstates         = 8;
options.hmm.nreps           = 5;
options.hmm.use_old_hmm_tbx = 0;
options.hmm.filename        = [options.concat.filename '_hmm_NK' num2str(options.hmm.nstates)];

options.output.method       = 'pcorr';
options.output.filename     = options.hmm.filename;

%hmmdir  = ['/Users/abaker/Scratch/oxford_resting2/HMMtest/results_parcel' num2str(use_parcels) '_env' num2str(options.prepare.envelope) '_' options.prepare.parcellation.method '_' options.prepare.parcellation.protocol '/'];

% Run HMM with the default settings:
[HMMresults,statemaps,epoched_statepath_sub] = osl_hmm_groupinference_parcels(BFfiles,hmmdir,todo,options);

load(HMMresults)

%% View state maps
fslview(statemaps);


%% Plot temporal statistics
stats = osl_hmm_stats(hmm,'do_plots');


%% Plot statepath
figure
osl_hmm_plotstatepath(hmm);

%% Plot epoched statepath
figure
option=[];
option.mode='overlay';
option.epoched_statepath_sub=epoched_statepath_sub;
option.win=0.5;
osl_hmm_plotstatepath(hmm,option);



