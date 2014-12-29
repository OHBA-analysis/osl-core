% Runs group HMM analysis of 4-30 Hz band-limited source-space amplitude 
% envelopes based on the data and methods described in:
%
% "Fast transient networks in spontaneous human brain activity", eLife 2014
% Baker, A.P., Brookes, M.J., Rezek, I.A., Smith, S.M., Behrens, T., 
% Probert Smith, P.J., Woolrich, M.W. 
% DOI:10.7554/eLife.01867

% Beamformed data (4-30 Hz)


spmdir      = '/Users/abaker/Data/notts_3state/spm12/';
spm_files = strcat([spmdir 'subj'],num2str(sort(repmat(1:10,1,3))','%0.2d'),'_',repmat({'ec','eo','mov'},1,10)');

%sessions2use = setdiff(1:30,[3:3:30,4:6,7:9,10:12]);
sessions2use = setdiff(1:3:30,[4:6,7:9,10:12]);
spm_files = spm_files(sessions2use);

BFfiles = prefix(spm_files,'fA');
                  
%hmmdir  = '/Users/abaker/Scratch/oxford_resting2/HMM/';
hmmdir  = '/Users/abaker/Data/notts_3state/spm12/HMM';

mkdir(hmmdir);

osldir = '/Users/abaker/Dropbox/Code/OSL2/osl2.0/';    
addpath(osldir);
osl_startup(osldir);

%% Run HMM with the default settings:
todo.envelope = 0;
todo.concat   = 0;
todo.infer    = 0;
todo.output   = 1;
    
options.envelope.windowsize = 0.1;
options.concat.log          = 1;
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

