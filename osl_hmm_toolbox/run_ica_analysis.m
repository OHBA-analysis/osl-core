% Runs group ICA analysis of 4-30 Hz band-limited source-space amplitude 
% envelopes based on the data and methods described in:
%
% "Fast transient networks in spontaneous human brain activity", eLife 2014
% Baker, A.P., Brookes, M.J., Rezek, I.A., Smith, S.M., Behrens, T., 
% Probert Smith, P.J., Woolrich, M.W. 
% DOI:10.7554/eLife.01867

% Beamformed data (4-30 Hz)

BFfiles = ABgetfilelist('/Users/abaker/Scratch/oxford_resting2/','fAOxRest*.mat');
           
        
icadir  = '/Users/abaker/Scratch/oxford_resting2/ICA/';
mkdir(icadir);

osldir = '/Users/abaker/Dropbox/Code/OSL2/osl2.0/';    
addpath(osldir);
osl_startup(osldir);

%% Run HMM with the default settings:
todo.envelope = 1;
todo.concat   = 1;
todo.infer    = 1;
todo.output   = 1;
    
options.envelope.windowsize = 1;
options.concat.norm_subs    = 1;
options.concat.pcadim       = 20;
options.ica = [];
options.output.method       = 'pcorr';

[ICAresults,statemaps] = AB_runICA2(BFfiles,icadir,todo,options);


%% View state maps
fslview(statemaps);


%% Plot temporal statistics
load(HMMresults)
stats = osl_hmm_stats(hmm,'do_plots');


%% Plot statepath
figure
osl_hmm_plotstatepath(hmm);

