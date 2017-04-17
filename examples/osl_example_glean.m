% GLEAN DEMO SCRIPT FOR INFERRING A GROUP HMM FROM SOURCE SPACE MEG DATA
% Adam Baker, Jan 2016
% Diego Vidaurre, Feb 2017
   
% Directory of the data
data_dir = fullfile(osldir,'example_data','glean_example');
% Name for this GLEAN analysis:
glean_name = fullfile(osldir,'example_data','glean_example','hmmmar_demo.mat');

do_analysis = 0; 

% List of data files. These data consist of 2 sessions of 10 minutes of 
% eyes open resting state data recorded from 6 subjects on an Elekta
% Neuromag system (using the 204 planar gradiometers only). The data were
% band-pass filtered between 4 and 30 Hz and projected onto a regular 8 mm
% grid in source space using an LCMV beamformer. The data are saved as
% SPM12 MEEG objects in sensor space, and the beamformer weights are saved
% as a virtual montage.
data = {'/data/s01_rest1.mat'
        '/data/s01_rest2.mat'
        '/data/s02_rest1.mat'
        '/data/s02_rest2.mat'
        '/data/s03_rest1.mat'
        '/data/s03_rest2.mat'
        '/data/s04_rest1.mat'
        '/data/s04_rest2.mat'
        '/data/s05_rest1.mat'
        '/data/s05_rest2.mat'
        '/data/s06_rest1.mat'
        '/data/s06_rest2.mat'};
data = fullfile(data_dir,data);

% Clear settings
settings = struct;

% Envelope settings:
settings.envelope.overwrite = 0;
settings.envelope.log       = 0;
settings.envelope.fsample   = 20;
settings.envelope.mask      = fullfile(data_dir,'MNI152_T1_8mm_brain.nii.gz');

% Subspace settings:
settings.subspace.overwrite                         = 0;
settings.subspace.normalisation                     = 'none';
settings.subspace.parcellation.file                 = fullfile(data_dir,'fMRI_parcellation_ds8mm.nii.gz');
settings.subspace.parcellation.orthogonalisation    = 'symmetric';
settings.subspace.parcellation.method               = 'spatialBasis';

% Model settings:
settings.model.overwrite   = 0;
settings.model.hmm.nstates = 8;
settings.model.hmm.nreps   = 1;

% Set up the GLEAN settings, data paths etc:
GLEAN = glean.setup(glean_name,data,settings);

% Run the analysis:
if do_analysis
    glean.run(GLEAN)
    save(glean_name,'GLEAN')
end

%% Load results and visualise them 

load(glean_name)

% show the state time courses
glean.plot_timecourse(GLEAN)

% compute temporal properties related to the estimation
settings = struct('plot',1);
GLEAN = glean.temporal_stats(GLEAN,settings);

% display the "temporal_stats"
GLEAN.results.temporal_stats
% and Fractional occupancy within it
GLEAN.results.temporal_stats.FractionalOccupancy

% open the plots for this statistic
open(GLEAN.results.temporal_stats.FractionalOccupancy.plots.stats)


%% Creating the spatial maps for each state 

% The function "glean.pcorr" will create spatial maps for each state, 
% by computing the partial correlation between each session's state time 
% courses and the envelope data at each voxel. These maps may be output 
% as .nii files (or alternatively as .mat files), and may be computed using 
% the voxelwise or parcelwise data:

settings            = [];
settings.format     = 'nii';
settings.space      = {'parcel'};
GLEAN = glean.pcorr(GLEAN,settings);


% GLEAN now contains a new field "temporal_stats" within GLEAN.results: 
GLEAN.results.pcorr

% This results field contains the settings, as well as spatial maps in 
% each subspace (voxel/parcel) for each session and the group average:
GLEAN.results.pcorr.parcel

fslview(GLEAN.results.pcorr.parcel.sessionmaps{1},[])

