%% OSL_example_oil
%% OHBA's ICA pipeLine (OIL) example script

% The OIL Pipeline is as follows:
% 1.) Feed in beamformed MEG data
% 2.) Envelope estimation and down-sampling foolowed by spatial smoothing.
% 3). Concatenation of all the data in time
% 4). Temporal ICA
% 5). 1st Level Stats
% 6). Group Stats
% 7.) Do Multiple Comparisons Corrections
% 8.) Generate Nifti Maps
%
% You must have setup osl to run this script

% Henry Luckhoo
% henry.luckhoo@trinity.ox.ac.uk
% 13.08.12

% References:
% Brookes et al. 2011 PNAS - Resting State Analysis
% Luckhoo et al. 2012 Neuroimage - Task-postive analysis


%%%%%%%%%%%%%%%%%%
%% INITIALISE GLOBAL SETTINGS FOR THIS ANALYSIS
% directory where the source recon data are - YOU NEED TO CHANGE THIS PATH
% files from this analysis will be nested within this directory.
datadir =  '/home/gileslc/data/eyesopen-eyesclosed/inverse-LCMV-4_30Hz-8mm/'; 

% load in the source recon information from an oat first stage
oatFile = fullfile(datadir, 'oat', 'oat.mat');
tmp = load(oatFile);
SourceRecon = tmp.oat.source_recon;

%% Setup exploratory ICA analysis
% The pipeline is controlled by the to_do flag. 
% Each element is a switch for steps:
%   1. Enveloping
%   2. Concatenation
%   3. ICA
%   4. Subject statistics and spatial maps
%   5. Group level statistics

% setup a name for where to save the OIL analysis settings structure.
oilSaveFile = fullfile(datadir, 'oil-test.mat');

oil              = struct();
oil.fname        = oilSaveFile;
oil.paradigm     = 'rest';
oil.source_recon = SourceRecon;

% Enveloping Parameters
oil.enveloping.gridstep      = 8; % mm
oil.enveloping.window_length = 2; % s

% Run enveloping and concatenation
oil.to_do = [1 1 0 0 0];
osl_save_oil(oil);
oil       = osl_run_oil(oil);
osl_save_oil(oil);

%% Part 3 - Concatenation and ICA for a larger data set - Set Up

% Running the beamforming and enveloping at sensible resolution for a group
% ICA would take too long for this practical. You will now load an OIL
% structure where the beamforming and enveloping have already been done for
% 6 subjects, each with a 10 min session of eyes open resting-state and eyes closed 
% resting-state. 

% You will now run the concatenation and temporal ICA for 12 session
% The beamforming and enveloping has already been done for a large group at 
% reasonable resolution. Here we will run the concatenation and temporal ICA. 
% In this cell we set up the analysis. Run the cell and then take two minutes 
% to look at the different setting. In most cases the default parameters should 
% be used. However the user must think carefully about whether to use
% temporal or
% spatial ICA, how many components to use and whether they want to use
% ICASSO.

oil.ica.num_ics      = 25;
oil.ica.temp_or_spat = 'temporal';
oil.to_do = [0 0 1 0 0];
oil       = osl_run_oil(oil);

%% Part 5 - Viewing the ICA maps.

% One way to look at the spatial extent of an independent component is to
% consider the covariance between each independent time course and each
% voxel of the concatenated data. OIL produces these normalised maps.

% Run this cell to load the map into FSLVIEW. Add the negative colours in
% "blue-light blue". This cannot be done until FSL has been launched. Note
% that OSL includes a wrapper function for launching FSL, loading multiple
% files, setting positive colour maps and thresholds.

ica_covariance_maps = oil.ica.results.maps;

fslview({[workingdir 'fmri_visual_network.nii.gz'];  ica_covariance_maps},[10 20; 2 5],{'Green';'Red-Yellow'})

%% Part 5 cont.

% Have a quick look through the 25 components. Pick the one or two
% components which correspond to the fMRI mask of a visual network. We will
% limit our group analysis to these components. Modify the vector below to
% include the "VOLUME NUMBER" from FSL for the visual components. 

visual_components = [8]+1; % Please identify any visual network components.

% Note that the "+1" is needed as FSL indexes volumes from 0 where as
% MatLab will index them from 1.

%% Part 6 - Statistical Analysis using OIL - Set Up

% Now we will try and detect any significant differences between the visual
% networks in the eyes closed and eyes open states. In principle, we expect
% increased alpha power/activity in the eyes closed state.

% For this analysis we don't actually need to modify anything for the first
% level stage. Stage 5 will automatically estimate the single subject
% component COPE maps using Dual Regression of the ICA maps we just
% estimated from the non-weights-normalised data.

% However, we do need to think carefully about what group comparisons we
% want to make. In this case, the most powerful test we can do is a paired
% t-test looking at the difference between the eyes open and eyes closed
% sessions in each subject. This is implemented with the following design
% matrix and contrast.

oil.ica_group_level.group_design_matrix = [
    1 -1 1 -1 1 -1 1 -1 1 -1 1 -1;...
    1  1 0  0 0  0 0  0 0  0 0  0;...
    0  0 1  1 0  0 0  0 0  0 0  0;...
    0  0 0  0 1  1 0  0 0  0 0  0;...
    0  0 0  0 0  0 1  1 0  0 0  0;...
    0  0 0  0 0  0 0  0 1  1 0  0;...
    0  0 0  0 0  0 0  0 0  0 1  1;...
    ]';

oil.ica_group_level.group_contrast = {};
oil.ica_group_level.group_contrast{1}=[1 0 0 0 0 0 0]';

% In addition, we will want to do some variance smoothing. Finally, to
% account for any multiple comparisons and inflated statistics we will use
% RANDOMISE to test for significant clusters. Note that we will limit the
% analysis to the components listed in "visual_components".

oil.ica_group_level.group_varcope_spatial_smooth_fwhm = 8;
oil.ica_group_level.use_randomise = 1;
oil.ica_group_level.Npermutations = 500;
oil.ica_group_level.cluster_threshold = 2.3;
oil.ica_group_level.comps2use = visual_components;

oil.to_do = [0 0 0 0 1 1];

disp('The first level parameters... oil.ica_first_level contains')
disp(oil.ica_first_level);

disp('The group level parameters... oil.ica_group_level contains')
disp(oil.ica_group_level);

%% Part 7 - Statistical Analysis using OIL - Running the analysis

oil= osl_run_oil(oil);

disp('The first level results... oil.ica_first_level.results contains')
disp(oil.concat_subs.results);

disp('The group level results... oil.ica_group_level.results contains')
disp(oil.ica.results);

%% Part 8 - Viewing the statistics results

% We will now view the results in FSLVIEW. THe following command will load
% in the ICA covariance maps, the super-threshold t-stat maps and the
% corresponding p-values.

% You need to set the negative values of the ICA Covariance maps to
% "blue light-blue" and set the thresholds to [2 5].
% You need to set the negative values of the t-stat maps to "yellow" and set the thresholds to [2.3 5].
% You need to set the p-value maps to be partially transparent.

fslview({oil.ica.results.maps; oil.ica_group_level.results.clustered_tstats_names{1}; oil.ica_group_level.results.corr_1minusp_names{1}},[2 5; 2.3 5; 0 1],{'Red-Yellow'; 'Green'; 'Pink'})

% Now look at the visual components.

% Do we see sensible network differences?
% Are they in a sensible direction?
% What needs to be improved in our analysis?

%% OPTIONAL Part 9 - Using an alternative basis set (e.g. fMRI).

% Here we will use the fMRI maps generated by a 20 component spatial ICA as
% a spatial basis set instead of our MEG ICA spatial maps. We will then run
% the same statistical analysis but only look at three visual network
% components.

oil.ica_first_level.spatial_basis_set =[workingdir 'fmri_20ICs.nii.gz'] ;

oil.ica_group_level.group_design_matrix = [
    1 -1 1 -1 1 -1 1 -1 1 -1 1 -1;...
    1  1 0  0 0  0 0  0 0  0 0  0;...
    0  0 1  1 0  0 0  0 0  0 0  0;...
    0  0 0  0 1  1 0  0 0  0 0  0;...
    0  0 0  0 0  0 1  1 0  0 0  0;...
    0  0 0  0 0  0 0  0 1  1 0  0;...
    0  0 0  0 0  0 0  0 0  0 1  1;...
    ]';

oil.ica_group_level.group_contrast = {};
oil.ica_group_level.group_contrast{1}=[1 0 0 0 0 0 0]';

% In addition, we will want to do some variance smoothing. Finally, to
% account for any multiple comparisons and inflated statistics we will use
% RANDOMISE to test for significant clusters. Note that we will limit the
% analysis to the components listed in "visual_components".

oil.ica_group_level.group_varcope_spatial_smooth_fwhm = 16;
oil.ica_group_level.use_randomise = 1;
oil.ica_group_level.Npermutations = 500;
oil.ica_group_level.cluster_threshold = 'TFC'; % Let's try threshold-free clustering for funsies!
oil.ica_group_level.comps2use = [5 6 9]; % We will limit ourselves to looking at the visual networks (3 of them).

oil.to_do = [0 0 0 0 1 1];

oil= osl_run_oil(oil);

fslview({[workingdir 'fmri_20ICs.nii.gz']; oil.ica_group_level.results.clustered_tstats_names{1}; oil.ica_group_level.results.corr_1minusp_names{1}},[10 20; 800 1500; 0.7 1],{'Red-Yellow'; 'Green'; 'Pink'})

% Adjust the colour bars and thresholds as necessary and look at volumes 4,5 and 8.