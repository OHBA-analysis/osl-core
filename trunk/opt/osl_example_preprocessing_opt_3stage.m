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
%
% We are going to run opt in three stages.
%
% 1: maxfilter, convert, downsample, africa
% 2: manual africa, component removal
% 3: epoching and artefact detection
%
% These stages will complete all of the subprocesses for every datafile defined
% in the inputs in the previous cell. The only point which required manual
% intervention is the ICA artefact identification in stage 2. The africa step
% in stage 1 will automatically identfy some components which might be
% considered artefactual. That said, the question of which artefacts reject or
% not to reject is a complex one and often dependent on the particular dataset
% and experimental paradigm in question, as such it is best to take a look
% yourself at this point.
%
% If you want to trust the automatically selected artefactual components or are
% going to skip africa alltogether you may want to combine the three stages
% into 1.


%%%%%%%%%%%%%%%%%
%% Set up stage 1
%
% In this stage we are going to load in a set of raw fif files from a Neuromag
% dataset, convert them to spm format, downsample the data to 250Hz and run a
% preliminary ICA decomposition on the converted data.
%
% The opt struct requires at least the datatype and input_files however, for
% this example we are going to specify that OPT should only perform the stages
% above

opt=[];

% required inputs
opt.input_files=input_files; %post-maxfilter fif files
opt.datatype='neuromag';

% optional inputs

opt.dirname=[datadir '/practical_badseg_africa.opt']; % directory opt settings and results will be stored

% maxfilter settings
opt.maxfilter.do=0; % here we are going to skip the double maxfilter call as this has been run already for us
 
% africa settings
%
% There are three stages to africa
% 1: ica decomposition
% 2: artefact component identification
% 3: bad component removal
%
% For this stage we only want to run the first two stages, decomponsition and
% automatic identification. We specify this in the africa.todo struct.
% Automatic classification is the default so we don't have to add anything
% extra
opt.africa.todo.ica=1;
opt.africa.todo.ident=1;
opt.africa.todo.remove=0;
opt.africa.ident.artefact_chans={'ECG','EOG'}; % artefact channels
opt.africa.ident.mains_kurt_thresh=0.5;

% turn the remaining options off
opt.highpass.do=0;
opt.bad_segments.do=0;
opt.epoch.do=0;
opt.outliers.do=0;
opt.coreg.do=0;

%% Check stage 1
opt = osl_check_opt( opt );

%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOK AT OPT SETTINGS
%
% osl_check_opt will populate the opt struct with many default settings that
% you may later want to change. Take a look at the contents of opt in more
% detail, particularly the downsampling and africa structs

disp('opt settings:');
disp(opt);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOK AT OPT SUB-SETTINGS
%

disp('opt.maxfilter settings:');
disp(opt.maxfilter);

disp('opt.downsample settings:');
disp(opt.downsample);

disp('opt.africa settings:');
disp(opt.africa);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run stage 1
opt = osl_run_opt( opt );

%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESULTS
%
% Take a look at the html page generated by osl_run_opt, this contains some of
% the information collected within the results section of opt.

disp('opt.results:');
disp(opt.results);

%%%%%%%%%%%%%%%%%
%% Set up stage 2
%
% Stage 1 should have generated a downsampled spm object and a
% _africa_results.mat file in the specified opt dir. The africa components
% which have been automatically identified as artefactual will have also been
% plotted out. Take a look at these images in your optdir
%
% In this stage we are going to use a GUI to manually decide whether individual
% components should be removed from the dataset or not. To do this we input the
% spm_files from opt stage 1 and specify that africa should use the
% identify_artefactual_components_manual tool so we can directly intervene

opt2 = [];

% required inputs
opt2.spm_files = opt.results.spm_files;
opt2.datatype='neuromag';

% optional inputs
opt2.dirname=[datadir '/practical_badseg_africa.opt']; % directory opt settings and results will be stored

% africa settings
%
% In contrast to stage 1, we now was to run a manual identification followed by
% bad component removal. For this we specify that ident and remove should be
% completed within africa.todo and that the identification function should be
% @identify_artefactual_components_manual.
opt2.africa.todo.ica=0;
opt2.africa.todo.ident=1;
opt2.africa.todo.remove=1;
opt2.africa.ident.artefact_chans={'ECG','EOG'}; % artefact channels
opt2.africa.ident.mains_kurt_thresh=0.5;
opt2.africa.ident.func = @identify_artefactual_components_manual;

% turn the remaining options off
opt2.maxfilter.do=0;
opt2.convert.spm_files_basenames = opt.results.spm_files_basenames;
opt2.downsample.do=0;
opt2.highpass.do=0;
opt2.bad_segments.do=0;
opt2.epoch.do=0;
opt2.outliers.do=0;
opt2.coreg.do=0;

%% Check stage 2
opt2 = osl_check_opt( opt2 );

%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOK AT OPT SETTINGS

disp('opt2 settings:');
disp(opt2);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOK AT OPT SUB-SETTINGS

disp('opt2.africa settings:');
disp(opt2.africa);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run stage 2
%
% Now, when we run opt, a GUI window will open. This will contain several
% sections. The drop down menu specifies that we are looking at the components
% ranked by their VARIANCE, the bars below indicate how much variance each
% component has. The top bar will be blue, this is the component you can see in
% the top time-series to the left, the power-spectrum in the bottom left and
% the topo-plots to the bottomw right.
%
% You can cycle through the different components by pressing up and down and
% change what data feature they are ranked by by selecting a different feature
% from the box in the top-right.
%
% Some of the components will be marked in RED in the bars on the left and
% time-series. These are the components the automatic idenfication in stage 1
% deemed as artefactual. If you want to mark an additional bad component, cycle
% though to it using the up and down keys and press the red cross on the top
% left, the component will turn red. You can change your mind by pressing the
% green tick.
%
% Finally, if you select EOG or ECG from the box in the top-right, you will see
% the COVARIATE plotted in the middle. This is the time-series which the ICA
% components are correlated with (to generate the bars on the left). This is a
% good sanity check to make sure that any components you reject do strongly
% resemble the covariate channel and that the covarate channel itself is clean
% enough to trust
%
% When you have finished making your selections, close the GUI and your bad
% channels will be removed from the data

opt2 = osl_run_opt( opt2 );

%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESULTS
%
% Take a look at the html page generated by osl_run_opt, this contains some of
% the information collected within the results section of opt.

disp('opt2.results:');
disp(opt2.results);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up stage 3
%
% Stage 2 should have generated an Adspm_meg1 file in your optdir. The A
% indicates that the components you selected in stage 2 have now been removed
% from the data.
%
% The final stage will perform a highpass filter before epoching the data and
% automatically identfy channels and epochs which are outliers and should be
% removed from the analysis

opt3 = [];

% required inputs
opt3.spm_files = opt2.results.spm_files;
opt3.datatype='neuromag';

% optional inputs
opt2.dirname=[datadir '/practical_badseg_africa.opt']; % directory opt settings and results will be stored

% Epoching settings
%
% Here the epochs are set to be from -1s to +2s relative to triggers
% in the MEG data.
opt3.epoch.do=1;
opt3.epoch.time_range = [-1 2]; % epoch end in secs   
opt3.epoch.trialdef(1).conditionlabel = 'StimLRespL';
opt3.epoch.trialdef(1).eventtype = 'STI101_down';
opt3.epoch.trialdef(1).eventvalue = 11;
opt3.epoch.trialdef(2).conditionlabel = 'StimLRespR';
opt3.epoch.trialdef(2).eventtype = 'STI101_down';
opt3.epoch.trialdef(2).eventvalue = 16;
opt3.epoch.trialdef(3).conditionlabel = 'StimRRespL';
opt3.epoch.trialdef(3).eventtype = 'STI101_down';
opt3.epoch.trialdef(3).eventvalue = 21;
opt3.epoch.trialdef(4).conditionlabel = 'StimRRespR';
opt3.epoch.trialdef(4).eventtype = 'STI101_down';
opt3.epoch.trialdef(4).eventvalue = 26;
opt3.epoch.trialdef(5).conditionlabel = 'RespLRespL'; %L but
opt3.epoch.trialdef(5).eventtype = 'STI101_down';
opt3.epoch.trialdef(5).eventvalue = 13;
opt3.epoch.trialdef(6).conditionlabel = 'RespLRespR';
opt3.epoch.trialdef(6).eventtype = 'STI101_down';
opt3.epoch.trialdef(6).eventvalue = 19;
opt3.epoch.trialdef(7).conditionlabel = 'RespRRespL'; % L but press
opt3.epoch.trialdef(7).eventtype = 'STI101_down';
opt3.epoch.trialdef(7).eventvalue = 23;
opt3.epoch.trialdef(8).conditionlabel = 'RespRRespR';
opt3.epoch.trialdef(8).eventtype = 'STI101_down';
opt3.epoch.trialdef(8).eventvalue = 29;

% coreg settings
opt.coreg.do=0; % not doing coreg here - although normally you would if you want to do subsequent analyses in source space 

% highpass filter
opt3.highpass.do=1;

% outliers settings
opt3.outliers.do=1;

% turn the remaining options off
opt3.maxfilter.do=0;
opt3.convert.spm_files_basenames = opt2.results.spm_files_basenames;
opt3.downsample.do=0;
opt3.africa.todo.ica=0;
opt3.africa.todo.ident=0;
opt3.africa.todo.remove=0;
opt3.bad_segments.do=0;
% We are skipping coreg here but would normally include it for subsequent source analysis
opt3.coreg.do=0; 


%% Check stage 3
opt3 = osl_check_opt( opt3 );

%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOK AT OPT SETTINGS

disp('opt settings:');
disp(opt3);


%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOK AT OPT SUB-SETTINGS

disp('opt3.highpass settings:');
disp(opt3.highpass);

disp('opt3.epoch settings:');
disp(opt3.epoch);

disp('opt3.outliers settings:');
disp(opt3.outliers);

disp('opt3.coreg settings:');
disp(opt3.coreg);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN OPT
%
% Take a look at the html page generated by osl_run_opt, this contains some of
% the information collected within the results section of opt.

opt3=osl_run_opt(opt3);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESULTS
%
% The outlier option in stage 3 will automatically assess each of your channels
% and epochs to see if they could be considered outliers and removed from the
% analysisi. This is done using the distribution of standard deviations or
% minimum values across channels or epochs. These distributions can be seen in
% top half of the OUTLIERS plots. The epochs/channels considered to be normal
% are circled in green, the remaining observations will be removed from the
% analysis. The 'cleaned' distributions can be seen in the bottom half of the
% OUTLIERS plots, if there were large artefacts these should now look closer to
% a normal distribution.

disp('opt3.results:');
disp(opt3.results);
