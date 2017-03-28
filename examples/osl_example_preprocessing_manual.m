% In this practical/template script we will work with a single subject's data from an
% emotional faces task (data courtesy of Susie Murphy). This can be
% downloaded from:  
% 
% www.fmrib.ox.ac.uk/~woolrich/faces_subject1_data.tar.gz
% 
% Note that this contains the fif file:
% fifs/sub1_face_sss.fif 
% that has already been SSS Maxfiltered and downsampled to 250 Hz.
% 
% In this example we will take this fif file and run it through a
% manual preprocessing pipeline
%
% MWW 2013


% ROBERT:
% COMMENT: at which point is there an overview of standard prefixes?
% COMMENT: for the two preproc practicals - just ONE data set maybe?
% artefact_OSL marker - only used within OSL?

% AFRICA: not for automated.
% visualize hierarchy of tools

%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path

%tilde='/home/mwoolrich/Desktop';
tilde='/Users/robert/';

osldir=['/Applications/osl/osl2'];    

%osldir = getenv('OSLDIR');
%data_dir = fullfile(osldir,'example_data','roinets_example');
%output_directory = fullfile(osldir,'practical','roinets_demo');


addpath(osldir);

osl_startup();

%%%%%%%%%%%%%%%%%%
%% SPECIFY DIRS FOR THIS ANALYSIS

% directory where the data is:
datadir=[tilde '/Documents/workshop/example/faces_subject1_data'];

fullfile(osldir,'example_data','roinets_data')

% this is the directory the analysis files will be stored in:
workingdir=[datadir]; 
cmd = ['mkdir ' workingdir]; unix(cmd); % make dir to put the results in
 
%%%%%%%%%%%%%%%%%%
%% Set up the list of subjects and their structural scans for the analysis  
clear fif_files spm_files_basenames;

% Specify a list of the existing fif files for subjects
% Note that here we only have 1 subject, but more generally there would be
% more than one, e.g.:
% fif_files{1}=[testdir '/fifs/sub1_face_sss.fif']; 
% fif_files{2}=[testdir '/fifs/sub2_face_sss.fif']; 
% etc...
fif_files{1}=[datadir '/fifs/sub1_face_sss.fif']; 

% Setup a list of SPM MEEG object file names to be created, in the same order as spm_files and fif_files:
% Note that here we only have 1 subject, but more generally there would be
% more than one, e.g.:
% spm_files{1}=[workingdir '/spm8_meg1.mat'];
% spm_files{2}=[workingdir '/spm8_meg1.mat'];
% etc...
spm_files_basenames{1}=['spm_meg1.mat'];

%%%%%%%%%%%%%%%%%%%%
%% CONVERT FROM FIF TO AN SPM MEEG OBJECT:
% The fif file that we are workingedi with is sub1_face_sss.fif. This has
% already been maxfiltered for you and downsampled to 250Hz.
% 
% This will produce a histogram plot showing the number of events detected
% for each code on the trigger channel. The codes used on the trigger
% channel for this experiment were: 
%
% 1 = Neutral face
% 2 = Happy face
% 3 = Fearful face
% 4 = Motorbike
% 18 = Introduction screen
% 11 = Break between blocks
% 19 = Midway break
% 12 = Green fixation cross (response trials)
% 13 = Red fixation cross (following green on response trials)
% 14 = Red fixation cross (non-response trials)
%
% For example, there should be 120 motorbike trials, and 80 of each of the
% face conditions
 for subnum = 1:length(fif_files), % iterates over subjects
    spm_files{subnum}=[workingdir '/' spm_files_basenames{subnum}];
end

if(length(fif_files)>0),
    S2=[];
    for i=1:length(fif_files), % loops over subjects
        
        S2.fif_file=fif_files{i};
        S2.spm_file=spm_files{i};       
        S2.trigger_channel_mask='0000000000111111'; % binary mask to use on the trigger channel
        
        % The conversion to SPM will show a histogram of the event codes
        % and correspond to those listed below in the epoching section
        [D spm_files{i}] = osl_convert_script(S2);
    end;
end;

% Note that this spmfile is the output from the conversion:
spm_files{1}

%%%%%%%%%%%%%%%%%%%%
%% LOAD THE SPM MEEG OBJECT

% Set filenames used in following steps
for subnum = 1:length(spm_files), % iterates over subjects
    spm_files{subnum}=[workingdir '/' spm_files_basenames{subnum}];
end

% load in the SPM MEEG object
subnum = 1;
D = spm_eeg_load(spm_files{subnum});

% look at the SPM object. Note that it is continuous data, with 232000 time
% points at 250Hz. We will epoch the data later.
D

%%%%%%%%%%%%%%%%%%%
%% DOWNSAMPLE
% (particularly important if movement compensation is used, as this stops 
% you from downsampling when running Maxfilter - but worth doing anyway)

% Set filenames used in following steps
for subnum = 1:length(spm_files), % iterates over subjects
    spm_files{subnum}=[workingdir '/' spm_files_basenames{subnum}];
end

S=[];
for subnum=1:length(spm_files), % iterates over subjects
    S.D=spm_files{subnum};
    S.fsample_new = 150; % in Hz
    D = spm_eeg_downsample (S);    
end
close all

%%%%%%%%%%%%%%%%%%%%
%% LOAD THE DOWNSAMPLED SPM MEEG OBJECT

% Set filenames used in following steps
for subnum = 1:length(spm_files), % iterates over subjects
    spm_files{subnum}=[workingdir '/d' spm_files_basenames{subnum}];
end

% load in the SPM MEEG object
subnum = 1;
D = spm_eeg_load(spm_files{subnum});

% look at the SPM object. Note that it is continuous data, with 139200 time
% points at 150Hz. We will epoch the data later.
D

%%%%%%%%%%%%%%%%%%%
%% OSLVIEW
% Note that there are some large artefacts. Use the oslview functionality
% to remove the bad epochs (see the osl_example_africa.m practical for how
% to do this). 

%set filenames used in following steps
for subnum = 1:length(spm_files), % iterates over subjects
    spm_files{subnum}=[workingdir '/d' spm_files_basenames{subnum}];
end

% first do a bit of high-pass filtering to remove the worst of the trends
S2=[];
S2.D=spm_files{1};
S2.band='high';
S2.freq=0.1;
D=spm_eeg_filter(S2);

% Now load oslview 
% This data has some bad artefacts in. Mark the epochs at around 325s,
% 380s and 600s as bad. Marking is done by right-clicking in the proximity
% of the event and click on 'Mark Event'. A first click (green dashed label)
% marks the begin of a bad period, a nother second click indicates the end 
% (in red). Also, mark the really bad artefacts at the end of the the
% experiment from about 650 secs to the end. This will mean that we are not
% using about half of the data. But with such bad artefacts this is the
% best we can do. We can still obtain good results with what remains.
D=oslview(D);


% AFRICA WITH MANUAL COMPONENT REJECTION

% % Run AFRICA ICA denoising. AFRICA uses independent component analysis to decompose sensor data into a set of maximally independent time courses. Using this framework, sources of interference such as eye-blinks, ECG artefacts and mains noise can be identified and removed from the data.
% % 
% % In this practical we will use manual artefact rejection by looking at the time courses and sensor topographies of each component and rejecting those that correlate with EOG and ECG measurements.
% % 
% % The user interface displays the time course, power spectrum and sensor topography for each component. These components are sorted based on one of a number of metrics, which you can toggle using the dropdown menu.

% Set new SPM MEEG object filenames to be used in following steps
for subnum = 1:length(spm_files), % iterates over subjects
    spm_files{subnum}=[workingdir '/fd' spm_files_basenames{subnum}];
end

% % AFRICA with manual artefact rejection: Scroll through components using the cursor keys. Identify the two components that correlate with the EOG and ECG measurements and mark them for rejection using the red cross. NB: Just close the window when finished to save your results.

%%% BE GENEROUS WITH YOUR COMP REJECTION

for subnum = 1:length(spm_files)

    [dirname,filename] = fileparts(spm_files{subnum});

    S = [];
    S.D           = spm_files{subnum};
    S.logfile         = 1;
    S.ica_file        = fullfile(dirname,[filename '_africa']);
    S.used_maxfilter  = 1;
    S.ident.func      = @identify_artefactual_components_manual;
    S.to_do           = [1 1 1];
    S.ident.artefact_chans = {'EOG','ECG'};

    osl_africa(S);

end














%%%%%%%%%%%%%%%%%%%%
% [NOTE: at this point you would normally do AFRICA denoising: but we do
% not do this as part of this practical]   

%% Run Coregistration
S = [];
S.D                 = D.fullfile;
S.mri               = [datadir '/structurals/struct1.nii'];
S.useheadshape      = 1;
S.use_rhino         = 0;
S.forward_meg       = 'Single Shell';
S.fid.label.nasion  = 'Nasion';
S.fid.label.lpa     = 'LPA';
S.fid.label.rpa     = 'RPA';
D = osl_headmodel(S);



% Check Coregistration
% RHINO HAS NOT BEEN RUN!
% rhino_display(D);

% EPOCHING
% 
% Does preliminary epoching for the purpose of finding outliers This is not the final epoching. Instead this sets up the epoch definitions, and performs a temporary epoching for the purpose of doing semi-automated outlier trial rejection (before running the fully automated OAT).
% 
% The epoch definitions and the continuous data will be kept and passed into OAT. This is so that things like temporal filtering (which is dones as part of OAT) can be done on the continuous data, before the data is epoched inside OAT.
% 
% Note that this will also remove those trials that overlap with the bad epochs identified using OSLview.
% 
% Here the epochs are set to be from -1000ms to +2000ms relative to the triggers in the MEG data. We also specify the trigger values for each of the 4 epoch types of interest (motorcycle images, neutral faces, fearful faces, happy faces).
% 
% The codes used on the trigger channel for this experiment were: 1 = Neutral face 2 = Happy face 3 = Fearful face 4 = Motorbike 18 = Introduction screen 11 = Break between blocks 19 = Midway break 12 = Green fixation cross (response trials) 13 = Red fixation cross (following green on response trials) 14 = Red fixation cross (non-response trials)
% 
% Note that we're only interested in the first 4 event codes listed here for today's workshop.


%set filenames used in following step
for subnum = 1:length(spm_files), % iterates over subjects
    spm_files{subnum}=[workingdir '/Afd' spm_files_basenames{subnum}];
end

for i=1:length(spm_files), % iterates over subjects

    %%%%
    % define the trials we want from the event information
    S2 = [];
    S2.D = spm_files{i};
    D_continuous=spm_eeg_load(S2.D);
    
    pretrig = -1000; % epoch start in ms
    posttrig = 2000; % epoch end in ms   
    S2.timewin = [pretrig posttrig];

    % event definitions
    S2.trialdef(1).conditionlabel = 'Neutral face';
    S2.trialdef(1).eventtype = 'STI101_down';
    S2.trialdef(1).eventvalue = 1;
    S2.trialdef(2).conditionlabel = 'Happy face';
    S2.trialdef(2).eventtype = 'STI101_down';
    S2.trialdef(2).eventvalue = 2;
    S2.trialdef(3).conditionlabel = 'Fearful face';
    S2.trialdef(3).eventtype = 'STI101_down';
    S2.trialdef(3).eventvalue = 3;
    S2.trialdef(4).conditionlabel = 'Motorbike';
    S2.trialdef(4).eventtype = 'STI101_down';
    S2.trialdef(4).eventvalue = 4;
    
    S2.reviewtrials = 0;
    S2.save = 0;
    S2.epochinfo.padding = 0;
    S2.event = D_continuous.events;
    S2.fsample = D_continuous.fsample;
    S2.timeonset = D_continuous.timeonset;
    
    [epochinfo.trl, epochinfo.conditionlabels, S3] = spm_eeg_definetrial(S2);        
    
    %%%%
    % do epoching
    S3=[];
    S3 = epochinfo;
    S3.D = D_continuous;     
    D = osl_epoch(S3);
            
end;



%%%%%%%%%%%%%%%%%%%%
%% LOAD THE EPOCHED SPM MEEG OBJECT

%set filenames used in following step
for subnum = 1:length(spm_files), % iterates over subjects
    spm_files{subnum}=[workingdir '/eAfd' spm_files_basenames{subnum}];
end

subnum = 1;             




D = spm_eeg_load(spm_files{subnum});

% look at the SPM object. Note that this is now EPOCHED data
D

% Display a list of trial types:
D.condlist

% Display time points (in seconds) per trial
D.time

% Identify trials of a certain type using the indtrial function. E.g.:

% SPM8 notation, got obsolete
% motorbike_trls = D.pickconditions('Motorbike')
motorbike_trls = indtrial(D,'Motorbike');


% get an overview about the available methods for the D object
methods(D)

%good_motorbike_trls = setdiff(motorbike_trls,D.badtrials);



% 
% Identify channels of certain types using the meegchannels function. E.g.
% identify the channel indices for the planar gradiometers and for the
% magnetometers (what you have available may depend on the actual MEG device)
% (Note that you can use 'MEGMAG' to get the gradiometers, and D.chantype gives you a list of all channel types by index).
% 
planars = D.indchantype('MEGPLANAR');
magnetos = D.indchantype('MEGMAG');


% We can access the actual MEG data using the syntax: D(channels, samples, trials). E.g. plot a figure showing all the trials for the motorbike condition in the 135th MEGPLANAR channel. Note that the squeeze function is needed to remove single dimensions for passing to the plot function, and D.time is used to return the time labels of the within trial time points in seconds.

figure; plot(D.time,squeeze(D(planars(135),:,good_motorbike_trls))); xlabel('time (seconds)');
% We can average over all the motorbike trials to get a rudimentary ERF (Event-Related Field):

figure; plot(D.time,squeeze(mean(D(planars(135),:,good_motorbike_trls),3))); xlabel('time (seconds)');

% Although we should bear in mind that this data is averaging over all data
% including noisy data segments, channels and trials. To do better than
% this we need to perform outlier rejection.







%%%%%%%%%%%%%%%%%%%
%% VISUAL ARTEFACT REJECTION 
% This runs a Fieldtrip interactive tool.
%
% - Pass over the first interactive figure as it is the EOG channel - so just
% press the "quit" button.
%
% This will bring up another interactive figure which will show the
% magnetometers. You can choose the metric to display - it is best to stick
% to the default, which is variance. This metric is then displayed for 
% the different trials (bottom left), the different channels (top
% right), and for the combination of the two (top left). You need to use
% this information to identify those trials and channels with high variance
% and remove them.  
%
% - Remove the worst channel (with highest variance) by drawing a box
% around it in the top right plot with the mouse. 
% - Now remove the trials with high variance by drawing a box
% around them in the bottom left plot.
% - Repeat this until you are happy that there are no more outliers.
%
% - Press "quit" and repeat the process for the gradiometers.

%set filenames used in following step
for subnum = 1:length(spm_files), % iterates over subjects
    spm_files{subnum}=[workingdir '/eAfd' spm_files_basenames{subnum}];
end

% RUN THE VISUAL ARTEFACT REJECTION:
for i=1:length(spm_files),
    S2=[];
    S2.D = spm_files{i};
    S2.time_range=[-0.2 0.4];
    D2=osl_rejectvisual(S2);
end;

% 
% 
% EXAMINE THE CLEANED EPOCHED DATA
% 
% We can now repeat the average over all the motorbike trials with the bad trials removed to get a rudimentary ERF (Event-Related Field).

% Set new SPM MEEG object filenames to be used in following steps
for subnum = 1:length(spm_files), % iterates over subjects
    spm_files{subnum}=[workingdir '/eAfd' spm_files_basenames{subnum}];
end

% load in SPM MEEG object
subnum = 1;
D = spm_eeg_load(spm_files{subnum});

% List the marked bad channels
D.badchannels

% List the marked bad trials
D.badtrials
% Identify the motorbike trials. Note that indtrial includes good AND bad trials, so bad trials need to be excluded.

motorbike_trls = setdiff(D.indtrial('Motorbike'),D.badtrials);

% Plot a cleaned rudimentary ERF
figure; subplot(1,2,1);plot(D.time,squeeze(mean(D(planars(135),:,motorbike_trls),3))); xlabel('time (seconds)');ylim([-8 8])
subplot(1,2,2);plot(D.time,squeeze(mean(D(magnetos(49),:,motorbike_trls),3))); xlabel('time (seconds)');ylim([-400 400])

% Plot a cleaned rudimentary ERF for all sensors
figure; 
subplot(1,2,1);imagesc(D.time,[],squeeze(mean(D([planars(:)],:,motorbike_trls),3))); xlabel('time (seconds)');
subplot(1,2,2);imagesc(D.time,[],squeeze(mean(D([magnetos(:)],:,motorbike_trls),3))); xlabel('time (seconds)');

% Plot a cleaned rudimentary ERF topography at N170
topo=squeeze(mean(D(:,188,motorbike_trls),3));

figure;
sensors_topoplot(D,topo,{'MEGPLANAR' 'MEGMAG'},1);

