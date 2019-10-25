%% Preproc - OPT
% This is an example for running the OHBA recommended preprocessing
% pipeline on Elekta-Neuromag MEG data (a very similar pipeline will work on
% CTF and EEG data as well) using OPT (OSL's preproscessing tool). It works through
% similar steps as you would for the manual preprocessing, but this time it is all automated.
%

%%
% We will work with two subjects' data from a button press
% experiment. The data should be available in your installation already.
% Note that this contains the fif files, e.g. fifs/loc_S01_sss1.fif, which
% have
% already been SSS Maxfiltered and downsampled to 250 Hz, and which we will
% use for this analysis.

%%
% *SPECIFY DIRECTORIES FOR THIS ANALYSIS*
% 
% First, we set the data directory where the
% data is, this should be the correct path already:
datadir = fullfile(osldir,'example_data','preproc_example','automatic_opt')
 
%%
% We now need to specify a list of the fif files that we want to work on. 
% We also need to specify a list of SPM file names, 
% which will be used to name the resulting SPM data files after the fif files have been 
% imported into SPM. It is important to make
% sure that the order of these lists is consistent. 

%
% fif_files{1}=[testdir '/fifs/sub1_face_sss.fif'];
% fif_files{2}=[testdir '/fifs/sub2_face_sss.fif']; 
% etc...

clear fif_files spm_files;

fif_files{1}=fullfile(datadir,'fifs','loc_S01_sss1.fif');
fif_files{2}=fullfile(datadir,'fifs','loc_S02_sss1.fif');

spm_files{1}=fullfile(datadir,'spm_files','loc_S01');
spm_files{2}=fullfile(datadir,'spm_files','loc_S02');

%% CONVERT FROM FIF TO AN SPM M/EEG OBJECT
% OPT takes in SPM MEEG objects, so we first need to convert the fif file
% for each subject into this format using |osl_import|.
%    
% The fif files that we are working with, e.g. loc_S02_sss1.fif, have
% already been max-filtered and downsampled to 250Hz.
%
% The SPM M/EEG object is the data structure used to store and manipulate
% MEG and EEG data in SPM.
% 
% Calling |osl_import| will also produce a histogram plot showing the number of events detected
% for each code on the trigger channel. 

if(~isempty(fif_files))
    S2=[];
    for i=1:length(fif_files) % loops over subjects
        S2.outfile = spm_files{i};       
        S2.trigger_channel_mask = '0000000000111111'; % binary mask to use on the trigger channel
        D = osl_import(fif_files{i},S2);

        % The conversion to SPM will read events and assign them values depending on the 
        % trigger channel mask. Use |report.events()| to plot a histogram showing the event
        % codes to verify that the events have been correctly read. 
        report.events(D);
    end
end


%%
% *SET UP STRUCTURALS*
% 
% OPT can also perform coregistration (estimating the spatial position of a 
% subject's head and brain relative to the MEG sensors). 
% Since coregistration requires a structural MRI for 
% each subject, we need to supply these to OPT. Coregistration and its role 
% in source reconstruction will be explained in a different practical. 
%
% Here, we setup a list of existing structural files, in the same
% order as in spm_files and fif_files.

clear structural_files;

structural_files{1}=fullfile(datadir,'structs','anat1.nii'); % leave empty if no structural available
structural_files{2}=fullfile(datadir,'structs','anat2.nii'); % leave empty if no structural available

%% SETTING UP AN OPT ANALYSIS
% This sets up an OPT struct to pass to |osl_check_opt|, by setting the
% appropriate fields and values in the OPT struct. Note that some fields
% are mandatory while others are optional (and will be automatically set to
% their default values). The osl_check_opt.m function should be used to
% setup the settings for OPT. This function will check the settings, and
% will throw an error if any required inputs are missing, and will fill
% other settings that are not passed in with their default values. The OPT
% structure can then be passed to osl_run_opt to do an OPT analysis. On the
% Matlab command line type |help osl_check_opt| to see what the mandatory
% fields are. 
% Note that you MUST specify:
%
% |opt.spm_files|: A list of the spm meeg files for input into SPM
%
% AND:
%
% |opt.datatype|: Specifies the datatype, i.e. 'neuromag', 'ctf', 'eeg'
% 
% For more information, see
%
% <html>
% <a href="https://ohba-analysis.github.io/osl-docs/pages/docs/opt.html" target="_blank">https://ohba-analysis.github.io/osl-docs/pages/docs/opt.html</a>
% </html>
 
%%
% *Specify required inputs*
% 
% We need to specify the list of input SPM MEEG object files and data type: In our case the input 
% files were already setup above in the variable |spm_files|, and we are 
% using data acquired by the MEGIN Neuromag system
% (same type as in manual preproc practical).
 
opt=[];
opt.spm_files=spm_files;
opt.datatype='neuromag';

%%
% *Specify opt directory name*
% 
% This is the name of the directory (full path) where the OPT results will be stored,
% and is given a ".opt" extension. Note that each OPT directory is
% associated with an OPT run - if you rerun OPT with the same |opt.dirname| then
% this will overwrite an old directory, and the old OPT results will be
% lost. Hence, you should ensure that you change |opt.dirname| for a new
% analysis, if you want to avoid overwriting an old one!

opt.dirname=fullfile(datadir,'practical_singlesubject.opt');

%% 
% *Highpass filtering and mains/line noise filtering settings*
%
% Here, we set both the highpass filter and mains/line noise filters to attenuate
% slow drifts and 50 Hz line noise. This corresponds to our filtering part
% during the manual preprocessing; now OPT takes care of it.

opt.highpass.do=true;
opt.highpass.cutoff=0.1; %Hz

% Mains/line noise notch filter settings, note that this will remove 50Hz plus harmonics
opt.mains.do=true;

%% 
% *Downsampling settings*
% 
% This is where we specify any downsampling to the desired sampling frequency.
% Here, we will downsample to 150Hz.

opt.downsample.do=true;
opt.downsample.freq=150; %Hz


%%
% *Identifying bad segements settings*
%
% This identifies bad segments in the continuous
% data (similar to using oslview in the manual practical, just automated).

opt.bad_segments.do=true;

%%
% *AFRICA denoising settings*
%
% Automatic AFRICA (ICA denoising) is available in OPT, but it is not 
% recommended, as the automatic labelling of artefacts is only a beta release.
% As a result, we turn this off (it will be turned off by default)
opt.africa.do=false;

%%
% *Epoching settings*
%
% Here the epochs are set to be from -1s to +2s relative to the stimulus onset in the
% MEG data.
opt.epoch.do=true;
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

%% 
% *Outlier trial detection settings*
%
% Instead of identifying bad segments in the continuous data, we will rely
% on opt to identify bad trials in the epoched data using the opt.outliers
% settings. This is roughly equivalent to using osl_reject_visual during
% the manual procedure.

opt.outliers.do=true;

%%
% *Coregistration settings*
%
% We're not doing coregistration here, but normally
% you would if you want to do subsequent analyses in source space. This
% requires structural scans.
opt.coreg.do=false; 
opt.coreg.mri=structural_files;

%% CHECK OPT SETTINGS
% Checking chosen settings: By calling |osl_check_opt| we will check the
% validity of the OPT parameters we have specified. Then, OPT will fill in any missing
% parameters with their default values for us.

opt=osl_check_opt(opt);

%% LOOK AT OPT SETTINGS AND SUBSETTINGS
%%
% *DISPLAY OPT SETTINGS*
%
% This gives an overview of the set parameters

disp('opt settings:');
printstruct(opt)

%%
% *LOOK AT OPT SUB-SETTINGS*
% 
% The OPT structure contains a number of subfields containing the settings
% for the relevant stages of the pipeline. Note that each of these has a
% "do" flag (e.g. |opt.downsample.do|), which indicates whether that part of
% the pipeline should be run or not. 

%%
disp('opt.downsample settings:');
disp(opt.downsample);

disp('opt.highpass settings:');
disp(opt.highpass);

disp('opt.mains settings:');
disp(opt.mains);

disp('opt.epoch settings:');
disp(opt.epoch);

disp('opt.outliers settings:');
disp(opt.outliers);

disp('opt.coreg settings:');
disp(opt.coreg);

%%
% Note the |opt.cleanup_files| flag. This is used to indicate whether SPM files generated by each 
% stage of opt should be cleaned up as the pipeline progresses:
%
% * 0 indicates nothing will be deleted
%
% * 1 indicates all files will be deleted apart from
% pre/post AFRICA files 
%
% * 2 indicates that all files will be deleted.
%
% By default this flag is set to 1, so that the files are cleaned up as OPT
% progresses.

opt.cleanup_files
    
%% RUNNING THE OPT ANALYSIS 
% We will now run the main OPT analysis:
opt=osl_run_opt(opt);

%% VIEWING OPT RESULTS 
% 
% There is several ways to look how OPT has been run:
%%
% *VIEWING OPT RESULTS IN MATLAB*
% 
% Running the OPT analysis will create an OPT output directory (whose name
% is the name set in opt.dirname with a ".opt" suffix added). This contains
% all you need to access the results of the analysis. Note that you can
% load these into Matlab using the call:

opt = osl_load_opt(opt.dirname);

disp('opt.results:');
disp(opt.results);

%%
% In particular, the OPT object contains a sub-struct named results, (i.e.
% |opt.results|), containing:
%
% * .logfile (a file containing the matlab output) 
% * .report (a file corresponding to a web page report with diagnostic plots) 
% * .spm_files (a list of SPM MEEG object files corresponding to the continuous data (before epoching), e.g. to pass into an OAT analysis) 
% * .spm_files_epoched (a list of SPM MEEG object files corresponding to the epoched data, e.g. to pass into an OAT analysis)
% 
% For example, the SPM MEEG object corresponding to the now preprocessed
% epoched data ouput by OPT for subject 1 is:

subnum=1;
D=spm_eeg_load(opt.results.spm_files_epoched{1})

%%
% It is recommended that you always inspect the 
% |opt.results.report|, to ensure that OPT has run
% in a sensible manner. We will look at this next.


%%
% *VIEWING OPT RESULTS BY CHECKING OPT REPORTS IN BROWSER*
%
% Open the web page report file indicated in |opt.results.report| ) in a web
% browser:

opt.results.report.html_fname

%%
% Note, that a link is also provided in the matlab output of |osl_run_opt|.

%%
% This displays diagnostic plots. At the top of the file
% is a link to |opt.results.logfile| (a file containing the matlab output) -
% you can check this for any errors or unusual warnings. Then there will be a list
% of session specific reports. Here, we have preprocessed data for two
% sessions, and so there is a link available for each of these.
%
% Below that, there are some summary plots showing the number of:
%
% * bad trials
% * bad channels 
% * bad segments 
%
% detected for each subject/session.

%%
% Click on the link for Session 1.
% 
% At the top of the webpage, you can see the names of the input and output (final preprocessed) SPM
% MEEG objects files is listed.
%
% Below that, there are a number of diagnostic plots for this session. In order:
%
% *opt-bad_segments: Bad Channels*
% 
% Shows histograms and scatterplots of the standard deviations of 
% channel data with bad channels indicated.
% The scatterplots show the channel number versus the metric (e.g.
% "std" for standard deviation) as red crosses before rejection and green crosses after rejection.
% Channels to be retained are indicated by green circles.
% These are shown for both sensor types.
%
% *opt-bad_segments: Bad segment timings*
% 
% Shows you a plot of bad segments found in continuous data (prior to epoching).
% Bad segments are shown in black.
%
% *opt-epoch: Event timings*
%
% Visualisation of the trial timings in continuous data (prior to epoching). 
% Shows the different trial types alongside the
% samples that have been marked as bad.
% The start of a trial is marked with an 'o', and the end is marked with an
% 'x'. Bad trials are shown in black.
%
% *opt-outliers: Bad Channels*
%
% Shows histograms and scatterplots of the standard deviations of 
% channel data with bad channels indicated, before and after outlier
% trial detection.
%
% *opt-outliers: Bad Trials*
%
% Shows histograms and scatterplots of the standard deviations of 
% trial data with bad trials indicated.
% Scatterplots show the channels/trial number versus the metric (e.g.
% "std") as red crosses before rejection and green crosses after rejection.
% Channels/trials to be retained are indicated by green circles.


%% 
% *CHECKING OPT RESULTS BY LOOKING AT THE DATA*
%
% Last but not least, you might want to look at your actual data to check
% whether OPT gives your good results: We will now load the M/EEG object created by OPT (analogous to our
% resulting |D| objects in the manual preproc practical) for subject 2.

subnum=2;
D=spm_eeg_load(opt.results.spm_files_epoched{subnum});

%%
% We can then define some trials to look at:
% good_stimresp_trls = [D.indtrial('StimLRespL','good') D.indtrial('StimLRespR','good')];
allconds=D.condlist;
good_stimresp_trls = [D.indtrial(allconds(5:8),'good')]; % takes button press conditions

%%
% Next, we get the sensor indices for the two different MEG acquisition
% modalities from the data:
planars = D.indchantype('MEGPLANAR');
magnetos = D.indchantype('MEGMAG');

%%
% Finally, as in the manual preprocessing practical, we are going to have a
% quick look at data quality by just doing some preliminary and rudimentary ERF
% analysis. We will use the loaded |D| object, all good stimulus response
% trials and average them to get an idea about the data quality after OPT.

figure('units','normalized','outerposition',[0 0 0.6 0.3]); 
subplot(1,3,1); % plots gradiometer ERF image
imagesc(D.time,[],squeeze(mean(D([planars(:)],:,good_stimresp_trls),3)));
xlabel('Time (seconds)','FontSize',20);
ylabel('Sensors','FontSize',20);colorbar
title('ERF, planar gradiometers','FontSize',20)
set(gca,'FontSize',20)
set(gca,'XLim',[-1 1])

subplot(1,3,2);  % plots magnetometer ERF image
imagesc(D.time,[],squeeze(mean(D([magnetos(:)],:,good_stimresp_trls),3))); 
xlabel('Time (seconds)','FontSize',20);
ylabel('Sensors','FontSize',20);colorbar
title('ERF, magnetometers','FontSize',20)
set(gca,'FontSize',20)
set(gca,'XLim',[-1 1])

subplot(1,3,3); % plots 1 chosen planar gradiometer (over motor cortex) time-course
plot(D.time,squeeze(mean(D(planars(135),:,good_stimresp_trls),3)));
xlabel('Time (seconds)','FontSize',20);ylim([-15 10])
set(gca,'FontSize',20)
ylabel(D.units(planars(1)),'FontSize',20);
set(gca,'XLim',[-1 1])
title('ERF at sensor 135','FontSize',20)

%%
% Note that the ERFs have been epoched from -1 to 2secs with 0secs corresponding to the onset
% of a hand movement. 
%
% These ERFs should look reasonable, i.e. both the ERF across sensors as
% well as the single-sensor ERF (taken from a sensor over the motor cortex)
% should look sufficiently smooth, you should
% basically see a candidate 'Bereitschaftspotential' close to 0secs!

%%
% Together, the three above checks should give you a sufficiently good idea about
% your data quality. As a rule of thumb, always check your data, especially
% after running long chains of automated analyses like OPT. Once you are in
% source-space it will be even harder to tell whether your data has
% sufficient data quality or is contaminated by artefacts.

%% EXERCISES
% 
% Now that you have seen the wonders of automated preprocessing, why not
% take a look at the really bad data from the manual preprocessing
% practical? Open the corresponding script to identify its location and try
% to adapt the OPT template script described here to run the problematic data set via OPT. Keep in mind that
% this data was exceptionally bad, so expect to have to test and tweak your settings until you
% reach a satisfying output (if at all).
