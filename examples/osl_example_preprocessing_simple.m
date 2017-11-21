%% Preproc - Basic source space pipeline
%
% This example script shows how to perform basic preprocessing manually. 
% The main input to this pipeline is an MEEG object that has already
% had the coregistration and forward model run - refer to the
% relevant practicals for examples (in particular, the CTF data import tutorial). 
% This pipeline goes from a raw sensorspace coregistered MEEG to a source reconstruction. 

%% Overview
% Having read the raw files into an SPM MEEG object, we are now ready to implement 
% a basic source reconstruction pipeline. The main steps are
%
% # Update |D.inv| paths
% # Initial filtering (high pass, AC notch filter)
% # downsampling to reduce file size and processing time
% # Artefact detection to identify bad channels and epochs
% # ICA artefact removal
% # Optionally re-run forward model if changing the model type
% # Perform beamforming
%
% although you may not need or want to include all of these steps, depending on 
% the analysis you are running.

%% 
% First, we will load in the MEEG file. This example is a continuous resting
% state CTF recording
data_dir = fullfile(osldir,'example_data','preproc_simple_example');
spm_file = fullfile(data_dir,'3006.mat');
D = spm_eeg_load(spm_file);

%%
% Note that at any point, we can find an MEEG object on disk by checking the
% |fullfile| property, so you don't really need to keep track of the file name
% separately.
D.fullfile

%% Update D.inv paths
% One issue you might encounter is that the paths to structural files in |D.inv|
% may be different if coregistration was performed on someone else's machine, or even on the
% same machine but in a different directory. 
D.inv{1}.mesh

%%
% Notice how the structural files for this example are stored in the
% |preproc_simple_example| folder, but this is not the location contained in
% the MEEG object. This can be a problem if you want to display the
% coregistration e.g. with |rhino_display|. You can update these paths using
% |osl_update_inv_dir()|. The inputs to this file are the MEEG object, and the
% new path of the folder containing the structural files. This assumes that
% all of the files share a common root - this is generally a safe assumption
% if the files were created with RHINO.
D = osl_update_inv_dir(D,data_dir);
D.inv{1}.mesh

%% Setting chantypes and labels
% It's important that you set the channel types and labels correctly in your file, for two reasons
%
% * Bad samples are identified on a per-chantype basis. For example, you
%   detect artefacts separately in MEGMAG, MEGPLANAR, EMG and EOG. This works
%   best if the channel types are correct
% * Non-MEG channel types are recognized as such in AFRICA, as long as the chantype is not 'OTHER'
%
% Setting the labels and chantype correctly makes it easier for you and other users to work with the
% data further down the line, so it is always worth setting them correctly. In this case, suppose we know that
% the following channel labels correspond to artefact channels:
%
% * |EEG060| - EMG
% * |EEG059| - ECG
% * |EEG057| - EOG1
% * |EEG058| - EOG2
%
% At this point, we should make sure that these channel types are set correctly, and that the channel labels are informative
D = D.chantype(find(strcmp(D.chanlabels,'EEG060')),'EMG');
D = D.chantype(find(strcmp(D.chanlabels,'EEG059')),'ECG');
D = D.chantype(find(strcmp(D.chanlabels,'EEG057')),'EOG');
D = D.chantype(find(strcmp(D.chanlabels,'EEG058')),'EOG');
D = D.chanlabels(find(strcmp(D.chanlabels,'EEG060')),'EMG');
D = D.chanlabels(find(strcmp(D.chanlabels,'EEG059')),'ECG');
D = D.chanlabels(find(strcmp(D.chanlabels,'EEG057')),'EOG1');
D = D.chanlabels(find(strcmp(D.chanlabels,'EEG058')),'EOG2');

%% Initial filtering
% For filtering, we will use the |osl_filter| function. This function applies a basic Butterworth filter. You specify both 
% the filter type and key frequencies using two numbers - a lower and upper frequency |[f1,f2]|. The possibilities are
%
% * Low pass filter - |f1=0| e.g. |[0,45]| will be a low pass filter
% * High pass filter - |f2=inf| e.g. |[1,inf]| will be a high pass filter
% * Bandpass filter - specify the lower and upper frequencies e.g. |[8,13]| will filter to select the alpha band
% * Bandstop filter - specify the frequencies as negative to suppress them e.g. |-[48 52]| will be a notch filter that suppresses line artefacts at 50Hz
%
% We would typically apply a high pass filter to remove very slow drift, and then notch filters to remove AC line artefacts and their harmonics (if the sampling frequency of your data is high enough that they may be present)
D = osl_filter(D,[0.1 inf]); % Remove slow drift
D = osl_filter(D,-1*(50+[-2 2])); % Remove 50Hz with notch filter
D = osl_filter(D,-1*(2*100+[-2 2])); % Remove 100Hz with notch filter

%% Downsampling
% The raw files may have an unnecessarily high sampling rate for the analysis you wish to perform, which increases storage requirements and processing time.
% You can downsample your data at this point if you wish
D = spm_eeg_downsample(struct('D',D,'fsample_new',300)); 


%% Remove artefacts
% You can remove artefacts using |osl_detect_artefacts|. Broadly, there are
% two kinds of artefacts this function identifies and classifies
%
% * Bad channels - where an entire channel should be rejected. Rejection is
%   performed by setting |D.badchannels|
% * Bad times - periods of time in the recording that should be rejected. For
%   continuous recordings, this is performed by setting |D.badsegments|. For
%   epoched recordings, this is performed by setting |D.badtrials|.
%
% It can be very important to run |osl_detect_artefacts| or to otherwise perform artefact rejection 
% before doing ICA, as the presence of known artefacts can significantly degrade the ICA decomposition. 
D = osl_detect_artefacts(D);

%%
% By default, |osl_detect_artefacts| will identify both types of artefacts, although you can optionally
% specify which ones you want to detect. You can also choose different artefact detection metrics and thresholds, although
% the function is intended to use sensible defaults. See the artefact detection example for more detailed usage of |osl_detect_artefacts|. 
% function. 
%
% Aside from examining the properties of the MEEG object, it can also be helpful to examine the bad epochs visually. You can do this
% by opening the MEEG object in |oslview|
%
%   oslview(D)
%

%% AFRICA - ICA artefact removal
% Most of the information on how to perform ICA artefact removal is provided
% in the AFRICA example. Here, we are mainly concerned with how to perform
% automatic artefact rejection.  We set |artefact_channels| to correspond to
% the |chantypes| of all the channels whose correlations we want to check.
% |precompute_topos| is slow and generates large files, and is only mainly
% necessary if you are performing manual artefact rejection. So we can save
% time by skipping that step here. We set the |ident_func| to the automatic
% classification function. Lastly, if you examine
% |identify_artefactual_components_auto.m| you can see what options are
% available to control the rejection process. One common setting you might
% want to change is the correlation threshold for rejection. You might also
% want to enable or disable other types of ICA component rejection. Finally,
% call |osl_africa|
D = osl_africa(D,'artefact_channels',{'EOG1','ECG','EMG','EOG2'},'precompute_topos',false);
has_montage(D)

%%
% As you can see, a new online montage is created corresponding to the ICA-denoised signals.

%% Forward model
% You can check which forward model was used by looking at the contents of |D.inv|
D.inv{1}.forward.voltype

%%
% If you want to change the forward model, you can do that now
D = osl_forward_model(D,'forward_meg','Single Shell');
D.inv{1}.forward.voltype

%% Beamforming
% Finally, you can run the beamformer. You may want to use a filter at this point as well - 1-45Hz is typical. 
D = osl_filter(D,[1 45]); 
p = parcellation(8); 
D = osl_inverse_model(D,p.template_coordinates,struct('pca_order',150));

%%
% Notice that there are now online montages associated with the source reconstruction
has_montage(D)

%%
% For more information about beamforming, refer to the beamforming example on this website. 
% Now that the MEEG object has a beamformed montage, it is ready to be used in source space analysis
% e.g. source space connectivity analysis. 
