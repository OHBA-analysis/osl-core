%% Preproc - Standard artefact rejection
%
% This example script shows how to perform basic artefact detection and rejection. 
% The main input required is an MEEG object.
%

%%
% A key part of most analysis pipelines is some form of artefact detection and rejection. Aside from ICA (and indeed before running ICA) it can be important
% to identify bad channels, epochs, or trials. This functionality is provided in OSL by |osl_detect_artefacts| - if you are already familiar with OSL, 
% this is essentially the same as |osl_detect_badevent| but a bit cleaner and more versatile. 
%%
% There are broadly two kinds of artefacts we might want to identify
%
% * Bad channels - where an entire channel should be rejected. Rejection is
%   performed by setting |D.badchannels|
% * Bad times - periods of time in the recording that should be rejected. For
%   continuous recordings, this is performed by setting |D.badsegments|. For
%   epoched recordings, this is performed by setting |D.badtrials|.
%
% To start with, let's load in an MEEG object. We will remove all existing bad channels
D = spm_eeg_load(fullfile(osldir,'example_data','roinets_example','subject_1'));
D = D.badchannels(1:D.nchannels,0);

%%
% First, we can run artefact detection for both bad channels and bad times
D2 = osl_detect_artefacts(D);

%%
% Under the hood, the continuous recording is epoched into dummy trials, outliers in those trials are used to identify bad segments, and then
% those bad segments are marked in the original recording (and the temporary epoched data is deleted). If you pass in an MEEG object that is already 
% epoched, then these epochs/trials will be used instead. This is determined by whether |D.type| is continuous or not. 
%
% |osl_detect_artefacts| has a number of different options. Firstly, you may only want to detect bad channels - for example, if you have 
% a set of trial data that hasn't been epoched yet, and want to just mark bad channels first (and reject bad trials later)
D2 = osl_detect_artefacts(D,'badtimes',false);

%%
% Alternatively, you may have already marked bad channels and now want to identify only bad times
D2 = osl_detect_artefacts(D,'badchannels',false);

%%
% Note that any bad segments, trials, or channels that are already marked will be carried forward in |osl_detect_artefacts| - that is, you
% will not lose anything that has already been marked.  
%
% Another important parameter you may wish to vary is the artefact rejection threshold, which controls the sensitivity of the artefact detection.
% The implementation is such that the sensitivity increases with the threshold - if you raise the threshold, you will get more artefacts. You 
% can set thresholds for bad channels (|channel_threshold|) and bad times (|event_threshold|) independently. To demonstrate the effect of the threshold, compare
% we can run the artefact detection for bad channels with different threshold values. 
D2 = osl_detect_artefacts(D,'badtimes',false,'channel_threshold',0.01); 

%%
% If we raise the threshold, we find more artefacts
D2 = osl_detect_artefacts(D,'badtimes',false,'channel_threshold',0.02);

%%
% And this trend continues
D2 = osl_detect_artefacts(D,'badtimes',false,'channel_threshold',0.05);

%%
% You can limit the maximum number of channels that are marked bad if you wish
D2 = osl_detect_artefacts(D,'badtimes',false,'channel_threshold',0.05,'max_bad_channels',2);

%%
% You can view the results of your artefact detection using |oslview|. Check |osl_detect_artefacts.m| for a full list of options. 