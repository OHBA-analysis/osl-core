%% Import - CTF resting state
%
% This example shows how to read in raw data from CTF scanners

%%
% CTF raw data files are contained in a folder with extension |.ds|. This
% folder contains a number of binary files. For this example, we will use a
% recording from the UK MEG dataset gathered at the University of Nottingham
%
% _Please note that example data for this practical is not included in the
% public OSL release as it contains individual structural information_
%
ds_folder = fullfile(osldir,'example_Data','ctf_preprocessing','3006','3006_Eyes_Open_Rest_PROC.ds');

%% 
% In order to perform a source space analysis, we also need a structural MRI
% scan for the subject
mri_scan = fullfile(osldir,'example_Data','ctf_preprocessing','3006','3006_CRG.nii.gz');

%%
% To begin with, we can read the raw data in at the lowest level using a
% FieldTrip function
d=ft_read_data(ds_folder);
class(d)

%%
% Note that CTF data may optionally be processed using a synthetic higher-
% order gradiometer. Some information about this is available here:
% <http://www.fieldtriptoolbox.org/faq/how_does_the_ctf_higher-order_gradiometer_work>.
% This transformation may have already been applied to the data in the |.ds|
% folder. By default, FieldTrip will automatically undo this 'balancing'
% procedure when the data is imported into Matlab. This behaviour can be
% overriden by disabling this unbalancing in FieldTrip, by overwriting the
% relevant files. This is done automatically by |initialise_spm| which is
% called by |osl_startup|. If you want to use the original FieldTrip
% unbalancing, in |initialise_spm.m| at the top of the file set
% |disable_undobalancing = false|. You will need to restart Matlab and restart
% OSL for this to take effect. Matlab will display a warning whenever the
% custom unbalancing code is run, so that you are always aware when this
% occurs.

%%
% Without the header information, this data matrix is not very useful. We can
% read in both the header and the data using |ft_preprocessing|
%
% Secondly, the continuous data in the |.ds| folder is actually stored as a
% series of 10s epochs.
d=ft_preprocessing(struct('dataset',ds_folder))

%%
% Notice that even though we have a continous resting state recording, it has
% been stored as a sequence of trials. Each of these trials is 10s in
% duration, which is typical for CTF resting state recordings. If you know
% that you have continuous data you can automatically have these epochs
% stictched together
d=ft_preprocessing(struct('dataset',ds_folder,'continuous','yes'))

%%
% Although the above commands will read the data into Matlab, for OSL it is
% more useful to read the data into an SPM MEEG object. This can be performed
% using the |osl_import| function. The |.ds| file is passed in as the
% first input and you can optionally provide a filename for the MEEG object
% that will be created on disk. After conversion, the MEEG object is
% automatically loaded and returned
spm_file = fullfile(osldir,'example_Data','ctf_preprocessing','3006','3006');
D = osl_import(ds_folder,'outfile',spm_file)

%%
% By default, only data channels (e.g. MEG, trigger, clock) will be imported.
% You can optionally specify other channels to retain in the conversion. There
% are two possibilities
%
% * |artefact_channels| - these channels will be retained based on their
%   chantype e.g. |EMG|
% * |other_channels| - these channels will be retained based on their
%   chanlabel e.g. |EEG060|
%
% In this case, we have been told by our collaborators that the following
% additional channels were recorded
%
% * |EEG060| - EMG
% * |EEG059| - ECG
% * |EEG057| - EOG1
% * |EEG058| - EOG2
%
% First, we need to make sure these channels are retained in the import. These
% channels are identified by their label, so we will use |other_channels|
D = osl_import(ds_folder,'outfile',spm_file,'other_channels',{'EEG060','EEG059','EEG057','EEG058'})

%%
% Notice that we now have 280 channels instead of 276, because these 4
% channels have now been retained. For convenience and data integrity, we
% should now also take the opportunity to set the channel type in the MEEG
% object so that these channels can be readily identified further on in the
% analysis
D = D.chantype(find(strcmp(D.chanlabels,'EEG060')),'EMG');
D = D.chantype(find(strcmp(D.chanlabels,'EEG059')),'ECG');
D = D.chantype(find(strcmp(D.chanlabels,'EEG057')),'EOG');
D = D.chantype(find(strcmp(D.chanlabels,'EEG058')),'EOG');
D = D.chanlabels(find(strcmp(D.chanlabels,'EEG060')),'EMG');
D = D.chanlabels(find(strcmp(D.chanlabels,'EEG059')),'ECG');
D = D.chanlabels(find(strcmp(D.chanlabels,'EEG057')),'EOG1');
D = D.chanlabels(find(strcmp(D.chanlabels,'EEG058')),'EOG2');

%%
% In order to perform coregistration, in addition to the structural MRI scan,
% we also need headshape points. These are often stored in |.pos| files. You
% can read these in using the |ft_read_headshape()| function
pos_file = fullfile(osldir,'example_Data','ctf_preprocessing','3006','3006.pos');
ft_read_headshape(pos_file)

%%
% This struct can be passed straight into the MEEG's |fiducials| method to
% store it within the SPM object
D = D.fiducials(ft_read_headshape(pos_file,'unit','mm'));
D.fiducials

%%
% _Note that the headshape points MUST be read in units of 'mm', because this
% is assumed by RHINO later on in the pipeline_
%
% Don't forget to save your updated channel types and fiducials to disk!
D.save()

%% Coregistration
%
% Finally, we will need to perform the coregistration and forward model
% computation in order to be able to later do a source space analysis. To do
% this, we will use the |osl_headmodel| function. This function takes in the
% MEEG object, the structural scan, and some additional options including
% which forward model you would like to use (you can always run
% |osl_forward_model| again later if you want to change the forward model
% without changing the coregistration)
S = struct;
S.D = D.fullfile;
S.mri = mri_scan;
S.useheadshape = true;
S.forward_meg = 'MEG Local Spheres';
S.use_rhino = true;

%%
% You can refer to the coregistration example on this site for more
% information about RHINO and other coregistration options.
% 
% Note that the MEEG object *must* be passed in a filename. So you need to
% make sure that you have commited any pending changes to disk (by running
% |D.save()|) before using |osl_headmodel|. Lastly, you might need to rename
% your fiducials, depending on what their labels are in the |.pos| file. You
% can check what the labels are by looking at |D.fiducials|:
D.fiducials.fid.label

%%
% And then you should set these labels - notice for example that the 'nasion'
% point has been imported as 'nas'
S.fid.label.nasion='nas';
S.fid.label.lpa='lpa';
S.fid.label.rpa='rpa';

%%
% We can then go ahead and run |osl_headmodel| to do the coregistration and
% forward model.
osl_headmodel(S);

%%
% Don't forget to reload the MEEG object to read in the new data that was
% saved on disk.
D = spm_eeg_load(D.fullfile);

%% 
% You can verify that the coregistration was done by checking that |D.inv| now
% exists
D.inv

%%
% and you can verify that the forward model was run by checking that
% D.inv{1}.forward is not empty
D.inv{1}

%%
% And it's always worth verifying the quality of the coregistration
rhino_display(D)

%%
% Now you have a fully imported MEEG file that is ready for preprocessing and
% beamforming.


