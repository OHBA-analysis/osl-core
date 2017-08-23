%% Data structure introduction
%
% There are a handful of useful data structures provided by SPM and OSL to facilitate working with
% MEG data. First up is the MEEG object. This is a container for the actual MEG data itself. These live on
% disk in two files
%
% * a |.dat| file, which contains the time series data
% * a |.mat| that contains header information
%
% You can load this object into Matlab using |spm_eeg_load|:
D = spm_eeg_load(fullfile(osldir,'example_data','roinets_example','subject_1'))

%%
% You can access the data by indexing this object
x = D(:,:,:);
size(x)

%%
% For continuous data, this matrix will be channels x time. For epoched data, there is an additional dimension corresponding to trial.
% You can access the time points easily as well, and use this to quickly plot your data
t = D.time;
plot(t,x(3,:))
xlabel('Time')
ylabel('Signal value')

%%
% How do you know what signal you are plotting? You can get information from the header about the channels
D.chanlabels(3) % channel name
D.chantype(3) % type of sensor
D.units(3) % data units

%%
% and there are other useful properties about the data in the MEEG object
D.fsample % Sampling rate
D.nchannels % Number of channels
D.fname % Name of the mat file

%%
% If you make changes to the MEEG object, these can be saved to disk using
D.save();

%%
% To copy the file, use
D2 = D.copy('./test');

%%
% This is one potentially confusing part of MEEG objects - most of their methods return new MEEG objects. So for example
D.fname
D2.fname

%%
% Now we will experiment with changing information in the header
clear D D2
D = spm_eeg_load('test.mat')

%%
% We have loaded the new copy so that we don't accidentally make changes to the original data. 
% Most properties of MEEG objects need to be set through their methods. For example, we can set the sampling rate using
D.fsample(200)

%%
% But there is a catch!
D.fsample

%%
% As you can see, |D.fsample(200)| returns a new MEEG object, without changing the old object. So we need to do
D = D.fsample(200);
D.fsample

%%
% This change is not automatically saved to disk. For example
D = spm_eeg_load('test.mat');
D.fsample

%%
% We need to save the updated header to disk manually
D = D.fsample(200);
D.save()
D = spm_eeg_load('test.mat');
D.fsample

%%
% This changes the contents of the '.mat' file but not the '.dat' file, because the raw data itself hasn't changed. 
% One extremely important thing to understand about MEEG objects is that they memory-map the .dat file. This means when you run
% |D(1,:,:)| you are indexing into the .dat file directly. While changing the header requires saving the MEEG object, any changes you 
% make to the raw data occur on disk immediately. So for example
D(1,1)
D(1,1) = 0

%%
% This immediately writes a value of 0 to disk. So if we reload the MEEG object:
D = spm_eeg_load('test.mat');
D(1,1)

%%
% You can see the change persists even though '.save()' was never called. This is the main reason why some OSL operations
% return a MEEG object without making changes on disk (because they change the header, which is only written when 'save()' is called)
% while others make copies of the MEEG object (e.g. filtering results in a new file with a prefix of 'f' because otherwise the original
% data would be instantly overwritten). It's important to keep this in mind, because if you happen to run
D(:,:,:) = 0

%%
% your data is now gone, and you would need to restore it from a separate backup (which you presumably made before starting!). In our case, we can just reload the 
% original file (before we copied it) - and make a new copy, just to be safe!
D = spm_eeg_load(fullfile(osldir,'example_data','roinets_example','subject_1'))
D.copy('./test');

%%
% The final important aspect of MEEG objects is the 'online montage'. An online montage is a linear combination of the original sensor data. 
% This can be used to represent any linear operation, including ICA artefact rejection, beamforming, parcellation, and leakage correction. 
% Writing these combinations as linear transformations rather than actual data enables a considerable saving of disk space. You can list
% the montages present using |has_montage|
has_montage(D)

%%
% To switch montages, you can use
D.nchannels
D = D.montage('switch',2);
D.nchannels

%%
% Montage 0 corresponds to the raw sensor data. You can also get information about the montage, including the linear transformation matrix
D.montage('getmontage',2)

%%
% To add an online montage of your own, use the |add_montage| function. For example, if we have a single channel that we want to correspond to 
% the sum of all of the sensors, we could use
D = D.montage('switch',0);
tmp = sum(D(:,1)) % The target value at the first timepoint
D = add_montage(D,ones(1,D.nchannels),'My new montage');
has_montage(D)
D(1,1) 

%%
% Note that because the montage information is stored in the header, the |save()| method must be used to save changes that affect montages e.g. 
D = spm_eeg_load(fullfile(osldir,'example_data','roinets_example','subject_1'))
has_montage(D)


%% NIFTI files


spatial_basis_file = fullfile(osldir,'parcellations','fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');

%%
% In general, an SPM MEEG object may have multiple online montages corresponding to 
% sensor space, source space, parcellated data, and orthogonalized data (refer to the
% Preprocessing tutorial for more information about online montages). When an MEEG
% object is passed to ROInets, the active montage must be in source space i.e. 
% with the same number of channels as there are voxels in the parcellation. 
% You can print a list of the montages stored in the MEEG object using the |has_montage()|
% function
has_montage(D)