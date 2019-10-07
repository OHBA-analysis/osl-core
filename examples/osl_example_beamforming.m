%% Beamforming in OSL
% This example shows how to perform beamforming

%%
% This practical will work with a single subject's data from an emotional
% faces experiment (Elekta Neuromag data).

%% SET UP ANALYSIS
% The only thing you need to do is to go into your OSL directory (i.e. type _cd /somedirectory/osl-core_ ) and then run the following.
osl_startup;

%%
% This will not be necessary after you have done this once, i.e. no need to
% repeat during one of the follow-up practicals.

% First, we will load in two example files, one continuous, and one epoched
% (to demonstrate that beamformer usage is the same for both types of data)

%D_continuous = spm_eeg_load(fullfile(osldir,'example_data','faces_singlesubject','spm_files','Aface_meg1'))

%%
%

D_epoched = spm_eeg_load(fullfile(osldir,'example_data','beamforming_example','ReBffdloc_S01'))

%%
% The MEEG object must be in a sensor-space montage - this could be the raw
% data (montage 0) or it could be an online montage obtained after running
% AFRICA
%D_continuous = D_continuous.montage('switch',0);

D_epoched = D_epoched.montage('switch',0);

%%
% It is worth noting here that each MEEG object keeps track of its filename
D_epoched.fullfile

%% 
% Some of the OSL functions result in the file being copied, and a prefix added. Typically this happens when filtering. You can update your
% |D| variable in place, and use the |fullfile| property to keep track of the updated filename instead of 
% dealing with the prefixes themselves. The only caveat is that you will end up with a lot of files on disk unless you delete them as you go

%% Filtering
% It is typical that we want the beamforming to focus on a particular frequency band.
% A frequency range of 1-45Hz is fairly typical for a standard ERF or induced response analysis, although you
% could beaform to more specific bands (e.g. alpha, beta). Here we will beamform the 1-45Hz band. You can do this
% using the |osl_filter| function, which we apply here to both the epoched
% and continuous data

D_epoched = osl_filter(D_epoched,[1 45]);
D_epoched.fullfile

%%
% Notice how the filename has automatically been updated. There is only a single variable (|D_epoched|) in Matlab, but there are two files on disk, one filtered and one unfiltered. Effectively, the unfiltered file is no longer loaded in Matlab. If you 
% didn't want to retain the original unfiltered data on disk, you could use
%
%   filtered_D_epoched = osl_filter(D_epoched);
%   D_epoched.delete()
%

%% Beamforming
%
% With MEGIN Neuromag data we have two sensor types, planar grads and magnetometers. 
% Before beamforming, if we have more than one sensor type, then we need
% to normalise the sensor types so that they can contribute equally to the
% beamformer weights. 

S=[];
S.D=D_epoched;
S.datatype='neuromag';
S.do_plots=true; 
[D_epoched pcadim] = osl_normalise_sensor_data(S);

%% 
%The main entry point for beamforming is |osl_inverse_model|. In order to
% perform beamforming, the MEEG object needs to have been coregistered and the
% forward model needs to have been run. You can determine if this is the case
% by examining the MEEG object's |inv| property
D_epoched.inv{1}

%%
% If the |forward| field is not empty, then you should be fine. Otherwise
%
% * If |forward| is empty but |D.inv| is not, then you need to run
%   |osl_forward_model|
% * If |D.inv| is empty or missing, then you need to run |osl_headmodel|.
% Note that osl_headmodel runs the coregistration AND osl_forward_model.
%
% The function that does the actual beamforming is called |osl_inverse_model|.
%
% This takes three arguments:
% 1) an MEEG object, 
% 2) an array of MNI coordinates to compute voxel timecourses at, and
% 3) an optional settings structure if you would like to override some of the defaults in
% |osl_inverse_model|. 
%
% You can get the set of MNI coordinates in two ways. If
% you are beamforming onto one of the OSL standard masks (i.e. something like
% 'MNI152_T1_8mm_brain.nii.gz') then you can create a parcellation object at
% that spatial resolution and use the template coordinates

p = parcellation(8);
mni_coords = p.template_coordinates;

%%
% This syntax might also be appealing if you plan to use a parcellation later
% in your analysis e.g.
%
%   p = parcellation('my_parcellation.nii.gz')
%   mni_coords = p.template_coordinates;
%
% Alternatively, you can get the set of MNI coordinates by passing a NIFTI file to |osl_mnimask2mnicoords|

mni_coords = osl_mnimask2mnicoords(fullfile(osldir,'std_masks','MNI152_T1_8mm_brain.nii.gz'));
size(mni_coords)

%%
% Either way, mni_coords should be an |n_voxels x 3| matrix. 
%
% Next, we will set up the setting structure. These fields are all optional
% but you may of course wish to override the defaults

S = struct;
S.modalities        = {'MEGPLANAR' 'MEGMAG'}; 
S.timespan          = [-0.3 0.5]; % in secs
S.pca_order         = 50;
S.type              = 'Scalar';
S.inverse_method    = 'beamform';
S.prefix            = '';
S.fuse              = 'meg';

%%
% One important option is |S.path| - the beamforming process creates a file 
% 'BF.mat' that is placed in a temporary directory. You can specify a particular
% directory if you like, but be aware that if you process multiple files in
% parallel, then you need to ensure that they each save BF.mat in a different
% location. By default, this will not be a problem because the temporary directory
% should be different for all processes.

%%
% The |S.timespan| option indicates the time range for which the data is
% used to estimate the source reconstruction weights. 

%%
% The |S.fuse='meg'| option means that the source reconstruction will fuse
% information from all MEG sensor types, in this case from both the MEG
% planar grads and the magnetometers.

%%
% The PCA order is a form of regularization and can help improve your results.
% The PCA order limits the rank of the output and thus sets an upper limit on
% the number of parcels you can compute orthogonal timecourses for e.g. if you
% set |pca_order=50| you will be unable to orthogonalize the parcel
% timecourses for a parcellation with 60 ROIs, and the PCA order must be less
% than or equal to the rank of the input data. For CTF data, a value of 150 is
% typical. For Elekta data that has been maxfiltered, a value of more like 50
% may be appropriate.
%
% Now we can actually run the beamformer on both the continuous and epoched
% data. Before running beamforming, we make sure that the normalised sensor space
% montage is the one being used. 

normalise_data_montage=1;
D_epoched = D_epoched.montage('switch',normalise_data_montage);

D_epoched = osl_inverse_model(D_epoched,mni_coords,S);

%%
% Note that |osl_inverse_model| writes changes to disk and thus will commit
% any unsaved changes you have made to the MEEG object to disk. This is
% something to be aware of if you are used to changing the montages in memory
% (e.g. by parcellation or orthogonalization) without saving to disk.
%

%%
% Note that the D_epoched (and correspondingly the D_continuous) has now
% got a number of channels equal to the number of MNI coordinates we
% specified in mni_coords:

D_epoched

disp('Number of channels:')
disp(D_epoched.nchannels)

%%
% You'll see that the result of running the beamformer is the addition of new online montages
% corresponding to the beamformed result

has_montage(D_epoched)

%%
% switch to the montage theat corresponds to the source recon

source_recon_montage=7;
D_epoched=D_epoched.montage('switch',source_recon_montage)

%% 
% We first identify 'RespRRespL' trials using the
% indtrial function

resp_trls = indtrial(D_epoched,'RespRRespL','good');

%%
% This gives all good trials for the RespRRespL condition in the data

%%
% We will focus on the activity at 0s, i.e. when the response was given. We
% expect a reduction of activity at that time on right motor cortex in
% response to moving the left hand. To do this we first find the time index for 0 secs

tpt=nearest(D_epoched.time,0); % find nearest time index to 0sec

%%
% We then compute the erf at 0s by averaging over all trials of this type,
% and find the voxel with the largest negative ERF amplitude at 0s

erf_zero_sec=squeeze(mean(D_epoched(:,tpt,resp_trls),3));
[tmp,voxel]=min(erf_zero_sec); % find dipole at which there is a maximum amplitude at 0secs


%%
% We can then plot the ERF the voxel with the largest negative ERF amplitude at 0s
figure;
plot(D_epoched.time,squeeze(mean(D_epoched(voxel,:,resp_trls),3))); 
xlabel('Time (seconds)');
ylabel('Magnetic field gradient (fT/mm)','FontSize',15);
set(gca,'FontSize',15)
title('Event-related field','FontSize',15)

%%
% We can then view the spatial map of the ERF at 0s.
% Remember, at 0sec we expect a negative response in right motor cortex.
% Make sure you turn on the negative color map in fsleyes, and threshold to
% see the largest negative value

fname_out=nii.quicksave(erf_zero_sec,'erf',8,8); % output spatial map of erf at 0sec
fsleyes(fname_out);
