%% Beamforming in OSL
% This example shows how to perform beamforming

%%
% This practical will work with a single subject's data from an emotional
% faces experiment (MEGIN Neuromag data).

%% SET UP ANALYSIS
% The only thing you need to do is to go into your OSL directory (i.e. type _cd /somedirectory/osl-core_ ) and then run the following.
osl_startup;

%%
% This will not be necessary after you have done this once, i.e. no need to
% repeat during one of the follow-up practicals.

% First, we will load in an SPM MEEG object that has been preprocessed,
% including epoching. This is from a subject doing a simple finger tapping
% task.
%%
%

%D_epoched = spm_eeg_load(fullfile(osldir,'example_data','beamforming_example','ReBffdloc_S01'))
D_epoched = spm_eeg_load(fullfile(osldir,'example_data','beamforming_example','sept2019_vml_session1_epoched'))

%%
% The MEEG object must be in a sensor-space montage - this could be the raw
% data (montage 0), or it could be an online montage obtained after running
% AFRICA

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
% could beamform to more specific bands (e.g. alpha, beta). Here we will beamform the 1-45Hz band. You can do this
% using the |osl_filter| function, which we apply here to both the epoched
% and continuous data

D_epoched = osl_filter(D_epoched,[1 45]);
D_epoched.fullfile

%%
% Notice how the filename has automatically been updated. There is only a 
% single variable (|D_epoched|) in Matlab, but there are two files on disk, 
% one filtered and one unfiltered. Effectively, the unfiltered file is no 
% longer loaded in Matlab. If you 
% didn't want to retain the original unfiltered data on disk, you could use
%
%   filtered_D_epoched = osl_filter(D_epoched);
%   D_epoched.delete()
%

%% Beamforming
%
% With MEGIN Neuromag data we have two sensor types, planar grads and magnetometers. 
% Before beamforming, we need
% to normalise these two sensor types so that they can contribute equally to the
% beamformer calculation. In brief, this is done by scaling the different sensor types 
% to so that their variances over time are equal.

S=[];
S.D=D_epoched;
S.datatype='neuromag';
S.datatype='ctf';
S.do_plots=true; 
[D_epoched pcadim] = osl_normalise_sensor_data(S);

%% 
%The main entry point for beamforming is |osl_inverse_model|. In order to
% perform beamforming, the MEEG object needs to have been coregistered and the
% forward model needs to have been run. You can determine if this is the case
% by examining the MEEG object's |inv| property
D_epoched.inv{1}

%%
% If the |forward| field is not empty, then you should be fine. Otherwise:
%
% * If |forward| is empty but |D.inv| is not, then you need to run
%   |osl_forward_model|
% * If |D.inv| is empty or missing, then you need to run |osl_headmodel|.
% Note that osl_headmodel runs the coregistration AND osl_forward_model.

%% 
% The function that does the actual beamforming is called |osl_inverse_model|.
%
% This takes three arguments:
% 1) an MEEG object, 
% 2) a set of MNI coordinates to compute voxel timecourses at, and
% 3) an optional settings structure, used to override the default settings in
% |osl_inverse_model|. 
%
% You can get the set of MNI coordinates in two ways. If
% you are beamforming onto one of the OSL standard masks (i.e. something like
% 'MNI152_T1_8mm_brain.nii.gz') then you can create a parcellation object at
% that spatial resolution and use the template coordinates

spatial_res=10; % Spatial resolution of the voxel grid to beamform to, in mm
p = parcellation(spatial_res);
mni_coords = p.template_coordinates;

%%
% This syntax might also be appealing if you plan to use a parcellation later
% in your analysis e.g.
%
%   p = parcellation('my_parcellation.nii.gz')
%   mni_coords = p.template_coordinates;
%
% Alternatively, you can get the set of MNI coordinates by passing a NIFTI file to |osl_mnimask2mnicoords|

spatial_res=10; % Spatial resolution of voxel grid to beamform to, in mm
mni_coords = osl_mnimask2mnicoords(fullfile(osldir,'std_masks',['MNI152_T1_', num2str(spatial_res), 'mm_brain.nii.gz']));
size(mni_coords)

%%
% Either way, mni_coords should be an |n_voxels x 3| matrix. 
%
% Next, we will set up the setting structure. These fields are all optional
% but you may of course wish to override the defaults

S = struct;
S.timespan          = [-2 4]; % in secs
S.pca_order         = 55;
S.type              = 'Scalar';
S.inverse_method    = 'beamform';
S.prefix            = '';
S.modalities        = {'MEGPLANAR' 'MEGMAG'}; 
%S.fuse              = 'meg';
S.modalities        = {'MEG'}; 

%%
% The |S.timespan| option indicates the time window to be
% used to estimate the source reconstruction weights. Note however that the
% whole time range in D_epoched.time will still effectively be reconstructed.

%%
% The |S.fuse='meg'| option means that the source reconstruction will fuse
% information from all MEG sensor types, in this case from both the MEG
% planar grads and the magnetometers.

%%
% The PCA order is a form of regularization and can help improve your results.
% For CTF data, a value of 150 is
% typical. For MEGIN data that has been maxfiltered, a value 
%  less than ~60 will be appropriate (as this reflects the fact that after default
% maxfiltering, the rank of the MEG data in sensor space is ~64).
%
% Now we can actually run the beamformer on both 
% data. Before calling |osl_inverse_model|, we make sure that the normalised sensor space
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
% You'll see that the result of running the beamformer is the addition of two new online montages
% corresponding to the beamformed result:
% 1) Source space (MEG) without weights normalisation
% 2) Source space (MEG) with weights normalisation
% Weights normalisation is used to correct the fact that, with beamforming, 
% voxels that are deeper in the brain tend to have higher variance. 

has_montage(D_epoched)

%%
% Switch to the montage that corresponds to the source recon
% with weights normalisation
source_recon_montage=3;
D_epoched=D_epoched.montage('switch',source_recon_montage)

%% 
% To see the conditions (i.e. trial types) in the data
% we can use:

D_epoched.condlist

%%
% Here we want to focus on 'RespRRespL' trials, which we can do using the
% indtrial function:

resp_trls = indtrial(D_epoched,'StimLRespL','good');
resp_trls = indtrial(D_epoched,'Stim_Onset','good');

%%
% This gives all good trials for this condition in the data.

%%
% We then compute the erf by averaging over all these trials. Note that
% because of the ambiguity of the sign of the activity following source
% reconstruction (see
% https://ohba-analysis.github.io/osl-docs/pages/docs/sign_ambiguity.html),
% we work with the absolute value of the ERF:

erf=squeeze(mean(D_epoched(:,:,resp_trls),3));
abs_erf=abs(erf);

%%
% We can then save out and view the spatio-temporal activity of the
% abs(ERF), as a 4D niftii file where the 4th dimension is time.

fname_out=nii.quicksave(abs_erf,'abs_erf',spatial_res,spatial_res);
fsleyes(fname_out);

%%
% Once FSLeyes is open, make sure you:
% * select the abs_erf image in the overlay list
% * turn on the negative color map in fsleyes
% * set the min to 5
% * select "View->Time-series"
% * find the evoked activity just after 0secs, NOTE: fsleyes currently
% shows the time-axis (i.e. the x-axis) using the time index (i.e. the volume-index in the
% 4D niftii file being shown). To know what
% this is in seconds, you need to look at the values in |D_epoched.time|. To find the
% index for 0 secs, you can use:

disp(nearest(D_epoched.time,0))

%% Time-Frequency event-related (de)synchronisations
% As well as viewing the ERF, we can also look at oscillatory power using
% the time-frequency TF transform (i.e. the induced response, which will
% reveal event-related synchronisations and de-synchronisations)

%%
% We start by cropping the data in time to the window in time we are
% interested in (this will speed things up too). 
% Note that we need to temporarily switch to no montage being applied to
% be able to use |spm_eeg_crop|

D_epoched=D_epoched.montage('switch',0)

S = struct;
S.D = D_epoched;
S.timewin = [-3 4]*1000; % in msecs
D_epoched_timecropped = spm_eeg_crop(S);

D_epoched_timecropped = D_epoched_timecropped.montage('switch',source_recon_montage);

%%
% Here we want to focus on 'Abduction' trials, which we can do using the
% indtrial function:

resp_trls = indtrial(D_epoched_timecropped,'Abduction','good');

%%
% We then use |osl_tf_transform| to do the time-frequency TF transform.
% Here we are
% using a Hilbert transform within the beta-band frequency range, i.e. 13-30Hz.

S = struct;
S.tf_method = 'hilbert';
S.tf_freq_range = [13,30];
S.tf_num_freqs = 1;
S.raw_times = D_epoched_timecropped.time;
S.ds_factor = 0.5; % smaller value means less time samples in result
dat = osl_tf_transform( S , D_epoched_timecropped(:,:,resp_trls) );

%%
% We can now compute the induced response in the beta band by averaging
% over trials:
induced_response_beta = mean(dat.dattf(:,:,:,:),3);

%%
% Basline correction, carried out on trial averaged data, separately for
% each frequency bin

S=struct;
S.data = induced_response_beta; % pass in trial averaged beta power data [nchannels x ntpts]
S.time = dat.tf_times; % vector or tpts for 2nd dimension of S.data
S.time_window = [-Inf -1.5]; % [start end] in secs
induced_response_beta_bc = osl_baseline_correct(S);

%%
% Save out and open fsleyes

fname_out=nii.quicksave(induced_response_beta_bc,'induced_response_beta_bc',spatial_res,spatial_res); % output spatial map of erf at 0sec
fsleyes(fname_out);

%%
% Using fsleyes, find a voxel with high value in the post movement beta
% rebound at about 1-2 sec post stimulus. To find the time index for 1.3 sec
% post stimulus, we can use:

nearest(dat.tf_times,1.3)

%%
% Find a voxel with high beta rebound (a high positive value, indicating a
% high beta power or beta ERS, and enter it here in MNI coordinates in mm:

beta_ers_mnicoord=[-30 -21 58]; % in mm

%%
% Find voxel index nearest to that:
beta_ers_voxel_index=nearest_vec(mni_coords,beta_ers_mnicoord);

%%
% We can use this next bit of code to find the time and voxel with maximum activity in the beta
% band

%%
% Plot time course of beta power at peak voxel:
figure;
plot(dat.tf_times,induced_response_beta(beta_ers_voxel_index,:));
xlabel('time (s)');
ylabel('beta power');

%% 
% We can also plot the TF transform across a range of freq bands at this peak ERS voxel

S = struct;
S.tf_method = 'morlet';
S.tf_freq_range = [1,40];
S.tf_num_freqs = 30;
S.raw_times = D_epoched_timecropped.time;
S.tf_morlet_factor=5;
S.ds_factor = 0.5;
dat = osl_tf_transform( S , D_epoched_timecropped(beta_ers_voxel_index,:,resp_trls) );

%%
% Basline correction, carried out on trial averaged data, separately for
% each frequency bin

S=struct;
S.data = mean(dat.dattf,3); % pass in trial averaged TF data [nchannels x ntpts x 1 x nfreq]
S.time = dat.tf_times; % vector or tpts for 2nd dimension of S.data
S.time_window = [-Inf -1.5]; % [start end] in secs
dat_bc = osl_baseline_correct(S);

%%
% Plot TF at peak ERS voxel

figure;
tf = squeeze(dat_bc(:,:,:,:))';
grid on;
contourf(dattf.tf_times,dattf.tf_freqs,tf,32,'linestyle','none' )
colorbar
xlabel('Time (seconds)');
ylabel('Power','FontSize',15);
set(gca,'FontSize',15)
title('Event-related field','FontSize',15)

