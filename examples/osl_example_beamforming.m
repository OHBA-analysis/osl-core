%% Beamforming
% This example shows how to perform beamforming

%%
% This practical will work with a single subject's data from an emotional
% faces experiment (Elekta Neuromag data).
% 
% First, we will load in two example files, one continuous, and one epoched
% (to demonstrate that beamformer usage is the same for both types of data)

D_continuous = spm_eeg_load(fullfile(osldir,'example_data','faces_singlesubject','spm_files','Aface_meg1'))

%%
%
D_epoched = spm_eeg_load(fullfile(osldir,'example_data','faces_singlesubject','spm_files','eAface_meg1'))

%%
% The MEEG object must be in a sensor-space montage - this could be the raw
% data (montage 0) or it could be an online montage obtained after running
% AFRICA
D_continuous = D_continuous.montage('switch',0);
D_epoched = D_epoched.montage('switch',0);

%%
% It is worth noting here that each MEEG object keeps track of its filename
D_epoched.fullfile

%% 
% Some of the OSL functions result in the file being copied, and a prefix added. Typically this happens when filtering. You can update your
% |D| variable in place, and use the |fullfile| property to keep track of the updated filename instead of 
% dealing with the prefixes themselves. The only caveat is that you will end up with a lot of files on disk unless you delete them as you go

%% Filtering
% Before beamforming, you may want to filter your data. A frequency range of 1-45Hz is typical. You can do this
% using the |osl_filter| function:
D_continuous = osl_filter(D_continuous,[1 45]);
D_epoched = osl_filter(D_epoched,[1 45]);
D_epoched.fullfile

%%
% Notice how the filename has automatically been updated. There is only a single variable (|D_continuous|) in Matlab, but there are two files on disk, one filtered and one unfiltered. Effectively, the unfiltered file is no longer loaded in Matlab. If you 
% didn't want to retain the original unfiltered data on disk, you could use
%
%   filtered_D_continuous = osl_filter(D_continuous);
%   D_continuous.delete()
%

%% Beamforming
%
% The main entry point for beamforming is |osl_inverse_model|. In order to
% perform beamforming, the MEEG object needs to have been coregistered and the
% forward model needs to have been run. You can determine if this is the case
% by examining the MEEG object's |inv| property
D_epoched.inv{1}

%%
% If the |forward| field is not empty, then you should be fine. Otherwise
%
% * If |forward| is empty but |D.inv| is not, then you need to run
%   |osl_forward_model|
% * If |D.inv| is empty or missing, then you need to run |osl_headmodel|
%
% |osl_inverse_model| takes three arguments - an MEEG object, an array of MNI
% coordinates to compute voxel timecourses at, and an optional settings
% structure if you would like to override some of the defaults in
% |osl_inverse_model|. You can get the set of MNI coordinates in two ways. If
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
% Alternatively, you could pass a NIFTI file to |osl_mnimask2mnicoords|
mni_coords = osl_mnimask2mnicoords(fullfile(osldir,'std_masks','MNI152_T1_8mm_brain.nii.gz'));
size(mni_coords)

%%
% Either way, mni_coords should be an |n_voxels x 3| matrix. 
%
% Next, we will set up the setting structure. These fields are all optional
% but you may of course wish to override the defaults
S = struct;
S.modalities        = {'MEGPLANAR'}; 
S.timespan          = [0 Inf];
S.pca_order         = 250;
S.type              = 'Scalar';
S.inverse_method    = 'beamform';
S.prefix            = '';

%%
% One important option is |S.path| - the beamforming process creates a file 
% 'BF.mat' that is placed in a temporary directory. You can specify a particular
% directory if you like, but be aware that if you process multiple files in
% parallel, then you need to ensure that they each save BF.mat in a different
% location. By default, this will not be a problem because the temporary directory
% should be different for all processes.

%%
% The PCA order is a form of regularization and can help improve your results.
% The PCA order limits the rank of the output and thus sets an upper limit on
% the number of parcels you can compute orthogonal timecourses for e.g. if you
% set |pca_order=50| you will be unable to orthogonalize the parcel
% timecourses for a parcellation with 60 ROIs, and the PCA order must be less
% than or equal to the rank of the input data. For CTF data, a value of 150 is
% typical. For Elekta data that has been maxfiltered, a value of more like 60
% may be appropriate.
%
% Now we can actually run the beamformer
D_continuous = osl_inverse_model(D_continuous,mni_coords,S);
D_epoched = osl_inverse_model(D_epoched,mni_coords,S);

%%
% Note that |osl_inverse_model| writes changes to disk and thus will commit
% any unsaved changes you have made to the MEEG object to disk. This is
% something to be aware of if you are used to changing the montages in memory
% (e.g. by parcellation or orthogonalization) without saving to disk.
%
% The result of running the beamformer is the addition of new online montages
% corresponding to the beamformed result
has_montage(D_epoched)
