%% ROInets 2 - Amplitude envelope connectivity analysis
%
% This example shows how to examine functional connectivity using amplitude 
% envelope correlations for a signal subject.
%
% This example shows how to use low-level ROInets functionality to 
% compute parcel timecourses and perform orthogonalization using
% SPM objects
%
% First, set up the input file locations
spatial_basis_file = fullfile(osldir,'parcellations','fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
data_dir = fullfile(osldir,'example_data','roinets_example');
output_directory = fullfile(osldir,'practical','roinets_demo');
if ~exist(output_directory)
	mkdir(output_directory)
end

%%
% When ROInets is used with SPM objects, the node timecourses can be stored
% as online montages. However, this results in new montages being written to 
% the MEEG file on disk. Therefore, we make a copy first to prevent modifying 
% the originals
D = spm_eeg_load(fullfile(data_dir,'subject_1'));
D = D.copy(fullfile(output_directory,'subject_1'));
D = D.montage('switch',2);

%%
% The input MEEG file has 2 online montages, and we select the second one which
% corresponds to the source-reconstructed data. Note that there are 3559 channels,
% one for each voxel in the 8mm grid
D

%%
% We can load the spatial basis file to confirm the grid size and number of parcels.
% 23x27x23 corresponds to an 8mm grid, and there are 38 brain regions
spatial_basis = read_avw(spatial_basis_file);
size(spatial_basis) 

%%
% This can also be done using an OSL Parcellation object. The resolution is 8 and 
% n_parcels confirms there are 38 brain regions. n_voxels corresponds to the number
% of voxels in the template mask - here 3559 agrees with the number of channels in the
% MEEG object
p = parcellation(spatial_basis_file)

%% 
% Our first task is to compute a timecourse for each parcel. This is done using 
% |ROInets.get_node_tcs| e.g. |D = ROInets.get_node_tcs(D,spatial_basis_file,'pca')|.
% However, the PCA method only works if the parcellation is binary and non-overlapping. 
% Instead of specifying a .nii file, we can specify a matrix of voxel assignments
% n_voxels x n_parcels. We can assemble this matrix by first constructing a binary
% parcellation mask (23x27x23x38) and converting this to 3559x38
size(p.binarize) % Binarize the voxel assignments
size(p.to_matrix(p.binarize)) % Reshape from volume to matrix representation
D = ROInets.get_node_tcs(D,p.to_matrix(p.binarize),'pca');

%% 
% Now our MEEG object has 38 channels. If we tried to run get_node_tcs again, 
% an error would be thrown because the wrong montage is selected - to repeat, 
% first run |D = D.montage('switch',2)|
D

%%
% We can extract the timeseries from the D object and compute the parcelwise distribution of power.
% This can then be displayed as a spatial map, using fslview or using the Parcellation object
D = D.montage('switch',3);
ts = D(:,:,:);
ts = ft_preproc_bandpassfilter(ts, D.fsample, [8 12], 4, 'but','twopass','no');
parcel_power = sum(abs(ts),2)/size(ts,2)/(D.time(end)-D.time(1));
p.plot_activation(parcel_power);
% fname = p.savenii(p.to_vol(parcel_power),fullfile(output_directory,'parcel_power'));
% fslview(fname)

%%
% This could be compared to the original voxel data
D = D.montage('switch',2);
ts = D(:,:,:);
ts = ft_preproc_bandpassfilter(ts, D.fsample, [8 12], 4, 'but','twopass','no');
voxel_power = sum(abs(ts),2)/size(ts,2)/(D.time(end)-D.time(1));
p.plot_activation(voxel_power);

%%
% We could also plot seed-based power differences. For example, the first parcel is 
% in the left occipital cortex. We can plot the power difference between this region
% and all others
p.plot_activation(parcel_power-parcel_power(1),0.1);

%%
% Switching back to the parcel montage, we can compute the Hilbert envelope timecourses
D = D.montage('switch',3);
ts = D(:,:,:);
Hen = hilbenv(ts);
figure
plot(D.time,ts');
xlabel('Time (s)')
ylabel('Raw signal')
figure
plot(D.time,Hen');
xlabel('Time (s)')
ylabel('Amplitude envelope value')

%%
% And then can plot the amplitude envelope correlation
figure
imagesc(corr(Hen')+diag(nan(38,1)))
axis square
colorbar
title('Envelope correlation before leakage correction')

%%
% However, spatial leakage is still an issue
figure
imagesc(corr(ts')+diag(nan(38,1)))
axis square
colorbar
set(gca,'CLim',[-1 1])
title('Raw correlation before leakage correction')

%%
% We can correct for spatial leakage using |ROInets.remove_source_leakage|
% operating on the D file directly. Now, the active montage needs to be
% the one with the parcel timecourses
D = D.montage('switch',3);
D = ROInets.remove_source_leakage(D,'symmetric');
ts_lc = D(:,:,:);

%%
% We could also operate on the vector of data directly - this can be useful if
% your data are not being stored in an MEEG object. The same syntax can be used 
% for |ROInets.get_node_tcs| if your original data are not MEEG objects. 
ts_lc = ROInets.remove_source_leakage(ts,'symmetric');

%%
% The orthogonalized signals are now uncorrelated
figure
imagesc(corr(ts_lc')+diag(nan(38,1)))
axis square
colorbar
set(gca,'CLim',[-1 1])
title('Raw correlation after leakage correction')

%%
% But amplitude envelope correlations are still present
Hen_lc = hilbenv(ts_lc);
figure
imagesc(corr(Hen_lc')+diag(nan(38,1)))
axis square
colorbar
title('Envelope correlation after leakage correction')




