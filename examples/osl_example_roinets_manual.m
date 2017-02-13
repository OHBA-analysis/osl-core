%% Example script for ROInets-based analysis
% This script shows how to run ROInets manually
%
osldir = getenv('OSLDIR');
data_dir = fullfile(osldir,'example_data','roinets_example');
subjects = 1:10;

%%
% Copy the demo files

D_files = {};
mkdir(fullfile(osldir,'practical','roinets_demo'));
for j = 1:length(subjects)
	fname = sprintf('subject_%d',j);
	D = spm_eeg_load(fullfile(data_dir,fname));
	D_files{j} = D.copy(fullfile('roinets_demo',fname));
end

%%
% Load a parcellation
p = nii_quickread(fullfile(osldir,'parcellations','fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'),8);

%%
% Each voxel can belong to at most one parcel
% Assign each voxel to a single parcel
[~,assignments] = max(p,[],2);
assignments(all(p==0,2)) = 0;

%%
%   PARCELFLAG is a logical array (nVoxels x nParcels) which identifies
%   the voxels in each parcel. Each column identifies the voxels 
%   making up each parcel. True entries indicate membership. 
parcelflag = zeros(size(p));
for j = 1:size(p,2)
	parcelflag(:,j) = assignments == j;
end

%%
% Need to make sure montages are correct
D_files{1}
D_files{j}.montage('switch',2)

%%
% Make PCA timecourse montages
for j = 1:length(subjects)
	D_files{j} = D_files{j}.montage('switch',2); % Make sure montage is correct!
	D_files{j} = ROInets.get_node_tcs(D_files{j},parcelflag,'PCA',[],true);
end

%%
% New timecourses are in a montage
D_files{1}


