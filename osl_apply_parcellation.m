function [D,parcellation,assignments] = osl_apply_parcellation(S)
% Computes node time series from a source space MEEG object using a
% parcellation
%
% [D,assignments] = osl_apply_parcellation(S)
%
% INPUTS:
%
% S.D                   - source space MEEG object
% S.parcellation        - .nii or .nii.gz file at the same gridstep as S.D
% S.orthogonalisation   - Orthogonalisation protocol to apply:
%                           ['none','symmetric','closest','householder']
% S.method              - method for reconstructing parcel time course:
%                           ['PCA','mean','peakVoxel']
% S.prefix              - filename prefix for new MEEG object (default 'p')
%
% OUTPUTS:
%
% D             - new MEEG object containing the parcel time courses
% parcellation  - the parcellation weights for each voxel
% assignments   - parcel assignments for each voxel from a winner takes all
%                 voting

global OSLDIR 

% Check SPM File Specification:
try
    if isa(S.D,'meeg')
        D = S.D;
    else
        S.D = char(S.D);
        [pathstr,filestr] = fileparts(S.D);
        S.D = fullfile(pathstr,[filestr '.mat']); % force .mat suffix
        D = spm_eeg_load(S.D);
    end;
    D.check;
catch
    error('SPM file specification not recognised or incorrect');
end

D = montage(D,'switch',2);

S.prefix            = ft_getopt(S,'prefix','p');
S.orthogonalisation = ft_getopt(S,'orthogonalisation','none');
S.method            = ft_getopt(S,'method','PCA');

stdbrain = read_avw([OSLDIR '/std_masks/MNI152_T1_' num2str(getmasksize(D.nchannels)) 'mm_brain.nii.gz']);
try
    parcellation = vols2matrix(read_avw(S.parcellation),stdbrain); %nVoxels x nSamples
catch
    error('Make sure the parcellation file and the data are compatible, including having the same spatial resolution.');
end

if size(parcellation,2) == 1
    % Currently nvoxels x 1 with a index at each voxel indicating
    % parcel membership. Instead, needs to be nvoxels x nparcels binary labelling
    num_parcels = max(parcellation);
    parcellation_new = zeros(size(parcellation,1),num_parcels);
    for voxel = 1:size(parcellation,1)
        if parcellation(voxel) > 0
            parcellation_new(voxel,parcellation(voxel)) = 1;
        end
    end
    parcellation = parcellation_new;
    clear parcelflag_new;
end

% create parcellation based on winner takes all voting
% (this will be used for outputting results in nii format)
assignments = false(size(parcellation));
for voxel = 1:size(parcellation,1)
    [a,b] = max(parcellation(voxel,:));
    if a > 0
        assignments(voxel,b) = 1;
    end
end


% get voxel data
data = D(:,:,:);

voxeldata_concat = reshape(data,[D.nchannels,D.nsamples*D.ntrials]);
good_samples = ~all(badsamples(D,':',':',':'));
good_samples = reshape(good_samples,1,D.nsamples*D.ntrials);
voxeldata_concat = voxeldata_concat(:,good_samples);

nodedata = ROInets.get_corrected_node_tcs(voxeldata_concat, parcellation, S.orthogonalisation, S.method);

data = zeros(size(nodedata,1),length(good_samples));
data(:,good_samples) = nodedata;

data = reshape(data,[size(data,1),D.nsamples,D.ntrials]);

clear voxeldata_concat nodedata_concat;

% Save data to new MEEG object
outfile = prefix(fullfile(D.path,D.fname),S.prefix);
Dnode = clone(montage(D,'switch',0),outfile,[size(data,1),D.nsamples,D.ntrials]);
Dnode = chantype(Dnode,1:Dnode.nchannels,'VE');
Dnode(:,:,:) = data;
Dnode.save;

D = Dnode; % For output

end