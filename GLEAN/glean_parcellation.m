function D = glean_parcellation(S)
% Computes node time series from a MEEG object using a parcellation.
%
% D = GLEAN_PARCELLATION(S)
%
% REQUIRED INPUTS:
%
%   S.D                 - SPM12 MEEG object
%   S.parcellation      - .nii or .nii.gz file at the same gridstep as S.D
%                         OR a .mat file containing a matrix of dimensions 
%                         voxels x parcels with the same number of voxels 
%                         as S.D
%   S.mask              - mask for this parcellation
%   S.orthogonalisation - Orthogonalisation protocol to apply:
%                         ['none','symmetric','closest','householder']
%   S.method            - method for reconstructing parcel time course:
%                         ['PCA','mean','peakVoxel','spatialBasis']
% OUTPUTS:
%   D                   - Newly created SPM12 MEEG object
%
% Adam Baker 2015

D = spm_eeg_load(S.D);

S.orthogonalisation = ft_getopt(S,'orthogonalisation','none');
S.method            = ft_getopt(S,'method','PCA');


% Parcellation - pass in as a P*V matrix or 1*V vector
if ~isempty(strfind(S.parcellation,'.nii'))
    parcellation = readnii(S.parcellation,S.mask);
elseif ~isempty(strfind(S.parcellation,'.mat'))
    parcellation = load(S.parcellation);
    if length(fieldnames(parcellation)) == 1
        parcellation = parcellation.(char(fieldnames(parcellation)));
    else
        error('.mat file should contain only one variable')
    end
end


nodedata = get_node_tcs(D, parcellation, S.method);
nodedata = remove_source_leakage(nodedata, S.orthogonalisation);

good_samples = ~all(badsamples(D,':',':',':'));
data = zeros(size(nodedata,1),length(good_samples));
data(:,good_samples) = nodedata;
data = reshape(data,[size(data,1),D.nsamples,D.ntrials]);

clear voxeldata_concat nodedata_concat;

% Save data to new MEEG object
outfile = fullfile(D.path,D.fname);
Dnode = clone(montage(D,'switch',0),outfile,[size(data,1),D.nsamples,D.ntrials]);
Dnode = chantype(Dnode,1:Dnode.nchannels,'VE');
Dnode(:,:,:) = data;
D = Dnode; % For output
D.save;
end