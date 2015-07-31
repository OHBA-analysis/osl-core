function D = glean_parcellation(S)
% Computes node time series from a source space MEEG object using a
% parcellation
%
% [D,assignments] = glean_parcellation(S)
%
% INPUTS:
%
% S.D                   - source space MEEG object
% S.parcellation        - .nii or .nii.gz file at the same gridstep as S.D
%                         OR a matrix of dimensions voxels x parcels
%                         with the same number of voxels as S.D
% S.mask                - mask for this parcellation
% S.orthogonalisation   - Orthogonalisation protocol to apply:
%                           ['none','symmetric','closest','householder']
% S.method              - method for reconstructing parcel time course:
%                           ['PCA','mean','peakVoxel','spatialBasis']


D = spm_eeg_load(S.D);

D = montage(D,'switch',2);

S.prefix            = 'p';
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


nodedata = ROInets.get_corrected_node_tcs(D, parcellation, S.orthogonalisation, S.method);


good_samples = ~all(badsamples(D,':',':',':'));



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