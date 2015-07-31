function HCP_nii_quicksave(map,fname,mask_fname,resample)
% Write data to nifti file, based on an input mask
%
% HCP_nii_quicksave(map,fname,mask_fname,resample)
%
% map           - voxels x tpts matrix
% fname         - filename of nifti file to save to
% mask_fname    - an appropriate HCP mask at the same grid spacing as map
% resample      - resample to 2mm MNI grid [0/1] (default 0)
%
% Adam Baker 2015

if nargin < 4
    resample = 0;
end

[~,~,ext] = fileparts(fname);
if ~strcmp(ext,'.gz')
    fname = [fname '.gz'];
end

[mask,~,scales] = read_avw(mask_fname);
save_avw(matrix2vols(map,mask),fname,'f',scales);

system(['fslcpgeom ' mask_fname ' ' fname ' -d']);

if resample == 1
    HCP_resample_nii(fname,fname);
end

end