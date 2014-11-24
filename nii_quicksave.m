%% nii_quicksave.m
%
% converts from a matrix to a volume and saves a nii.gz file.
%
% fname_out = nii_quicksave(mat,fname,spat_res,resamp,interptype)
% fname is output file name (inc. path).
% spat_res in the resolution of the volume in mm.
% mat is the 2D matrix input.
%
% HL 2012

function fname_out = nii_quicksave(mat,fname,spat_res,resamp,interptype)

global OSLDIR;
stdbrain=read_avw([OSLDIR '/std_masks/MNI152_T1_' num2str(spat_res) 'mm_brain.nii.gz']);
save_avw(matrix2vols(mat,stdbrain),fname,'f',[spat_res spat_res spat_res size(mat,2)]);
if nargin <5
   interptype = 'trilinear';
end
if nargin >=4
    fname_rs=[fname '_ds' num2str(resamp) 'mm'];
    fname_out = osl_resample_nii(fname,fname_rs,resamp,interptype);
else
    fname_out=fname;
end
end
