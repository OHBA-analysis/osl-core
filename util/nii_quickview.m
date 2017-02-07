function nii_quickview(mat,spat_res,rs)
% View nii file without creating (permanently) a new file
% AB 2012

OSLDIR = getenv('OSLDIR');

fname_tmp = tempname;

nii_quicksave(mat,fname_tmp,spat_res);

if exist('rs','var')
  stdbrain = [OSLDIR '/std_masks/MNI152_T1_' num2str(rs) 'mm_brain.nii.gz'];
  osl_resample_nii(fname_tmp, fname_tmp, rs, 'trilinear',stdbrain);
end
  
fslview(fname_tmp);

pause(5)
dos(['rm ' fname_tmp '.nii.gz']);

end