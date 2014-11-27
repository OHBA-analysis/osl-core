function fname_out = nii_quicksave(mat,fname,spat_res,resamp,interp)
% nii_quicksave.m
%fname_out = nii_quicksave(mat,fname,spat_res,resamp)

global OSLDIR;
stdbrain=read_avw([OSLDIR '/std_masks/MNI152_T1_' num2str(spat_res) 'mm_brain.nii.gz']);
save_avw(matrix2vols(mat,stdbrain),fname,'f',[spat_res spat_res spat_res size(mat,2)]);

if ~exist('interp','var')
  interp = 'trilinear';
end

fname = strrep(fname,'.gz','');
fname = strrep(fname,'.nii','');


if exist('resamp','var')
  fname_rs=[fname '_ds' num2str(resamp) 'mm'];
  stdbrain = [OSLDIR '/std_masks/MNI152_T1_' num2str(resamp) 'mm_brain.nii.gz'];
  osl_resample_nii(fname, fname_rs, resamp,interp,stdbrain);
  dos(['rm ' fname '.nii.gz']);
  dos(['mv ' fname_rs '.nii.gz ' fname '.nii.gz']);
end

  fname_out=[fname '.nii.gz'];


end
