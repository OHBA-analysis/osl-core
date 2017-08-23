
% To reduce dependencies, and improve cross-platform compatibility, OSL attempts work within Matlab 
% as far as possible. To this end, reading and writing NIFTI files is done via a Matlab toolbox
% <LINK TO FILE EXCHANGE> 

% How to save
% How to load
% How to deal with xform


%% NIFTI FILES WITH UNUSUAL DATA TYPES
% The catch-22 is that loading with load_nii() will perform all rescalings AND
% also perform flipping. Subsequently saving with save_nii() could cause problems
% with orientation. Hence why is it safer to use read_avw
% [~,res,xform] = nii.load('original.nii.gz');
% vol = read_avw('original.nii.gz');
% nii.save(vol,res,xform,'new.nii.gz');