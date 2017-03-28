function HCP_resample_nii(input_fname,output_fname)
% Resamples nifti file to 2mm resolution on the MNI space grid. Requires
% that the original nifti was saved with proper orientation (e.g. using
% HCP_nii_quicksave with an appropriate HCP template nifti).
%
% HCP_resample_nii(input_fname,output_fname)
%
% input_fname   - input nii file at HCP gridstep
% output_fname  - output nii file in 2mm MNI space
%
% Adam Baker 2015


std_brain = [getenv('FSLDIR'),'/data/standard/MNI152_T1_2mm_brain.nii'];
std_mask  = [getenv('FSLDIR'),'/data/standard/MNI152_T1_2mm_brain_mask.nii'];

% Transformation from HCP to MNI space:
T = eye(4);
T(1,4) = 14;
T(2,4) = 14;
trans = [tempname '.mat'];
dlmwrite(strrep(trans,'.mat','.txt'),T,' ');
runcmd(['mv ' strrep(trans,'.mat','.txt') ' ' trans]);

% Dilate to avoid edge effects:
[~,kernelsize] = runcmd(['fslval ' input_fname ' pixdim1']);
kernelsize = str2double(kernelsize);
runcmd(['fslmaths ' input_fname ' -kernel gauss ' kernelsize ' -dilM ' output_fname]);

% Flirt to MNI space:
runcmd(['flirt -in ' input_fname ' -applyxfm -init ' trans ' -out ' output_fname ' -paddingsize 0.0 -interp ' 'trilinear' ' -ref ' std_brain]);
runcmd(['rm ' trans]);

% Mask to remove edges:
runcmd(['fslmaths ' output_fname ' -mas ' std_mask ' ' output_fname]);

end