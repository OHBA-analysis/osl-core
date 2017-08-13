function output_fname = osl_resample_nii(input_fname, output_fname, output_ref, interp, output_mask_fname, dilate)
% Resample a nii file using FSL's flirt
%
% output_fname = osl_resample_nii(input_fname, output_fname, output_ref, interp)
% output_fname = osl_resample_nii(input_fname, output_fname, output_ref, interp, output_mask_fname)
%
% Uses FSL's flirt to resample a nifti file onto a reference grid
% 
% INPUTS
% - input_fname : The original image file. All pixdims in this file must be the same
% - output_fname : File that will contain the resampled image
% - output_ref : Specifies the output resolution. If a number (e.g. 8) then this will specify one 
%                of the OSL standard masks. Alternatively, you can provide a .nii file that will
%                be passed as the reference image to flirt
% - interp : Specify interpolation method (e.g. 'trilinear' or 'nearestneighbour')
% - output_mask_fname : Optionally mask the image using fslmaths using this .nii file
% - dilate : Avoid edge artifacts by applying a dilation filter prior to resampling
%
% RA 2017
% MWW 2012

if(nargin<6) || isempty(dilate)
    dilate=1;
end;

if(nargin<5) || isempty(output_mask_fname)
    output_mask_fname='';
end;

if(nargin<4) || isempty(interp)
    interp='trilinear';
end;

if isnumeric(output_ref)
    output_ref = [osldir '/std_masks/MNI152_T1_' num2str(output_ref) 'mm_brain.nii.gz'];
else
    output_ref = output_ref;
end

res = runcmd(['fslval ' input_fname ' pixdim1']);
try
    pixdim(1)=str2num(res);
catch 
    error(['Invalid file: ' input_fname]);
end;

res = runcmd(['fslval ' input_fname ' pixdim2']);pixdim(2)=str2num(res);
res = runcmd(['fslval ' input_fname ' pixdim3']);pixdim(3)=str2num(res);

if(pixdim(1)~=pixdim(2) | pixdim(2)~=pixdim(3))
    error('all pixel dimensions should be the same for the input nii');
end;

%runcmd(['fslcreatehd ' num2str(round([91 109 91]*2/output_ref)) ' 1 ' num2str(output_ref) ' ' num2str(output_ref) ' ' num2str(output_ref) ' 1 0 0 0 4 ' input_fname '_tmp.nii.gz']); 
%runcmd(['flirt -in ' input_fname ' -applyxfm -init ' getenv('FSLDIR') '/etc/flirtsch/ident.mat -out ' output_fname ' -paddingsize 0.0 -interp ' interp ' -ref ' input_fname '_tmp']);

if(dilate),
    if(~isempty(output_mask_fname)) % AB - dilate volume to avoid edge effects
      runcmd(['fslmaths ' input_fname ' -kernel gauss ' num2str(pixdim(1)) ' -dilM ' output_fname]);
      input_fname = output_fname;
    end
end;

runcmd(['flirt -in ' input_fname ' -applyxfm -init ' getenv('FSLDIR') '/etc/flirtsch/ident.mat -out ' output_fname ' -paddingsize 0.0 -interp ' interp ' -ref ' output_ref]);

%runcmd(['rm -f ' input_fname '_tmp.nii.gz']);

if(~isempty(output_mask_fname))
    runcmd(['fslmaths ' output_fname ' -mas ' output_mask_fname ' ' output_fname]);
end;


%runcmd(['fslhd ' output_fname '.nii.gz >& ' output_fname '_image_format.txt']);
%runcmd(['fslhd ' output_ref ' >& ' input_fname '_mni_format.txt']);
%runcmd(['diff ' input_fname '_mni_format.txt ' output_fname '_image_format.txt']);
