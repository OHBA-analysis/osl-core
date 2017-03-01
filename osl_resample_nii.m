function output_fname = osl_resample_nii(input_fname, output_fname, out_gridstep, interp, output_mask_fname, dilate)

% output_fname = osl_resample_nii(input_fname, output_fname, out_gridstep, interp)
% output_fname = osl_resample_nii(input_fname, output_fname, out_gridstep, interp, output_mask_fname)
%
% Uses FSL's flirt to resample a niftii file input_fname to the pixdims
% [out_gridstep out_gridstep out_gridstep]
% interp is the interpolation used (e.g. 'trilinear' or
% 'nearestneighbour');
% the ouput can be optionally masked using niftii file output_mask_fname
%
% MWW 2012
OSLDIR = getenv('OSLDIR');

if(nargin<6)
    dilate=1;
end;

if(nargin<5)
    output_mask_fname='';
end;

if(nargin<4)
    interp='trilinear';
    output_mask_fname='';
end;

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

%runcmd(['fslcreatehd ' num2str(round([91 109 91]*2/out_gridstep)) ' 1 ' num2str(out_gridstep) ' ' num2str(out_gridstep) ' ' num2str(out_gridstep) ' 1 0 0 0 4 ' input_fname '_tmp.nii.gz']); 
%runcmd(['flirt -in ' input_fname ' -applyxfm -init ' getenv('FSLDIR') '/etc/flirtsch/ident.mat -out ' output_fname ' -paddingsize 0.0 -interp ' interp ' -ref ' input_fname '_tmp']);

if(dilate),
    if(~isempty(output_mask_fname)) % AB - dilate volume to avoid edge effects
      runcmd(['fslmaths ' input_fname ' -kernel gauss ' num2str(pixdim(1)) ' -dilM ' output_fname]);
      input_fname = output_fname;
    end
end;

runcmd(['flirt -in ' input_fname ' -applyxfm -init ' getenv('FSLDIR') '/etc/flirtsch/ident.mat -out ' output_fname ' -paddingsize 0.0 -interp ' interp ' -ref ' OSLDIR '/std_masks/MNI152_T1_' num2str(out_gridstep) 'mm_brain.nii.gz']);

%runcmd(['rm -f ' input_fname '_tmp.nii.gz']);

if(~isempty(output_mask_fname))
    runcmd(['fslmaths ' output_fname ' -mas ' output_mask_fname ' ' output_fname]);
end;


%runcmd(['fslhd ' output_fname '.nii.gz >& ' output_fname '_image_format.txt']);
%runcmd(['fslhd ' OSLDIR '/std_masks/MNI152_T1_' num2str(out_gridstep) 'mm_brain >& ' input_fname '_mni_format.txt']);
%runcmd(['diff ' input_fname '_mni_format.txt ' output_fname '_image_format.txt']);
