function [ mni_res ] = get_nii_spatial_res( mask_fname )

% [ mni_res ] = get_nii_spatial_res( mask_fname )
%
% MWW 

res=[];

[status,res]=dos(['fslval ' mask_fname ' pixdim1']);
try
    mni_res(1)=str2double(res);
    assert(isfinite(mni_res(1)));
catch
    if(strcmp(getenv('FSLDIR'),'')),
        error('FSLDIR environmental variable is not set, please make sure that FSL is fully installed.');
    end;
    
    disp(res);
    
    error(['Invalid file: ' mask_fname]);
end;

[status,res]=dos(['fslval ' mask_fname ' pixdim2']);mni_res(2)=str2num(res);
[status,res]=dos(['fslval ' mask_fname ' pixdim3']);mni_res(3)=str2num(res);

if(status>1)
    disp(res);
    error('FSL is not installed properly. Perhaps check that the $FSLDIR/bin directory is in your PATH before starting Matlab');
end;

if(mni_res(1)~=mni_res(2) | mni_res(2)~=mni_res(3))
    error('all pixel dimensions should be the same for the input nii');
end;  

end

