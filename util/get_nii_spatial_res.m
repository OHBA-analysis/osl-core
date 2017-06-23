function [ mni_res ] = get_nii_spatial_res( mask_fname )
    % [ mni_res ] = get_nii_spatial_res( mask_fname )
    %
    % MWW 
    % RA - runcmd() will gracefully display a 'command not found' error if FSLDIR is not 
    % correctly set, and this is also checked by osl2_startup.m, so no need to check here
    % as well

    res=[];

    assert(exist(mask_fname,'file')~=0,'File %s not found',mask_fname);

    res = runcmd(['fslval ' mask_fname ' pixdim1']);
    mni_res(1)=str2double(res);
    assert(isfinite(mni_res(1)));

    res=runcmd(['fslval ' mask_fname ' pixdim2']);
    mni_res(2)=str2num(res);

    res=runcmd(['fslval ' mask_fname ' pixdim3']);
    mni_res(3)=str2num(res);

    if(mni_res(1)~=mni_res(2) | mni_res(2)~=mni_res(3))
        error('all pixel dimensions should be the same for the input nii');
    end 



