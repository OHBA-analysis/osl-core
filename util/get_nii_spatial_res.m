function [ mni_res ] = get_nii_spatial_res( mask_fname )
    % [ mni_res ] = get_nii_spatial_res( mask_fname )
    %
    % MWW 
    % RA - runcmd() will gracefully display a 'command not found' error if FSLDIR is not 
    % correctly set, and this is also checked by osl2_startup.m, so no need to check here
    % as well
    
    assert(exist(mask_fname,'file')~=0,'File %s not found',mask_fname);
    nii = load_untouch_nii(mask_fname);
    mni_res = nii.hdr.dime.pixdim(2:4);





