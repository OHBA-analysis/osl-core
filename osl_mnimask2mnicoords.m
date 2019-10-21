function [ mni_coords, xform ] = osl_mnimask2mnicoords(fname)
    % [ mni_coords xform ] = osl_mnimask2mnicoords(mask)
    %
    % converts an MNI standard brain mask into a list of mni_coords.
    % 
    % INPUTS
    % - fname - Name of a nii file with volume data. Spatial resolution stored in header
    %
    % Romesh Abeysuriya 2017
    % MW (pre-2014)
    
    [vol,~,xform] = nii.load(fname); 
    [x,y,z] = ind2sub(size(vol),find(vol)); % Convert nonzero entries in mask to array indices
    mni_coords = xform*[(x-1)';(y-1)';(z-1)';ones(size(x'))]; % Use xform to convert indices to the MNI coordinates
    mni_coords = mni_coords(1:3,:)';
