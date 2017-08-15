function [ mask_indices_in_lower_level, current_level_mni_coord, new_mask ] = setup_mask_indices( Sin )

% [ mask_indices_in_lower_level ] = setup_mask_indices( Sin )
%
%
% sets up current level mask and returns mask_indices_in_lower_level, which is the voxel index in the lower level mask 
% of the voxels in the current level mask 
%
% MWW 2011

OSLDIR = getenv('OSLDIR');

if(isfield(Sin,'lower_level_results')),
    Sin.lower_level_gridstep=Sin.lower_level_results.gridstep;    
    Sin.lower_level_mni_coord=Sin.lower_level_results.mni_coord;
end;

new_mask=1;

current_level=Sin.current_level;

if(isfield(current_level,'mask_fname')) % current level mask also provided

    gridstep=Sin.lower_level_gridstep;

    % setup lower level mask coords
    [ lower_level_mni_coord xform ] = osl_mnimask2mnicoords(Sin.lower_level_mask_fname);

    disp(['Basing mask on high res mask ' current_level.mask_fname]);

    % setup current level mask at grid resolution 
    % NOTE - update from osl_resample_nii to nii.resample is not tested!
    % If there are problems, try nii.resample_flirt which should be the same as the old osl_resample_nii
    
    % current_level_mask_fname_lowres = osl_resample_nii(current_level.mask_fname, [current_level.mask_fname '_' num2str(gridstep) 'mm.nii.gz'], gridstep, 'sinc', [osldir '/std_masks/MNI152_T1_' num2str(gridstep) 'mm_brain.nii.gz']);
    current_level_mask_fname_lowres = nii.resample(current_level.mask_fname, [current_level.mask_fname '_' num2str(gridstep) 'mm.nii.gz'], gridstep, 'interptype','linear','enforce_mask',true);

    % only use mask that intersects with the lower level mask
    runcmd(['fslmaths ' current_level_mask_fname_lowres ' -mas ' Sin.lower_level_mask_fname ' -thr 0.05 ' current_level_mask_fname_lowres]);
    
    % setup mni coords using std space mask
    [ current_level_mni_coord, xform ] = osl_mnimask2mnicoords(current_level_mask_fname_lowres);

     disp(['Using low res mask ' current_level_mask_fname_lowres]);
     
    [c, iA, mask_indices]=intersect(current_level_mni_coord, lower_level_mni_coord,'rows');

    % sort indices to be in the order of the mask mni_coords:
    [ff gg]=sort(iA);

    mask_indices_in_lower_level=mask_indices(gg);

    if(length(mask_indices_in_lower_level)==0)
      error('No dipoles fall in the provided ROI mask, using nearest to centre of mask');          
    end;

    str=['Using dipoles at MNI coordinates ' ];

    for vox=1:length(mask_indices_in_lower_level),
        str=[str ', [' num2str(lower_level_mni_coord(mask_indices_in_lower_level(vox),:)) ']'];
    end;

    disp(str);
    
    str=['Using dipoles at MNI coordinates ' ];

    for vox=1:size(current_level_mni_coord,1),
        str=[str ', [' num2str(current_level_mni_coord(vox,:)) ']'];
    end;

    disp(str);
    
    % save current level mask to disk at low res
    mask=nii.load(current_level_mask_fname_lowres); 
    
    % WAS: save_avw(mask,Sin.current_level.mask_fname, 'f',
    % [gridstep,gridstep,gridstep,1]);  % causes error in first level when
    % using an ROI mask (trys to overwrite 2mm
    % Sin.current_level.mask_fname)!!

    nii.save(mask,[gridstep,gridstep,gridstep,1],[],Sin.current_level_mask_fname); 

else % no current level mask provided, so use lower level mask (if provided)

    gridstep=Sin.lower_level_gridstep;

    disp(['Using same indices as lower_level']);
    mask_indices_in_lower_level=1:size(Sin.lower_level_mni_coord,1); % indices for the output of this level that correspond to the niftii mask file 

    current_level_mni_coord=Sin.lower_level_mni_coord;
    
    new_mask=0;
    
    if isfield(Sin,'lower_level_mask_fname') && exist(Sin.lower_level_mask_fname)
        % save current level mask as lower level mask to disk if it exists
        copyfile(Sin.lower_level_mask_fname,Sin.current_level_mask_fname);
    else
        disp('No current level mask, just working with MNI coords.');
    end
end

