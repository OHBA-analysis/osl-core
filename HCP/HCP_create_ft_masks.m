function HCP_create_ft_masks(res)
    % Iterate over 
    % Interpolate an MNI standard brain onto an HCP

    res = [4 5 6 8 10];

    for j = 1:length(res)
        input_sourcemodel = fullfile(osldir,'std_masks','fieldtrip_masks',sprintf('standard_sourcemodel3d%dmm.mat',res(j)));
        output_mask = fullfile(osldir,'std_masks',sprintf('ft2_%dmm_brain_mask.nii.gz',res(j)));
        output_brain = fullfile(osldir,'std_masks',sprintf('ft2_%dmm_brain.nii.gz',res(j)));
        std_brain =  fullfile(osldir,'std_masks',sprintf('MNI152_T1_%dmm_brain.nii.gz',res(j)));
        
        d = load(input_sourcemodel);
        HCP_sourcemodel_to_nii(d.sourcemodel,output_mask);
        nii.resample(std_brain,output_brain,output_mask,'interptype','nearest','enforce_mask',true,'force_positive',true);
        
        [vol,niires,xform] = nii.load(output_brain);
        mask = nii.load(output_mask);
        nii.save(vol + mask,niires,xform,output_brain); % Ensure all nonzero vals in mask have nonzero brain image values
    end

