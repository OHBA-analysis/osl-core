function HCP_create_ft_masks(res)
	% Iterate over 
	% Interpolate an MNI standard brain onto an HCP

	res = [4 5 6 8 10];

	for j = 1:length(res)
		input_sourcemodel = fullfile(osldir,'std_masks','fieldtrip_masks',sprintf('standard_sourcemodel3d%dmm.mat',res(j)));
		output_mask = fullfile(osldir,'std_masks',sprintf('ft_%dmm_brain_mask.nii.gz',res(j)));
		output_brain = fullfile(osldir,'std_masks',sprintf('ft_%dmm_brain.nii.gz',res(j)));
		std_brain =  fullfile(osldir,'std_masks',sprintf('MNI152_T1_%dmm_brain.nii.gz',res(j)));
		
		d = load(input_sourcemodel);
		HCP_sourcemodel_to_nii(d.sourcemodel,output_mask);
		osl_resample_nii_matlab(std_brain,output_mask,output_brain,'interptype','nearest','enforce_mask',true,'force_positive',true);
		
		[vol,niires,xform] = osl_load_nii(output_brain);
		mask = osl_load_nii(output_mask);
		osl_save_nii(vol + mask,niires,xform,output_brain); % Ensure all nonzero vals in mask have nonzero brain image values
	end

