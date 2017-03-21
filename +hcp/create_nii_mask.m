function create_nii_mask(sourcemodel,mask_out)
	% Creates an nii mask of an HCP sourcemodel at the same resolution.
	%
	% HCP_create_nii_mask(sourcemodel,mask_out)
	%
	% sourcemodel   - HCP sourcemodel struct
	% mask_out      - filename of a nifti file to save to
	%
	% Romesh Abeysuriya 2017
	% Adam Baker 2015

	% Compute spatial resolution
	res = 10*[diff(sourcemodel.xgrid(1:2)) diff(sourcemodel.ygrid(1:2)) diff(sourcemodel.zgrid(1:2))]
	offset = 10*min(sourcemodel.pos(sourcemodel.inside,:));

	% Convert position to grid resolution offset, and then shift to integer values
	idx = bsxfun(@rdivide,10*sourcemodel.pos(sourcemodel.inside,:),res);
	idx = round(bsxfun(@minus,idx+1,min(idx)));

	% Index into a 3D mask
	vol = zeros(sourcemodel.dim);
	for j = 1:size(idx,1)
		vol(idx(j,1),idx(j,2),idx(j,3)) = 1;
	end  
	
	% Save the nii file
	save_avw(vol,mask_out,'i',[res 1]);

	% Add in xform data - also enables recovery of coordinates with osl_mnimask2mnicoords
	xform = [res(1)  0       0       offset(1),...
	         0      res(2)   0       offset(2),...
	         0      0       res(3)   offset(3),...
	         0      0       0       1];

	runcmd(['fslorient -setsformcode 0 ' mask_out])
	runcmd(['fslorient -setqformcode 2 ' mask_out])
	runcmd(['fslorient -setsform ' num2str(xform) ' ' mask_out])
	runcmd(['fslorient -setqform ' num2str(xform) ' ' mask_out])

