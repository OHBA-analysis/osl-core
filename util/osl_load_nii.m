function [vol,res,xform] = osl_load_nii(fname)
	% Load a nii volume together with the spatial resolution and xform matrix
	%
	% Romesh Abeysuriya 2017
	assert(ischar(fname),'Input file name must be a string')
	assert(logical(exist(fname,'file')),sprintf('Requested file "%s" could not be found',fname))

	nii = load_untouch_nii(fname);
	vol = nii.img;
	res = nii.hdr.dime.pixdim(2:4);
	xform = eye(4);

	if nii.hdr.hist.sform_code == 0 &&  nii.hdr.hist.qform_code == 0
		fprintf(2,'Warning - *NO* qform/sform header code is provided! Your NIFTI file is *not* usable without this information!\n')
	elseif nii.hdr.hist.sform_code == 0
		fprintf(2,'Warning - osl_load_nii() currently only supports sform (not qform)!\n')
	end

	xform(1,:) = nii.hdr.hist.srow_x;
	xform(2,:) = nii.hdr.hist.srow_y;
	xform(3,:) = nii.hdr.hist.srow_z;




