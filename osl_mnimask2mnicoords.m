function [ mni_coords xform ] = osl_mnimask2mnicoords(mask_fname)
	% [ mni_coords xform ] = osl_mnimask2mnicoords(mask)
	%
	% converts an MNI standard brain mask, mask, into a list of mni_coords.
	% 
	% INPUTS
	% - mask_fname - Name of a nii file with volume data. Spatial resolution stored in header
	% - xstart - Optionally specify the MNI coordinates of first entry in each dimension
	%
	% Romesh Abeysuriya 2017
	% MW (pre-2014)

	% Extract the xform matrix from the nii header
	xform = runcmd('fslorient -getsform %s',mask_fname);
	xform = reshape(str2num(xform),4,4)';

	% Read in the mask
	mask = read_avw(mask_fname);

	% Convert nonzero entries in mask to array indices
	[x,y,z] = ind2sub(size(mask),find(mask));

	% Use xform to convert indices to the MNI coordinates
	mni_coords = [(x-1)*xform(1,1)+xform(1,4) (y-1)*xform(2,2)+xform(2,4) (z-1)*xform(3,3)+xform(3,4)];
