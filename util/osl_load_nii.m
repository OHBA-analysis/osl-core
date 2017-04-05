function [vol,res,xform] = load_nii(fname)
	% Load a nii volume together with the spatial resolution and xform matrix
	%
	% Romesh Abeysuriya 2017
	assert(ischar(fname),'Input file name must be a string')
	assert(logical(exist(fname,'file')),sprintf('Requested file "%s" could not be found',fname))
	[vol,~,scales] = read_avw(fname);
	res = scales(1:3).';
	xform = reshape(str2num(runcmd('fslorient -getsform %s',fname)),4,4)';
