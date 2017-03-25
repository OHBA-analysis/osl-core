function [vol,res,xform] = load_nii(fname)
	% Load a nii volume together with the spatial resolution and xform matrix
	%
	% Romesh Abeysuriya 2017
	[vol,~,scales] = read_avw(fname);
	res = scales(1:3).';
	xform = reshape(str2num(runcmd('fslorient -getsform %s',fname)),4,4)';
