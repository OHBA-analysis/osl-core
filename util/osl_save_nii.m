function save_nii(vol,res,xform,fname)
	% Save a nii file together with a given xform matrix
	%
	% INPUTS
	% vol - volume matrix to save to nii file
	% fname - File name of nii file to save
	% xform - 4x4 matrix
	%
	% Romesh Abeysuriya 2017
	
	save_avw(vol,fname,'d',[res 1]);

	runcmd(['fslorient -setsformcode 0 ' fname])
	runcmd(['fslorient -setqformcode 2 ' fname])
	runcmd(['fslorient -setqform ' num2str(reshape(xform',1,16)) ' ' fname])
	runcmd(['fslorient -setsform ' num2str(reshape(xform',1,16)) ' ' fname])