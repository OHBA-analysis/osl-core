function save_nii(vol,res,xform,fname)
	% Save a nii file together with a given xform matrix
	%
	% INPUTS
	% vol - volume matrix to save to nii file
	% fname - File name of nii file to save
	% xform - 4x4 matrix
	%
	% Romesh Abeysuriya 2017
	
    if length(res)==1
        save_avw(vol,fname,'d',[res res res 1]);
    else
        save_avw(vol,fname,'d',[res 1]); %MWW
    end
    
	runcmd(['fslorient -setsformcode 0 ' fname])
	runcmd(['fslorient -setqformcode 2 ' fname])
	runcmd(['fslorient -setqform ' num2str(reshape(xform',1,16)) ' ' fname])
	runcmd(['fslorient -setsform ' num2str(reshape(xform',1,16)) ' ' fname])