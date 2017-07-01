function osl_save_nii(vol,res,xform,fname)
	% Save a nii file together with a given xform matrix
	%
	% INPUTS
	% vol - volume matrix to save to nii file
	% res - Spatial resolution e.g. '2' if same in all dimensions, or '[2 2 2]' to specify independently
	% fname - File name of nii file to save
	% xform - 4x4 matrix. 
	%
	% Could use osl_load_nii to get res and xform from a standard mask
	%
	% Romesh Abeysuriya 2017
	
    if length(res)==1
        save_avw(vol,fname,'d',[res res res 1]);
    else
        save_avw(vol,fname,'d',[res 1]); %MWW
    end

    runcmd(['fslorient -setqform ' num2str(reshape(xform',1,16)) ' ' fname])
    runcmd(['fslorient -setsform ' num2str(reshape(xform',1,16)) ' ' fname])
	runcmd(['fslorient -setsformcode 0 ' fname])
	runcmd(['fslorient -setqformcode 2 ' fname])
