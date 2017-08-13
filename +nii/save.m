function save(vol,res,xform,fname)
	% Save a nii file together with a given xform matrix
	%
	% INPUTS
	% vol - volume matrix to save to nii file
	% res - Spatial resolution
	% fname - File name of nii file to save
	% xform - 4x4 matrix. 
	%
	% Could use osl_load_nii to get res and xform from a standard mask
	%
	% Romesh Abeysuriya 2017
	
	% Resolution can be specified in 3 ways
	% - Single number = same in all dimensions, with time resolution of 1
	% - 3 numbers, different in all dimensions, time resolution of 1
	% - 4 numbers, complete resolution specification
	%
	% To check that the use of load_untouch_nii and save_nii is consistent:
	%
	%	x = read_avw(fullfile(osldir,'std_masks','MNI152_T1_8mm_brain.nii.gz'));
	%	[y,res,xform] = nii.load(fullfile(osldir,'std_masks','MNI152_T1_8mm_brain.nii.gz'));
	%	nii.save(y,res,xform,fullfile(osldir,'std_masks','test.nii.gz'))
	%	z = read_avw(fullfile(osldir,'std_masks','test.nii.gz'));
	%	all(x(:)==z(:))
	
    assert(all(size(xform)==[4 4]),'xform must be a 4x4 matrix')
    
	switch length(res)
		case 1
			r = [res res res];
		case 3
			r = [res];
		case 4
			r = res(1:3);
		otherwise
			error('Unknown resolution - should be 1, 3, or 4 elements long');
    end

    
    nii = make_nii(vol,res);
    nii.hdr.hist.qform_code = 0;
    nii.hdr.hist.sform_code = 4;
    nii.hdr.hist.srow_x = xform(1,:);
    nii.hdr.hist.srow_y = xform(2,:);
    nii.hdr.hist.srow_z = xform(3,:);
    save_nii(nii,fname);
