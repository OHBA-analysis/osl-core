function D = osl_update_inv_dir(D,newpath)
	% Update folder containing files for inverse model
	%
	% >> D.inv{1}.mesh.tess_ctx
	% ans =
	%     '/home/disk3/ajquinn/Projects/drugface/structurals/M10/structcortex_8196.surf.gii'
	% D2 = osl_update_inv_dir(D,'testdir')
	% >> D2.inv{1}.mesh.tess_ctx
	% ans =
	%     'testdir/structcortex_8196.surf.gii'
		
	for j = 1:length(D.inv)
		field = fields(D.inv{j}.mesh);
		for k = 1:length(field)
			if isstr(D.inv{j}.mesh.(field{k}))
				[~,fname,ext] = fileparts(D.inv{j}.mesh.(field{k}));
				D.inv{j}.mesh.(field{k}) = fullfile(newpath,[fname ext]);
			end
		end
	end

