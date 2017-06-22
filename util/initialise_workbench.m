function workbench_initialize()
	% Add workbench to system path
	
	s = osl_conf.read();

	if ~isfield(s,'WORKBENCH')
		s.WORKBENCH = '';
		osl_conf.write(s);
	end

	if isempty(s.WORKBENCH) && exist(s.WORKBENCH) && isempty(strfind(getenv('PATH'),s.WORKBENCH))
		if ismac
			workbench_path = fullfile(s.WORKBENCH,'bin_macosx64');
		else
			workbench_path = fullfile(s.WORKBENCH,'bin_linux64');
		end

		setenv('PATH',sprintf('%s%s%s',workbench_path,pathsep,getenv('PATH')));
		if system('hash wb_command')
			error('Workbench path is specified in osl.conf but it does not appear to be correct');
		end
	end




