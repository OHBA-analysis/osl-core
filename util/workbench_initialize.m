function workbench_initialize(conf_file)
	% Add workbench to system path
	%
	% Minimal testing for now - much more fragile than fsl_initialize

	% Set default conf file location
	if nargin < 1 || isempty(conf_file) 
		conf_file = fullfile(getenv('OSLDIR'),'workbench_location.txt');
	end

	if ~exist(conf_file) % If conf file doesn't exist, write default locations AND use them immediately if they exist
		[workbenchdir] = guess_workbenchdir();
		f = fopen(conf_file,'w');
		fprintf(f,'# Specify location of WORKBENCH here\n');
		fprintf(f,'WORKBENCH=%s\n',workbenchdir);
		fclose(f);

		if ~exist(workbenchdir)
			fprintf(2,sprintf('Workbench is not present. Edit the contents of %s to specify the location of Workbench\n',conf_file))
		end
	else % Read conf file
		f = fopen(conf_file,'r');
		fgetl(f); % Skip header
		workbenchdir = fgetl(f); workbenchdir = regexp(workbenchdir,'.*= *','split'); workbenchdir = strtrim(workbenchdir{2});
		fclose(f);
	end

	if exist(workbenchdir) && isempty(strfind(getenv('PATH'),workbenchdir))
		if ismac
			workbench_path = fullfile(workbenchdir,'bin_macosx64');
		else
			workbench_path = fullfile(workbenchdir,'bin_linux64');
		end
		setenv( 'PATH', strjoin([ {workbench_path}, strsplit(getenv('PATH'),pathsep) ],pathsep) );
	end

function [workbenchdir] = guess_workbenchdir()
	workbenchdir = '/Applications/workbench/';



