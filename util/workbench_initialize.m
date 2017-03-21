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
	else % Read conf file
		f = fopen(conf_file,'r');
		fgetl(f); % Skip header
		workbenchdir = fgetl(f); workbenchdir = regexp(workbenchdir,'.*= *','split'); workbenchdir = workbenchdir{2};
		fclose(f);
	end

	if ~exist(workbenchdir)
		fprintf(2,sprintf('Workbench is not present. Edit the contents of %s to specify the location of Workbench\n',conf_file))
	end

	workbench_path = fullfile(workbenchdir,'bin_macosx64');

	setenv( 'PATH', strjoin([ {workbench_path}, strsplit(getenv('PATH'),pathsep) ],pathsep) );

function [workbenchdir] = guess_workbenchdir()
	workbenchdir = '/Applications/workbench/';


