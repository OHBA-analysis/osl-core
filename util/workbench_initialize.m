function workbench_initialize(conf_file)
	% Add workbench to system path
	%
	% Minimal testing for now - much more fragile than fsl_initialize

	% Set default conf file location
	if nargin < 1 || isempty(conf_file) 
		conf_file = fullfile(getenv('OSLDIR'),'workbench_location.txt');
	end

	if ~exist(conf_file) % If conf file doesn't exist, write default locations AND use them immediately if they exist
		[workbenchdir,caretdir] = guess_workbenchdir();
		f = fopen(conf_file,'w');
		fprintf(f,'# Specify location of WORKBENCH and CARET here\n');
		fprintf(f,'WORKBENCH=%s\n',workbenchdir);
		fprintf(f,'CARET=%s\n',caretdir);
		fclose(f);
	else % Read conf file
		f = fopen(conf_file,'r');
		fgetl(f); % Skip header
		workbenchdir = fgetl(f); workbenchdir = regexp(workbenchdir,'.*= *','split'); workbenchdir = workbenchdir{2};
		caretdir = fgetl(f); caretdir = regexp(caretdir,'.*= *','split'); caretdir = caretdir{2};
		fclose(f);
	end

	if ~exist(workbenchdir)
		fprintf(2,sprintf('Workbench is not present. Edit the contents of %s to specify the location of Workbench\n',conf_file))
	end

	if ~exist(caretdir)
		fprintf(2,sprintf('Caret is not present. Edit the contents of %s to specify the location of Caret\n',conf_file))
	end

	workbench_path = fullfile(workbenchdir,'bin_macosx64');
	caret_path = fullfile(caretdir,'bin_macosx64');

	setenv( 'PATH', strjoin([ {workbench_path, caret_path}, strsplit(getenv('PATH'),pathsep) ],pathsep) );

function [workbenchdir,caretdir] = guess_workbenchdir()
	workbenchdir = '/Applications/workbench/';
	caretdir = '/Applications/caret/';



