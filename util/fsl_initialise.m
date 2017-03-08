function fsl_initialise(conf_file)
	% Adds FSL directories to system path, sets FSL environment variables
	% and adds Matlab utilities to the path
	%
	% Reads FSL directories from a configuration file to provide support for multiple
	% versions
	% 
	% Inputs:
	% 	conf_file - Configuration file, defaults to OSLDIR/fsl_location.txt
	%
	% If config file doesn't exist, one will be written with the guessed FSL location
	% 
	% If this file doesn't exist, try and guess t
	% Set up FSL environment variables reading in from a file
	%
	% If FSLDIR is already set and config file is empty, test FSL and return
	% This means nothing will change for users who have already set up the 
	% environment externally
	if nargin == 0 && ~isempty(getenv('FSLDIR'))
		addpath(fullfile(getenv('FSLDIR'),'etc','matlab')); % Make sure FSL matlab functions are on the path

		try
			test_fsl()
		catch ME
			error('fsl_initialize:error',sprintf('FSLDIR was already set, but the call to FSL failed. Try running ''setenv(''FSLDIR'')'' and running ''osl_startup'' again\nMessage: %s',ME.message))
		end
		return
	end

	% Set default conf file location
	if nargin < 1 || isempty(conf_file) 
		conf_file = fullfile(getenv('OSLDIR'),'fsl_location.txt');
	end

	if ~exist(conf_file) % If conf file doesn't exist, write default locations AND use them immediately if they exist
		[fsldir,bindir,libdir] = guess_fsldir();
		f = fopen(conf_file,'w');
		fprintf(f,'# Specify location of FSL here. If environment variable FSLDIR is set prior to starting Matlab, this file is ignored\n');
		fprintf(f,'FSLDIR=%s\n',fsldir);
		fprintf(f,'BIN=%s\n',bindir);
		fprintf(f,'LIB=%s\n',libdir);
		fclose(f)
	else % Read conf file
		f = fopen(conf_file,'r');
		fgetl(f); % Skip header
		fsldir = fgetl(f); fsldir = regexp(fsldir,'.*= *','split'); fsldir = fsldir{2};
		bindir = fgetl(f); bindir = regexp(bindir,'.*= *','split'); bindir = bindir{2};
		libdir = fgetl(f); libdir = regexp(libdir,'.*= *','split'); libdir = libdir{2};
		fclose(f);
	end

	if ~exist(fsldir)
		fprintf(2,'FSLDIR=%s\nFSLLIB=%s\n',fsldir,libdir)
		error(sprintf('FSL directory does not exist. Check FSL is installed, then edit the contents of %s to specify the location of FSL',conf_file))
	end

	setenv('FSLDIR',fsldir)
	setenv('FSLOUTPUTTYPE','NIFTI_GZ')
	setenv('PATH',sprintf('%s:%s',bindir,getenv('PATH')));
	setenv('LD_LIBRARY_PATH',sprintf('%s:%s',libdir,getenv('LD_LIBRARY_PATH')));
	addpath(fullfile(fsldir,'etc','matlab'));

	test_fsl()

function test_fsl()
	% Does a test call to FSL and will throw an error if it doesn't work
	assert(~isempty(getenv('FSLDIR')),'FSLDIR was not set')
	assert(exist(getenv('FSLDIR'))~=0,'FSLDIR does not exist')
	assert(~isempty(which('read_avw')),'FSL matlab utilities not found')

	[status,res] = system('fslval');
	if status~=1
	    error('FSL is not installed properly. Perhaps check that the $FSLDIR/bin directory is in your PATH before starting Matlab. See the Prerequisites section at https://sites.google.com/site/ohbaosl/download');
	end

function [fsldir,bindir,libdir] = guess_fsldir()
	% If the user hasn't specified a conf file, try and guess where FSL is

	fsldir = '/usr/local/fsl';
	bindir = '/usr/local/fsl/bin';
	libdir = '/usr/local/fsl/lib';

	% Debian FSL package
	if ~exist(fsldir) &&  exist('/usr/share/fsl/5.0')
		fsldir = '/usr/share/fsl/5.0';
		bindir = '/usr/share/fsl/5.0/bin';
		libdir = '/usr/lib/fsl/5.0';
	end


