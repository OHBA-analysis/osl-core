function fsl_initialise()
	% Adds FSL directories to system path, sets FSL environment variables
	% and adds Matlab utilities to the path
	%
	% Reads FSL directories from the OSL configuration file to provide support for multiple
	% versions
	% 
	% Specifically, it uses the following variables in osl.conf:
	%
	% FSLDIR
	% FSLBIN
	% FSLLIB
	
	s = osl_conf.read();

	% Are the directories set in the conf file? If not, write defaults
	[default_FSLDIR,default_FSLBIN,default_FSLLIB] = guess_fsldir();
	
	if ~isfield(s,'FSLDIR') || isempty(s.FSLDIR)
		s.FSLDIR = default_FSLDIR;
	end

	if ~isfield(s,'FSLBIN') || isempty(s.FSLBIN)
		s.FSLBIN = default_FSLBIN;
	end
	
	if ~isfield(s,'FSLLIB') || isempty(s.FSLLIB)
		s.FSLLIB = default_FSLLIB;
	end

	osl_conf.write(s);

	% Initialize FSL
	setenv('FSLDIR',s.FSLDIR)
	setenv('FSLOUTPUTTYPE','NIFTI_GZ')
	if isempty(strfind(getenv('PATH'),s.FSLBIN))
		setenv('PATH',sprintf('%s%s%s',s.FSLBIN,pathsep,getenv('PATH')));
	end
	if isempty(strfind(getenv('LD_LIBRARY_PATH'),s.FSLLIB))
		setenv('LD_LIBRARY_PATH',sprintf('%s%s%s',s.FSLLIB,pathsep,getenv('LD_LIBRARY_PATH')));
	end
	addpath(fullfile(s.FSLDIR,'etc','matlab'));

	% Test FSL
	if ~exist(getenv('FSLDIR'))
        fprintf(2,'FSLDIR does not exist. Check that it is set correctly in osl.conf\n');
    end
    
	if isempty(which('read_avw'))
        fprintf(2,'FSL matlab utilities not found\n')
    end
    
	[status,res] = system('fslval');
	if status~=1
	    fprintf(2,'FSL is not installed properly. Have you installed the correct version and set it in osl.conf?\n');
	end


function [fsldir,bindir,libdir] = guess_fsldir()
	% Guess where the FSL directories maybe

	fsldir = '/usr/local/fsl';
	bindir = '/usr/local/fsl/bin';
	libdir = '/usr/local/fsl/lib';

	% Debian FSL package
	if ~exist(fsldir) &&  exist('/usr/share/fsl/5.0')
		fsldir = '/usr/share/fsl/5.0';
		bindir = '/usr/share/fsl/5.0/bin';
		libdir = '/usr/lib/fsl/5.0';
	end


