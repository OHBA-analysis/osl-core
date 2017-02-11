function ok = osl_check_installation(do_log)
	% Run diagnostic tests

	if nargin < 1 || isempty(do_log) 
		do_log = true;
	end
	
	if do_log
		f = fopen('osl_debug_log.txt','w');
		fprintf(1,'Writing diagnostic record to %s\n',fullfile(pwd,'osl_debug_log.txt'));
		log = @(x) fprintf(f,'%s\n',x); % Write a string to the log file
		section = @(x) log(sprintf('\n------------------ %s ------------------',upper(x)));
	else
		log = @(x) 0;
		section = @(x) 0;
	end


	ok = true;

	

	log('OSL DEBUG LOG');
	log(datetime('now'));

	% Check OS
	section('System Information');
	try
		import	java.lang.*;
		import	com.sun.security.auth.module.*;
		import	java.net.*;
		log(char(System.getProperty('os.arch')));
		log(char(System.getProperty('os.name')));
		log(char(System.getProperty('os.version')));
	catch ME
		log('Error retrieving system information');
		log(ME.identifier);
	end

	% Check Matlab version
	section('Matlab Version');
	log(sprintf('Matlab version: %s',version));
	arrayfun(@(x) log(sprintf('%s - %s',x.Name,x.Version)),ver);

	% Check that expected directory structure is present
	section('Directory structure');
    osldir = fileparts(fileparts(mfilename('fullpath')));
    dirs = {'osl2','spm12','GLEAN','HMM-MAR','layouts','MEG-ROI-nets','ohba-external'};
    for j = 1:length(dirs)
    	d = fullfile(osldir,dirs{j});
    	if ~exist(d)
    		log(sprintf('MISSING: %s',d));
    	else
    		log(sprintf('Located: %s',d));
    	end
    end

	% Check SPM version

	% Check that FSL is installed
	section('FSL');
	fsldir = getenv('FSLDIR');
	if isempty(fsldir)
		log('FSL environment variable not set')
		if exist('/usr/local/fsl/')
			log('Found /usr/local/fsl/, attemping automatic configuration');
			fsldir = '/usr/local/fsl/';
			setenv('FSLDIR',fsldir)
			setenv('FSLOUTPUTTYPE','NIFTI_GZ')
			curpath = getenv('PATH');
			setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),curpath));
		end
	else
		log(sprintf('FSLDIR: %s',fsldir));
	end

	[status,res] = system('fslval');

	log(sprintf('FSL return status = %d',status));
	log(sprintf('FSL result = %s',res));



	% Check if any mex things are required

	if do_log
		fclose(f);
	end



