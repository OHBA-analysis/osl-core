function osl_check_installation(do_log)
	% Run diagnostic tests
	% Some tests copied from spm_check_installation (this file uses \b control characters which makes it hard to log to txt)
	% TODO - should this return an OK status?

	if nargin < 1 || isempty(do_log) 
		do_log = true;
	end
	
	% Set up printing functions
	log = @(x) fprintf(1,'%s\n',x); % Format message
	log_error = @(x,ME) fprintf(1,'%s\n%s\n%s',x,ME.identified,ME.message);
	section = @(x) log(sprintf('\n------------------ %s ------------------',upper(x)));

	if do_log
		try
			delete('osl_debug_log.txt');
		end
		diary('osl_debug_log.txt')
		fprintf(1,'Writing diagnostic record to %s\n',fullfile(pwd,'osl_debug_log.txt'));
	end

	log('OSL DEBUG LOG');
	log(sprintf('Date: %04d-%02d-%02d %02d:%02d:%02d',fix(clock())))

	% Check OS
	section('System Information');
	try
		spm_system_tests()
	catch ME
		log_error('Error getting system information',ME)
	end

	% Check Matlab version
	section('Matlab Version');
	log(sprintf('Matlab version: %s',version));
	log(sprintf('Matlab location: %s',matlabroot));

	arrayfun(@(x) log(sprintf('%s - %s',x.Name,x.Version)),ver);

	% Check that expected directory structure is present
	section('Directory structure');
    osldir = fileparts(fileparts(mfilename('fullpath')));
    dirs = {'osl2','spm12','GLEAN','HMM-MAR','layouts','MEG-ROI-nets','ohba-external','example_data','parcellations'};
    for j = 1:length(dirs)
    	d = fullfile(osldir,dirs{j});
    	if ~exist(d)
    		log(sprintf('MISSING: %s',d));
    	else
    		log(sprintf('Located: %s',d));
    	end
    end

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

	input_mask = fullfile(osldir,'std_masks','MNI152_T1_8mm_brain.nii.gz');

	try
		runcmd('fslmaths %s -thr 100 osl_fslmaths_test.nii.gz',input_mask);
		log('PASS - fslmaths')
	catch ME
		log(sprintf('FAIL - fslmaths\n%s',ME.message));
	end

	try
		m = read_avw(input_mask);
	catch ME
		log('PASS - read_avw')
		log(sprintf('FAIL - read_avw\n%s',ME.message));
	end

	try
		fslview(input_mask)
		log('fslview successful - window should have appeared');
	catch ME
		log(sprintf('FAIL - fslview\n%s',ME.message));
	end
	
	% try
	% 	log('FSLERRORREPORT OUTPUT')
	% 	system('fslerrorreport');
	% catch ME
	% 	log_error('Could not run fslerrorreport()',ME)
	% end

	% Check SPM 
	section('SPM')

	if exist(fullfile(osldir,'spm12'))
		try
			addpath(fullfile(osldir,'spm12'))
			spm_get_defaults('cmdline',true);
			spm('Defaults','EEG');
			spm_check()
		catch ME
			log_error('Error testing SPM',ME)
		end
	else
		log('CANNOT TEST SPM12 BECAUSE DIRECTORY IS MISSING');
	end

	% Check if any mex things are required
	if do_log
		diary off
	end

function spm_system_tests()
	% from spm_check_installation.m

	%-Detect Platform and Operating System
	%--------------------------------------------------------------------------
	[C, maxsize] = computer;
	fprintf('Platform: %s (maxsize=%d)\n', C, maxsize);
	if ispc
	   platform = [system_dependent('getos'),' ',system_dependent('getwinsys')];
	elseif ismac
	    [fail, input] = unix('sw_vers');
	    if ~fail
	    platform = strrep(input, 'ProductName:', '');
	    platform = strrep(platform, sprintf('\t'), '');
	    platform = strrep(platform, sprintf('\n'), ' ');
	    platform = strrep(platform, 'ProductVersion:', ' Version: ');
	    platform = strrep(platform, 'BuildVersion:', 'Build: ');
	    else
	        platform = system_dependent('getos');
	    end
	else    
	   platform = system_dependent('getos');
	end
	fprintf('OS: %s\n', platform);

	%-Detect Java
	%--------------------------------------------------------------------------
	fprintf('%s\n', version('-java'));
	fprintf('Java support: ');
	level = {'jvm', 'awt', 'swing', 'desktop'};
	for i=1:numel(level)
	    if isempty(javachk(level{i})), fprintf('%s ',level{i}); end
	end
	fprintf('\n');

	%-Detect Monitor(s)
	%--------------------------------------------------------------------------
	M = get(0,'MonitorPositions');
	fprintf('Monitor(s):');
	for i=1:size(M,1)
	    fprintf(' [%d %d %d %d]',M(i,:));
	end
	fprintf(' (%dbit)\n', get(0,'ScreenDepth'));

	%-Detect OpenGL rendering
	%--------------------------------------------------------------------------
	S =  opengl('data');
	fprintf('OpenGL version: %s',S.Version);
	if S.Software, fprintf('(Software)\n'); else fprintf('(Hardware)\n'); end
	fprintf('OpenGL renderer: %s (%s)\n',S.Vendor,S.Renderer);

	%-Detect MEX setup
	%--------------------------------------------------------------------------
	fprintf('MEX extension: %s\n',mexext);
	try
	    cc = mex.getCompilerConfigurations('C','Selected');
	    if ~isempty(cc)
	        cc = cc(1); % can be C or C++
	        fprintf('C Compiler: %s (%s).\n', cc.Name, cc.Version);
	        fprintf('C Compiler settings: %s (''%s'')\n', ...
	            cc.Details.CompilerExecutable, cc.Details.OptimizationFlags);
	    else
	        fprintf('No C compiler is selected (see mex -setup)\n');
	    end
	end
	try
	    [sts, m] = fileattrib(fullfile(SPMdir,'src'));
	    m = [m.UserRead m.UserWrite m.UserExecute ...
	         m.GroupRead m.GroupWrite m.GroupExecute ...
	         m.OtherRead m.OtherWrite m.OtherExecute];
	    r = 'rwxrwxrwx'; r(~m) = '-';
	    fprintf('C Source code permissions: dir %s, ', r);
	    [sts, m] = fileattrib(fullfile(SPMdir,'src','spm_resels_vol.c'));
	    m = [m.UserRead m.UserWrite m.UserExecute ...
	         m.GroupRead m.GroupWrite m.GroupExecute ...
	         m.OtherRead m.OtherWrite m.OtherExecute];
	    r = 'rwxrwxrwx'; r(~m) = '-';
	    fprintf('file %s\n',r);
	end

function spm_check()
	% from spm_check_installation.m

	SPMdir = which('spm.m','-ALL');
	if isempty(SPMdir)
	    fprintf('SPM is not in your MATLAB path.\n');
	    return;
	elseif numel(SPMdir) > 1
	    fprintf('SPM seems to appear in several different folders:\n');
	    for i=1:numel(SPMdir)
	        fprintf('  * %s\n',SPMdir{i});
	    end
	    fprintf('Remove all but one with ''pathtool'' or ''spm_rmpath''.\n');
	    return;
	else
	    fprintf('SPM is installed in: %s\n',fileparts(SPMdir{1}));
	end
	SPMdir = fileparts(SPMdir{1});

	%-Detect SPM version and revision number
	%--------------------------------------------------------------------------
	v = struct('Name','','Version','','Release','','Date','');
	try
	    fid = fopen(fullfile(SPMdir,'Contents.m'),'rt');
	    if fid == -1
	        fprintf('Cannot open ''%s'' for reading.\n',fullfile(SPMdir,'Contents.m'));
	        return;
	    end
	    l1 = fgetl(fid); l2 = fgetl(fid);
	    fclose(fid);
	    l1 = strtrim(l1(2:end)); l2 = strtrim(l2(2:end));
	    t  = textscan(l2,'%s','delimiter',' '); t = t{1};
	    v.Name = l1; v.Date = t{4};
	    v.Version = t{2}; v.Release = t{3}(2:end-1);
	catch
	    fprintf('Cannot obtain SPM version & revision number.\n');
	    return;
	end
	fprintf('SPM version is %s (%s, %s)\n', ...
	    v.Release,v.Version,strrep(v.Date,'-',' '));

