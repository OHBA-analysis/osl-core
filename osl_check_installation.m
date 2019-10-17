function osl_check_installation(do_log,test_fslview)
    % Run diagnostic tests
    % Some tests copied from spm_check_installation (this file uses \b control characters which makes it hard to log to txt)
    % TODO - should this return an OK status?

    if nargin < 2 || isempty(test_fslview) 
        test_fslview = true;
    end
    
    if nargin < 1 || isempty(do_log) 
        do_log = true;
    end
    
    osldir = fileparts(fileparts(mfilename('fullpath')));

    % Set up printing functions
    log = @(x) fprintf(1,'%s\n',x); % Format message
    log_error = @(x,ME) fprintf(1,'%s\n%s\n%s\n',x,ME.identifier,ME.message);
    section = @(x) log(sprintf('\n------------------ %s ------------------',upper(x)));

    if do_log
        logfile = fullfile(osldir,'osl_debug_log.txt');
        if exist(logfile)
            delete(logfile);
        end
        diary(logfile)
        fprintf(1,'Writing diagnostic record to %s\n',logfile);
    end

    log('OSL DEBUG LOG');
    log(sprintf('Date: %04d-%02d-%02d %02d:%02d:%02d',fix(clock())))

    section('OSL version')
    log(sprintf('osldir = %s',osldir))
    version_file = fullfile(osldir,'version.txt');
    if osl_util.isfile(version_file)
        type(version_file)
    else
        log('Version file missing')
    end

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

    if verLessThan('matlab','8.4')
        log('FAIL - Matlab versions earlier than R2014b may have graphics errors and have limited support for OSL');
    else
        log('PASS - Matlab version is supported')
    end

    if verLessThan('matlab','7.14')
        log('FAIL - Matlab versions earlier than R2012a are not officially supported by OSL');
    end

    % Check that expected directory structure is present
    section('Directory structure');
    dirs = {'osl-core','spm12','GLEAN','HMM-MAR','layouts','MEG-ROI-nets','ohba-external','example_data','parcellations','std_masks'};
    for j = 1:length(dirs)
        d = fullfile(osldir,dirs{j});
        if ~exist(d)
            log(sprintf('MISSING: %s',d));
        else
            log(sprintf('Located: %s',d));
        end
    end

    input_mask = fullfile(osldir,'std_masks','MNI152_T1_8mm_brain.nii.gz');

    section('NIFTI tools');
    try
        [a,b,c] = nii.load(input_mask);
        log('PASS - nii.load')
    catch ME
        log(sprintf('FAIL - nii.load\n%s',ME.message));
    end

    try
        nii.save(a,b,c,'test.nii.gz');
        log('PASS - nii.save')
    catch ME
        log(sprintf('FAIL - nii.save\n%s',ME.message));
    end

    % Check that FSL is installed
    section('FSL');
    fsldir = getenv('FSLDIR');
    if isempty(fsldir)
        log('FSL environment variable not set')
        log('Attempting automatic configuration')
        try
            addpath(fullfile(osldir,'osl-core','util'));
            fsl_initialise;
        catch ME
            log_error('Error initialising FSL',ME)
        end
    else
        log(sprintf('FSLDIR: %s',fsldir));
    end

    try
        [status,res] = system('fslval');
        if status == 1
            log('PASS - fslval')
        else
            log(sprintf('FAIL - fslval. Return status = %d',status));
            log(sprintf('FAIL - fslval. Result = %s',res));
        end
    catch ME
        log(sprintf('FAIL - fslval\n%s',ME.message));
    end

    try
        runcmd('fslmaths %s -thr 100 osl_fslmaths_test.nii.gz',input_mask);
        log('PASS - fslmaths')
        delete('osl_fslmaths_test.nii.gz')
    catch ME
        log(sprintf('FAIL - fslmaths\n%s',ME.message));
    end

    if test_fslview
        try
            fslview(input_mask)
            log('PASS - fslview. Check that fslview window has appeared');
        catch ME
            log(sprintf('FAIL - fslview\n%s',ME.message));
        end
    else
        log('Skipping fslview')
    end
        
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

    section('Fieldtrip')

    try
        a = ft_getopt(struct('a',1),'a');
        log('PASS - ft_getopt()');

    catch ME
        log_error('Error testing ft_getopt',ME)
    end

    section('Folder permissions')

    try
        [~,a]=fileattrib(spm('dir'));
        if a.UserWrite == 1
            log('PASS - SPMDIR is writable');
        else
            log(sprintf('FAIL - cannot write to SPMDIR %s',spm('dir')));
        end
    catch ME
        log_error('Error testing SPMDIR write attribute',ME)
    end

    try
        [~,a]=fileattrib(fullfile(osldir,'osl.conf'));
        if a.UserWrite == 1
            log('PASS - osl.conf is writable');
        else
            log('FAIL - cannot write to osl.conf');
        end
    catch ME
        log_error('Error testing osl.conf write attribute',ME)
    end

    section('MAC OS XQUARTZ');
    if ~ismac
        log('PASS - Not required on this system');
    else
        try
            if exist('/Applications/Utilities/XQuartz.app','dir')
                log('PASS - XQuartz present');
            else
                log('FAIL - Expected /Applications/Utilities/XQuartz.app to exist');
            end
        catch ME
            log_error('Error checking XQuartz',ME)
        end
    end

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

