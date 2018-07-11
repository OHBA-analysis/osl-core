function osl_startup( osl_root )
    % Initialize OSL
    % osl_root is the folder containing 'osl-core'
    
    osl_core = fileparts(mfilename('fullpath')); % folder where this script is
    if nargin < 1 || isempty(osl_root) 
        osl_root = fileparts(osl_core);
    else
        osl_core = fullfile(osl_root,'osl-core');
    end
    
    assert( isdir(osl_root), 'Specified OSL directory does not exist: %s', osl_root );
    assert( isdir(osl_core), 'Could not find OSL core directory: %s', osl_core );

    % Check that OSL hasn't already started
    if ~isempty(getenv('OSLDIR'))
        warning('Found OSLDIR environment variable; shutting down before starting up again...');
        osl_shutdown();
    end
    setenv('OSLDIR',osl_root);
    fprintf(1,'[OSL] Starting up from folder: %s\n',getenv('OSLDIR'));

    % Check for manually specified OSLCONF, otherwise set default
    setenv('OSLCONF', find_oslconf() );
    fprintf(1,'[OSL] Using configuration file: %s\n',getenv('OSLCONF'));

    % Save current path
    % JH: use separate file for path backup
    % JH: DO NOT move this above
    oldpaths = backup_path();

    % does no path-changing if running in deployed mode (gw '13).
    if ~isdeployed 

        % Check and remove toolboxes that are supplied internally as part of OSL
        % TODO - add ROInets etc. to this list
        checklist={'fieldtrip', 'spm', 'osl', 'mne', 'netlab', 'fsl', 'fmt'};
        restoredefaultpath;

        % If anything goes wrong, osl_startup.m should still be left on the path
        addpath(osl_core);
        addpath(fullfile(osl_core,'util'));
        addpath(fullfile(osl_core,'util','jh'));

        % add old paths back, as long as they don't conflict with internal dependencies
        mpath = matlabroot;
        for j = 1:length(oldpaths)
            % skip matlab paths
            if ~isempty(strfind(oldpaths{j},mpath))
                continue; 
            end
            % if none of the words in the checklist is found in the current token, add old path back
            if all(cellfun( @(x) isempty(strfind(oldpaths{j},x)), checklist ))
                addpath(oldpaths{j});
            else
                if ~strfind(oldpaths{j},'osl') % Don't warn about OSL
                    fprintf(2,'Found and removed conflicting toolbox: %s\n',oldpaths{j});
                end
            end
        end

        % Check/add FSL binaries to the underlying system path, and Matlab functions to Matlab
        initialise_fsl();
        
        % Add Workbench
        initialise_workbench();
       
    end

    initialise_spm();

    % Add OHBA shared libraries
    extpath = fullfile(osl_root,'ohba-external');
    if ~isdir(extpath)
        fprintf(2,'Could not find "%s"\n',extpath);
        error('ohba-external is missing. Clone https://github.com/OHBA-analysis/ohba-external into the same directory as osl-core.');
    end

    addpath(extpath);
    ohba_external_startup();

    addpath(fullfile(osl_root,'GLEAN'));
    addpath(genpath_exclude(fullfile(osl_root,'HMM-MAR'),{'.git','.svn'}));
    addpath(fullfile(osl_root,'MEG-ROI-nets'));

    % Ensure osl-core directories get priority in path by adding it last
    addpath(genpath_exclude(osl_core,{'.git','.svn','spm-changes'}));
    addpath(osl_root);

end

function pathstr = genpath_exclude(pathstr,excludes)
    % Take in list of strings to exclude from path

    if ischar(excludes)
        excludes = {excludes};
    end

    paths = genpath(pathstr);
    paths = strsplit(paths,pathsep);
    retain = true(size(paths));

    for j = 1:length(excludes)
        retain = retain & cellfun( @(x) isempty(strfind(excludes{j},x)), paths );
    end

    paths = paths(retain);
    pathstr = strjoin(paths,pathsep);

end

function p = backup_path()

    p = path(); 

    fname = fullfile( getenv('OSLDIR'), 'tmp_path-backup' );
    setenv('OSL_PATH_BACKUP', fname);
    fh = fopen(fname,'w+');
    fwrite(fh,p);
    fclose(fh);

    p = strsplit(p,pathsep);

end

function f = find_oslconf()

    p = fullfile( getenv('OSLDIR'), 'osl.conf' );
    c = fullfile( getenv('OSLDIR'), 'osl-core', 'osl.conf' );

    % manually set
    if ~isempty(getenv('OSLCONF'))
        f = getenv('OSLCONF');
        return;
    end

    % in the osl-core directory
    if exist(c,'file') == 2
        f = c; 
        return;
    end 

    % otherwise, default to osl/osl.conf
    f = p;

end
