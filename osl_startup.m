function osl_startup( osl_root )
%
% Initialize OSL
% Optionally accepts a single argument osl_root, which is the folder CONTAINING 'osl-core'.
%
    
    osl_core = fileparts(mfilename('fullpath')); % folder where this script is
    if nargin < 1 || isempty(osl_root) 
        osl_root = fileparts(osl_core);
    end
    
    assert( osl_util.isdir(osl_root), 'Specified OSL directory does not exist: %s', osl_root );
    assert( osl_util.isdir(osl_core), 'Could not find OSL core directory: %s', osl_core );

    % Check that OSL hasn't already started
    if osl_isactive()
        warning('Found OSLDIR environment variable; shutting down before starting up again...');
        osl_shutdown();
    end
    setenv('OSLDIR',osl_root);
    fprintf(1,'[OSL] Starting up from folder: %s\n',getenv('OSLDIR'));

    % Check for manually specified OSLCONF, otherwise set default
    setenv('OSLCONF', osl_conf.find() );
    fprintf(1,'[OSL] Using configuration file: %s\n',getenv('OSLCONF'));

    % Save current path
    % JH: use separate file for path backup
    % JH: DO NOT move this above
    oldpaths = osl_conf.path_backup();

    % does no path-changing if running in deployed mode (gw '13).
    if ~isdeployed 

        % Check and remove toolboxes that are supplied internally as part of OSL
        % TODO - add ROInets etc. to this list
        checklist={'fieldtrip', 'spm', 'osl', 'mne', 'netlab', 'fsl', 'fmt'};
        restoredefaultpath;

        % If anything goes wrong, osl_startup.m should still be left on the path
        addpath(osl_core);
        addpath(fullfile(osl_core,'util'));

        % add old paths back, as long as they don't conflict with internal dependencies
        mpath = matlabroot;
        for j = 1:length(oldpaths)
            % skip matlab paths
            if osl_util.contains(oldpaths{j},mpath)
                continue; 
            end
            % if none of the words in the checklist is found in the current token, add old path back
            if all(cellfun( @(x) ~osl_util.contains(oldpaths{j},x), checklist ))
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
    if ~osl_util.isdir(extpath)
        error([ ...
            'Folder ohba-external is missing (not found: %s).' newline ...
            'Clone https://github.com/OHBA-analysis/ohba-external into %s.'], extpath, osl_root );
    end

    addpath(extpath);
    ohba_external_startup();

    addpath(fullfile(osl_root,'GLEAN'));
    addpath(osl_util.genpath_exclude(fullfile(osl_root,'HMM-MAR'),{'.git','.svn'}));
    addpath(fullfile(osl_root,'MEG-ROI-nets'));

    % Ensure osl-core directories get priority in path by adding it last
    local = {'','africa','examples','HCP','oat','oil','opt','osl_hmm_toolbox','oslview','rhino','util'};
    local = cellfun( @(x) fullfile(osl_core,x), local, 'UniformOutput', false );
    addpath(local{:});

end
