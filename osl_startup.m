function osl_startup( osl_root )
    % Initialize OSL
    % osl_root is the folder containing 'osl-core'
    
    if nargin < 1 || isempty(osl_root) 
        f = fileparts(mfilename('fullpath'));
        osl_root = fileparts(f);
    end
    
    if ~exist(osl_root,'dir')
        error(sprintf('Specified OSL directory does not exist: %s',osl_root));
    end

    % Back up original path
    path_backup = path;

    setenv('OSLDIR',osl_root)

    % does no path-changing if running in deployed mode (gw '13).
    if ~isdeployed 

        % Check and remove toolboxes that are supplied internally as part of OSL
        % TODO - add ROInets etc. to this list
        checklist={'fieldtrip', 'spm', 'osl', 'mne', 'netlab', 'fsl', 'fmt'};
        oldpaths = strsplit(path,pathsep);
        restoredefaultpath;

        % If anything goes wrong, osl_startup.m should still be left on the path
        addpath(fullfile(osl_root,'osl-core'))
        addpath(fullfile(osl_root,'osl-core','util'))
        addpath(fullfile(osl_root,'osl-core','util','jh'))

        for j = 1:length(oldpaths)
            if strfind(oldpaths{j},matlabroot)
                continue
            else
                % if none of the words in the checklist is found in the current token, add it back
                if all(cellfun( @(x) isempty(strfind(oldpaths{j},x)), checklist ))
                    addpath(oldpaths{j});
                else
                    if ~strfind(oldpaths{j},'osl') % Don't warn about OSL
                        fprintf(2,'Found and removed conflicting toolbox: %s\n',oldpaths{j});
                    end
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

    % Ensure osl-core directories gets priority in path by adding it last
    addpath(genpath_exclude(fullfile(osl_root,'osl-core'),{'.git','.svn','spm-changes'}))
    addpath(osl_root)

    % Save backed up path
    s = osl_conf.read();
    s.PATH_BACKUP = path_backup;
    osl_conf.write(s);

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


