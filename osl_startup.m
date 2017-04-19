function osl_startup( osldir )
    % Initialize OSL
    % osldir is the folder containing 'osl-core'
    
    if nargin < 1 || isempty(osldir) 
        f = fileparts(mfilename('fullpath'));
        osldir = fileparts(f);
    end
    
    if ~exist(osldir,'dir')
        error(sprintf('Specified OSL directory does not exist: %s',osldir));
    end

    setenv('OSLDIR',osldir)

    % does no path-changing if running in deployed mode (gw '13).
    if ~isdeployed 

        % Check and remove toolboxes that are supplied internally as part of OSL
        % TODO - add ROInets etc. to this list
        checklist={'fieldtrip', 'spm', 'osl', 'mne', 'netlab', 'fsl', 'fmt'};
        oldpaths = regexp(path,':','split');
        restoredefaultpath;

        % If anything goes wrong, osl_startup.m should still be left on the path
        addpath(fullfile(osldir,'osl-core'))
        addpath(fullfile(osldir,'osl-core','util'))

        for j = 1:length(oldpaths)
            if strfind(oldpaths{j},matlabroot)
                continue
            else
                if ~any(cellfun(@(x) ~isempty(strfind(oldpaths{j},x)),checklist))
                    addpath(oldpaths{j});
                else
                    if ~strfind(oldpaths{j},'osl') % Don't warn about OSL
                        fprintf(2,'Found and removed conflicting toolbox: %s\n',oldpaths{j});
                    end
                end
            end
        end

        % Check/add FSL binaries to the underlying system path, and Matlab functions to Matlab
        fsl_initialise() 
        
        % Add Workbench
        workbench_initialize()
       
    end

    if ~exist(fullfile(osldir,'spm12'),'dir')
        error(sprintf('SPM12 was not found at the required location: %s',fullfile(osldir,'spm12')))
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Copy changes to SPM code from osl
    filelist={};targetdir={};

    filelist{end+1}='osl-core/spm-changes/bf_save.m';
    targetdir{end+1}='spm12/toolbox/spm-beamforming-toolbox';

    filelist{end+1}='osl-core/spm-changes/bf_output_montage_osl.m';
    targetdir{end+1}='spm12/toolbox/spm-beamforming-toolbox';

    filelist{end+1}='osl-core/spm-changes/bf_write_spmeeg_osl.m';
    targetdir{end+1}='spm12/toolbox/spm-beamforming-toolbox';

    filelist{end+1}='osl-core/spm-changes/bf_inverse_mne_adaptive.m';
    targetdir{end+1}='spm12/toolbox/spm-beamforming-toolbox';

    filelist{end+1}='osl-core/spm-changes/private/ft_read_event_4osl.m';
    targetdir{end+1}='spm12/external/fieldtrip/fileio';

    filelist{end+1}='osl-core/spm-changes/private/read_trigger_4osl.m';
    targetdir{end+1}='spm12/external/fieldtrip/fileio/private';

    filelist{end+1}='osl-core/spm-changes/private/badsamples.m';
    targetdir{end+1}='spm12/@meeg';

    filelist{end+1} = 'osl-core/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/forward/private/';

    filelist{end+1} = 'osl-core/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/plotting/private/';

    filelist{end+1} = 'osl-core/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/fileio/private/';

    filelist{end+1} = 'osl-core/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/utilities/private/';

    filelist{end+1} = 'osl-core/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/private/';

    filelist{end+1} = 'osl-core/spm-changes/private/ft_headmodel_localspheres.m';
    targetdir{end+1}='spm12/external/fieldtrip/forward/';

    filelist{end+1}='osl-core/spm-changes/private/path.m';
    targetdir{end+1}='spm12/@meeg';

    filelist{end+1}='osl-core/spm-changes/private/spm_eeg_montage.m';
    targetdir{end+1}='spm12';

    filelist{end+1} ='osl-core/spm-changes/private/spm_eeg_inv_mesh_ui.m';
    targetdir{end+1}='spm12';

    filelist{end+1} ='osl-core/spm-changes/private/subsref1.m';
    targetdir{end+1}='spm12/@meeg';

    filelist{end+1} ='osl-core/spm-changes/private/ft_getopt.c';
    targetdir{end+1}='spm12/external/fieldtrip/src/';

    for k=1:length(filelist),
        copyfile( fullfile(osldir,filelist{k}), fullfile(osldir,targetdir{k}), 'f' );
    end

    % Add and initialize SPM
    addpath(fullfile(osldir,'spm12'))
    spm_get_defaults('cmdline',true);
    spm('Defaults','EEG')
    addpath(fullfile(spm('dir'),'matlabbatch')); % Required for cfg_getfile

    addpath(fullfile(osldir, 'spm12','external','fieldtrip','src'));
    addpath(fullfile(osldir, 'spm12','toolbox','spm-beamforming-toolbox'));

    % Add OHBA shared libraries
    if ~exist(fullfile(osldir,'ohba-external'))
        fprintf(2,'Could not find ''%s''\n',fullfile(osldir,'ohba-external'));
        error('ohba-external is missing. Clone https://github.com/OHBA-analysis/ohba-external into the same directory as osl-core');
    end

    addpath(fullfile(osldir,'ohba-external'));
    ohba_external_startup

    addpath(fullfile(osldir,'GLEAN'));
    addpath(genpath_exclude(fullfile(osldir,'HMM-MAR'),{'.git','.svn'}));
    addpath(fullfile(osldir,'MEG-ROI-nets'));

    % Ensure osl-core directories gets priority in path by adding it last
    addpath(genpath_exclude(fullfile(osldir,'osl-core'),{'.git','.svn','std_masks','docs'}))
    addpath(osldir)

    rmpath(fullfile(osldir,'osl-core','spm-changes')); % These are already copied into spm

    % Remote fieldtrip substitutes for Matlab toolboxes if the toolbox is installed
    if license('test', 'Statistics_Toolbox')
        rmpath(fullfile(osldir,'spm12/external/fieldtrip/external/stats'))
    end

    if license('test','Signal_Toolbox')
        rmpath(fullfile(osldir,'spm12/external/fieldtrip/external/signal'))
    end
end

function pathstr = genpath_exclude(pathstr,excludes)
    % Take in list of strings to exclude from path

    if ischar(excludes)
        excludes = {excludes};
    end

    paths = genpath(pathstr);
    paths = regexp(paths,':','split');

    retain = ones(size(paths));

    for j = 1:length(excludes)
        retain = retain & cellfun(@(x) isempty(regexp(x,excludes{j})),paths);
    end

    paths = paths(retain);
    pathstr = sprintf('%s:',paths{:});
    pathstr = pathstr(1:end-1); % Remove trailing delimiter
end


