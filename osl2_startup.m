function osl2_startup( osldir )
    % Initialize OSL
    % osldir is the folder containing 'osl2'
    
    if nargin < 1 || isempty(osldir) 
        f = fileparts(mfilename('fullpath'));
        osldir = fileparts(f);
    end
    
    setenv('OSLDIR',osldir)

    % does no path-changing if running in deployed mode (gw '13).
    if ~isdeployed 

        % Check and remove toolboxes that are supplied internally as part of OSL
        checklist={'fieldtrip', 'spm', 'osl', 'mne', 'netlab', 'fsl', 'fmt'};
        oldpaths = strsplit(path,':');
        restoredefaultpath;

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
        
        % Check fsl has been set up
        if isempty(getenv('FSLDIR'))
            % Try and detect default location
            if exist('/usr/local/fsl') 
                fsldir = '/usr/local/fsl/';
                setenv('FSLDIR',fsldir)
                setenv('FSLOUTPUTTYPE','NIFTI_GZ')

                % Add it to the shell path from within MATLAB 
                curpath = getenv('PATH');
                setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),curpath));
            else
                error('FSL is not installed properly. Perhaps check that the $FSLDIR/bin directory is in your PATH before starting Matlab. See the Prerequisites section at https://sites.google.com/site/ohbaosl/download');
            end
        end

        % Try a dummy call to an fsl tool to make sure FSL is properly installed
        [status,res]=dos(['fslval']);
        if status~=1
            error('FSL is not installed properly. Perhaps check that the $FSLDIR/bin directory is in your PATH before starting Matlab. See the Prerequisites section at https://sites.google.com/site/ohbaosl/download');
        end;
       
        addpath(sprintf('%s/etc/matlab',getenv('FSLDIR')));
       
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Copy changes to SPM code from osl
    filelist={};targetdir={};

    filelist{end+1}='osl2/spm-beamforming-toolbox-osl-addons/bf_output_montage_osl.m';
    targetdir{end+1}='spm12/toolbox/spm-beamforming-toolbox';

    filelist{end+1}='osl2/spm-beamforming-toolbox-osl-addons/bf_write_spmeeg_osl.m';
    targetdir{end+1}='spm12/toolbox/spm-beamforming-toolbox';

    filelist{end+1}='osl2/spm-beamforming-toolbox-osl-addons/bf_inverse_mne_adaptive.m';
    targetdir{end+1}='spm12/toolbox/spm-beamforming-toolbox';

    filelist{end+1}='osl2/spm-changes/private/ft_read_event_4osl.m';
    targetdir{end+1}='spm12/external/fieldtrip/fileio';

    filelist{end+1}='osl2/spm-changes/private/read_trigger_4osl.m';
    targetdir{end+1}='spm12/external/fieldtrip/fileio/private';

    filelist{end+1}='osl2/spm-changes/private/badsamples.m';
    targetdir{end+1}='spm12/@meeg';

    filelist{end+1} = 'osl2/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/forward/private/';

    filelist{end+1} = 'osl2/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/plotting/private/';

    filelist{end+1} = 'osl2/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/fileio/private/';

    filelist{end+1} = 'osl2/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/utilities/private/';

    filelist{end+1} = 'osl2/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/private/';

    filelist{end+1} = 'osl2/spm-changes/private/ft_headmodel_localspheres.m';
    targetdir{end+1}='spm12/external/fieldtrip/forward/';

    filelist{end+1}='osl2/spm-changes/private/path.m';
    targetdir{end+1}='spm12/@meeg';

    filelist{end+1}='osl2/spm-changes/private/spm_eeg_montage.m';
    targetdir{end+1}='spm12';

    filelist{end+1} ='osl2/spm-changes/private/spm_eeg_inv_mesh_ui.m';
    targetdir{end+1}='spm12';

    filelist{end+1} ='osl2/spm-changes/private/subsref1.m';
    targetdir{end+1}='spm12/@meeg';

    filelist{end+1} ='osl2/spm-changes/private/ft_getopt.c';
    targetdir{end+1}='spm12/external/fieldtrip/src/';

    for k=1:length(filelist),
        copyfile( fullfile(osldir,filelist{k}), fullfile(osldir,targetdir{k}), 'f' );
    end

    % Add and initialize SPM
    addpath(fullfile(osldir,'spm12'))
    spm_get_defaults('cmdline',true);
    spm('Defaults','EEG')
    addpath(fullfile(osldir, 'spm12','external','fieldtrip','src'));

    % Add OHBA shared libraries
    if ~exist(fullfile(osldir,'osl-external'))
        fprintf(2,'Could not find ''%s''\n',fullfile(osldir,'osl-external'));
        error('osl-external is missing. Clone https://github.com/OHBA-analysis/osl-external into the same directory as osl2');
    end

    addpath(fullfile(osldir,'osl-external'));
    osl_external_startup

    % Ensure osl2 dir gets priority in path by adding it last
    addpath(genpath(fullfile(osldir,'osl2')))
    addpath(osldir)

    rmpath(fullfile(osldir,'osl2','spm-changes')); % These are already copied into spm

    % Remote fieldtrip substitutes for Matlab toolboxes if the toolbox is installed
    if license('test', 'Statistics_Toolbox')
        rmpath(fullfile(osldir,'spm12/external/fieldtrip/external/stats'))
    end

    if license('test','Signal_Toolbox')
        rmpath(fullfile(osldir,'spm12/external/fieldtrip/external/signal'))
    end
