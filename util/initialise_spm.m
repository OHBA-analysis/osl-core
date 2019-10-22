function initialise_spm

    disable_undobalancing = true; % Default is true, to disable undobalancing and preserve third order gradients
    
    s = osl_conf.read();

    if isempty(s.SPMDIR)
        SPMDIR = osldir('spm12');
    else
        SPMDIR = s.SPMDIR;
    end

    if ~osl_util.isfile(fullfile(SPMDIR,'spm.m'))
        error('SPM12 was not found at the required location: %s',fullfile(SPMDIR))
    end

    if disable_undobalancing
        undobalancing_fname = 'undobalancing_giles.m';
    else
        undobalancing_fname = 'undobalancing_original.m';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Copy changes to SPM code from osl
    filelist={};targetdir={};

    % Insert the beamforming toolbox - only do this once
    % The alternative would be to overwrite this directory every time
    % The current code makes it possible that someone has an old copy of this toolbox
    % However, if the folder is deleted every time, this makes it very unreliable to use
    % OSL on a cluster with a shared OSL folder
    bf_path = fullfile(SPMDIR,'toolbox','spm-beamforming-toolbox');
    if isdir(bf_path)
        rmdir(bf_path,'s');
    end
    copyfile(fullfile(osldir,'osl-core','spm-changes','spm-beamforming-toolbox'),bf_path);

    filelist{end+1}  = fullfile('osl-core','spm-changes','ft_read_event_4osl.m');
    targetdir{end+1} = fullfile('external','fieldtrip','fileio');

    filelist{end+1}  = fullfile('osl-core','spm-changes','read_trigger_4osl.m');
    targetdir{end+1} = fullfile('external','fieldtrip','fileio','private');

    filelist{end+1}  = fullfile('osl-core','spm-changes',undobalancing_fname);
    targetdir{end+1} = fullfile('external','fieldtrip','forward','private','undobalancing.m');

    filelist{end+1}  = fullfile('osl-core','spm-changes',undobalancing_fname);
    targetdir{end+1} = fullfile('external','fieldtrip','plotting','private','undobalancing.m');

    filelist{end+1}  = fullfile('osl-core','spm-changes',undobalancing_fname);
    targetdir{end+1} = fullfile('external','fieldtrip','fileio','private','undobalancing.m');

    filelist{end+1}  = fullfile('osl-core','spm-changes',undobalancing_fname);
    targetdir{end+1} = fullfile('external','fieldtrip','utilities','private','undobalancing.m');

    filelist{end+1}  = fullfile('osl-core','spm-changes',undobalancing_fname);
    targetdir{end+1} = fullfile('external','fieldtrip','private','undobalancing.m');

    filelist{end+1}  = fullfile('osl-core','spm-changes','ft_headmodel_localspheres.m');
    targetdir{end+1} = fullfile('external','fieldtrip','forward');

    % filelist{end+1}='osl-core/spm-changes/spm_eeg_montage.m';
    % targetdir{end+1}='';

    filelist{end+1}  = fullfile('osl-core','spm-changes','spm_eeg_inv_mesh_ui.m');
    targetdir{end+1} = '';

    % filelist{end+1} ='osl-core/spm-changes/ft_getopt.c';
    % targetdir{end+1}='external/fieldtrip/src/';

    % filelist{end+1} ='osl-core/spm-changes/ft_select_range.m';
    % targetdir{end+1}='external/fieldtrip/plotting/';

    % filelist{end+1} ='osl-core/spm-changes/topoplot_common.m';
    % targetdir{end+1}='external/fieldtrip/private/';

    % filelist{end+1} ='osl-core/spm-changes/ft_select_range.m';
    % targetdir{end+1}='external/fieldtrip/plotting/';

    % filelist{end+1} ='osl-core/spm-changes/ft_singleplotER.m';
    % targetdir{end+1}='external/fieldtrip/';

    % filelist{end+1} ='osl-core/spm-changes/ft_singleplotTFR.m';
    % targetdir{end+1}='external/fieldtrip/';

    filelist{end+1} = fullfile('osl-core','spm-changes','functionSignatures.json');
    targetdir{end+1}= '';

    for k=1:length(filelist)
        copyfile( fullfile(osldir,filelist{k}), fullfile(SPMDIR,targetdir{k}), 'f' );
    end

    % Add and initialize SPM
    addpath(SPMDIR)
    spm_get_defaults('cmdline',true);
    spm('Defaults','EEG')
    addpath(fullfile(spm('dir'),'matlabbatch')); % Required for cfg_getfile

    addpath(fullfile(SPMDIR,'external','fieldtrip','src'));
    addpath(fullfile(SPMDIR,'toolbox','spm-beamforming-toolbox'));

    % Remove fieldtrip substitutes for Matlab toolboxes if the toolbox is installed
    if license('test', 'Statistics_Toolbox') && osl_util.contains(path,fullfile(SPMDIR,'external','fieldtrip','external','stats'))
        rmpath(fullfile(SPMDIR,'external','fieldtrip','external','stats'))
    end

    if license('test','Signal_Toolbox') && osl_util.contains(path,fullfile(SPMDIR,'external','fieldtrip','external','signal'))
        rmpath(fullfile(SPMDIR,'external','fieldtrip','external','signal'))
    end
    
end