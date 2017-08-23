function initialise_spm

	s = osl_conf.read();

	if ~isfield(s,'SPMDIR')
		s.SPMDIR = '';
		osl_conf.write(s);
	end

	if isempty(s.SPMDIR)
		SPMDIR = fullfile(osldir,'spm12');
    else
        SPMDIR = s.SPMDIR;
	end

	if ~exist(SPMDIR,'dir') || ~exist(fullfile(SPMDIR,'spm.m'),'file')
	    error(sprintf('SPM12 was not found at the required location: %s',fullfile(SPMDIR)))
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Copy changes to SPM code from osl
	filelist={};targetdir={};

	% Insert the beamforming toolbox
	rmdir(fullfile(SPMDIR,'toolbox','spm-beamforming-toolbox'),'s');
	copyfile(fullfile(osldir,'osl-core','spm-changes','spm-beamforming-toolbox'),fullfile(SPMDIR,'toolbox'));

	filelist{end+1}='osl-core/spm-changes/private/ft_read_event_4osl.m';
	targetdir{end+1}='external/fieldtrip/fileio';

	filelist{end+1}='osl-core/spm-changes/private/read_trigger_4osl.m';
	targetdir{end+1}='external/fieldtrip/fileio/private';

	% filelist{end+1}='osl-core/spm-changes/private/badsamples.m';
	% targetdir{end+1}='@meeg';

	filelist{end+1} = 'osl-core/spm-changes/private/undobalancing.m';
	targetdir{end+1}='external/fieldtrip/forward/private/';

	filelist{end+1} = 'osl-core/spm-changes/private/undobalancing.m';
	targetdir{end+1}='external/fieldtrip/plotting/private/';

	filelist{end+1} = 'osl-core/spm-changes/private/undobalancing.m';
	targetdir{end+1}='external/fieldtrip/fileio/private/';

	filelist{end+1} = 'osl-core/spm-changes/private/undobalancing.m';
	targetdir{end+1}='external/fieldtrip/utilities/private/';

	filelist{end+1} = 'osl-core/spm-changes/private/undobalancing.m';
	targetdir{end+1}='external/fieldtrip/private/';

	filelist{end+1} = 'osl-core/spm-changes/private/ft_headmodel_localspheres.m';
	targetdir{end+1}='external/fieldtrip/forward/';

	filelist{end+1}='osl-core/spm-changes/private/path.m';
	targetdir{end+1}='@meeg';

	filelist{end+1}='osl-core/spm-changes/private/spm_eeg_montage.m';
	targetdir{end+1}='';

	filelist{end+1} ='osl-core/spm-changes/private/spm_eeg_inv_mesh_ui.m';
	targetdir{end+1}='';

	filelist{end+1} ='osl-core/spm-changes/private/subsref1.m';
	targetdir{end+1}='@meeg';

	filelist{end+1} ='osl-core/spm-changes/private/ft_getopt.c';
	targetdir{end+1}='external/fieldtrip/src/';

	filelist{end+1} ='osl-core/spm-changes/private/ft_select_range.m';
	targetdir{end+1}='external/fieldtrip/plotting/';

	filelist{end+1} ='osl-core/spm-changes/private/topoplot_common.m';
	targetdir{end+1}='external/fieldtrip/private/';

	filelist{end+1} ='osl-core/spm-changes/private/ft_select_range.m';
	targetdir{end+1}='external/fieldtrip/plotting/';

	filelist{end+1} ='osl-core/spm-changes/private/ft_singleplotER.m';
	targetdir{end+1}='external/fieldtrip/';

	filelist{end+1} ='osl-core/spm-changes/private/ft_singleplotTFR.m';
	targetdir{end+1}='external/fieldtrip/';

	filelist{end+1} ='osl-core/spm-changes/private/functionSignatures.json';
	targetdir{end+1}='';

	for k=1:length(filelist),
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
	if license('test', 'Statistics_Toolbox')
	    rmpath(fullfile(osldir,'spm12/external/fieldtrip/external/stats'))
	end

	if license('test','Signal_Toolbox')
	    rmpath(fullfile(osldir,'spm12/external/fieldtrip/external/signal'))
	end