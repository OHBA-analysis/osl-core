function fif_out = osl_maxfilter(fif_in,method,varargin)

	if nargin < 2 || isempty(method) 
		method = '';
	end
	
	arg = inputParser;
	arg.addParameter('maxfilter_path','/neuro/bin/util');% Where to find maxfilter
	arg.addParameter('movecomp',true); % Always movecomp by default
	arg.addParameter('autobad',false); % Never autobad by default
	arg.addParameter('st_win',4)
	arg.addParameter('st_corr',0.98)
	arg.addParameter('badchannels',{},@iscell); % Cell array of channel names e.g. D.chanlabels(D.badchannels)
	arg.addParameter('downsample_factor',4)
	arg.addParameter('cal_file',[]) % :  full path to a fine calibration file to use
	arg.addParameter('ctc_file',[]) % ctc_file (optional):  full path to a cross-talk matrix file to use
	arg.addParameter('trans_ref_file',[]) % transformation to specified reference
	arg.addParameter('extra',[]) % Just stick stuff at the end
	arg.parse(varargin{:});

	% Handle the input file
	[dirpath,fname,ext] = fileparts(fif_in);
	if isempty(ext)
		ext = '.fif';
	end

	maxfilter_call = fullfile(arg.Results.maxfilter_path,'maxfilter');

	% Nested function to accumulate maxfilter call string
	function add_to_call(str,varargin)
		str = sprintf(str,varargin{:}); % Handle substitutions
		maxfilter_call = sprintf('%s %s',maxfilter_call,str); % Add to accumulator
	end

	switch method
		case ''
			fif_suffix = '_nosss'
			add_to_call('-nosss');
			if arg.Results.downsample_factor == 1
				fprintf(2,'Not running SSS but also not downsampling - HPI signal may still be present in output\n');
			end
			fprintf('Calling maxfilter WITHOUT SSS\n');
		case 'sss'
			fif_suffix = '_sss'
			fprintf('Calling maxfilter with SSS\n');
		case 'tsss'
			fif_suffix = '_tsss'
			add_to_call('-st %d -corr %g',arg.Results.st_win,arg.Results.st_corr);
			fprintf('Calling maxfilter with tSSS\n');
		otherwise
			error('Unrecognized method')
	end

	base_path = fullfile(dirpath,[fname fif_suffix]); % Output, with no extension
	fif_out = [base_path '.fif'];
	log_out = [base_path '_log.txt'];
	head_out = [base_path '_headpos.txt'];

	add_to_call('-f %s',fif_in); % Specify input
	add_to_call('-o %s',fif_out); % Specify output
	add_to_call('-hp %s',head_out); % Always output HPI 


	if ~isempty(arg.Results.badchannels)
		% Remove MEG from channel names if present
		badchans = cellfun(@(x) strrep(x,'MEG',''),arg.Results.badchannels,'UniformOutput',false);
		badchans = sprintf('%s ',badchans{:}); % Join them
		add_to_call('-bad %s',badchans);
	end

	if arg.Results.movecomp
		add_to_call('-movecomp');
	end

	if ~arg.Results.autobad
		add_to_call('-autobad off');
	end

	if arg.Results.downsample_factor ~= 1
		add_to_call('-ds %d',arg.Results.downsample_factor);
	end

	if ~isempty(arg.Results.cal_file)
		add_to_call('-cal %s',arg.Results.cal_file);
	end

	if ~isempty(arg.Results.ctc_file)
		add_to_call('-ctc %s',arg.Results.ctc_file);
	end

	if ~isempty(arg.Results.trans_ref_file) % ICG ADDITION, ALLOW TRANFORM TO REFERENCE FILE
		add_to_call('-trans %s',arg.Results.trans_ref_file);
	end

	if ~isempty(arg.Results.extra)
		add_to_call(arg.Results.extra);
	end

	if exist(fif_out)
		fprintf('Deleting existing file: %s\n',fif_out);
		delete(fif_out);
	end

	runcmd('%s &> %s',maxfilter_call,log_out) % Capture log output
	maxfilter_output = fileread(log_out);

end
