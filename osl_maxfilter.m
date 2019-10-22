function [fif_out,bad_times,head_out,final_offset] = osl_maxfilter(fif_in,method,varargin)
    % Run MaxFilter
    % 
    % INPUTS
    % - fif_in : File name of an Elekta fif file
    % - method : One of 'nosss' (performs downsampling for bad channel detection), 'sss', or 'tsss'
    % - varargin : Additional command line options, see inputParser below for details
    %
    % This function will automatically create a head position file and log file in the
    % same directory as the output file.
    %
    % OUTPUTS
    % - fif output file
    % - bad times matrix extracted from log file to 
    % - final offset - factor to add to raw times e.g. from osl_headpos to match SPM
    %
    % Modified by Romesh Abeysuriya 2017
    
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
    arg.addParameter('downsample_factor',[])
    arg.addParameter('fif_out',[]) % Output file name, default is in the same directory with method prepended
    arg.addParameter('cal_file',[]) % :  full path to a fine calibration file to use
    arg.addParameter('ctc_file',[]) % ctc_file (optional):  full path to a cross-talk matrix file to use
    arg.addParameter('trans_ref_file',[]) % transformation to specified reference
    arg.addParameter('extra',[]) % Just stick stuff at the end
    arg.addParameter('verbose',false) % If true, print log to stdout while maxfilter is running

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

    assert(any(strcmp(method,{'nosss','sss','tsss'})),'Method must be one of - nosss, sss, tsss')
    
    if isempty(arg.Results.fif_out)
        base_path = fullfile(dirpath,sprintf('%s_%s',method,fname));
    else
        % Need to remove the fif extension, if present
        [outdir,outfname] = fileparts(arg.Results.fif_out);
        base_path = fullfile(outdir,outfname);
    end

    fif_out = [base_path '.fif'];
    head_out = [base_path '_headpos.txt'];

    add_to_call('-v'); % Enable verbose output
    add_to_call('-f %s',fif_in); % Specify input
    add_to_call('-o %s',fif_out); % Specify output
    add_to_call('-hp %s',head_out); % Output HPI 

    if ~isempty(arg.Results.badchannels)
        % Remove MEG from channel names if present
        badchans = regexp(arg.Results.badchannels,'[0-9]*','match','once'); % Extract the channel numbers
        badchans = sprintf('%s ',badchans{:}); % Join them
        add_to_call('-bad %s',badchans);
    end

    if arg.Results.movecomp
        add_to_call('-movecomp');
    end

    if ~arg.Results.autobad
        add_to_call('-autobad off');
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

    switch method
        case 'nosss'
            add_to_call('-nosss');
            fprintf('Calling maxfilter WITHOUT SSS\n');
            if isempty(arg.Results.downsample_factor)
                add_to_call('-ds %d',4); % By default, do downsampling with nosss
            end
        case 'sss'
            fprintf('Calling maxfilter with SSS\n');
        case 'tsss'
            add_to_call('-st %d -corr %g',arg.Results.st_win,arg.Results.st_corr);
            fprintf('Calling maxfilter with tSSS\n');
        otherwise
            error('Unrecognized method')
    end

    if arg.Results.downsample_factor ~= 1
        add_to_call('-ds %d',arg.Results.downsample_factor);
    end
    
    if osl_util.isfile(fif_out)
        fprintf('Deleting existing file: %s\n',fif_out);
        delete(fif_out);
    end

    % Redirect output to log files and optionally stdout as well (via tee)
    stdout_log = [base_path '_log.txt'];
    stderr_log = [base_path '_err.txt'];
    add_to_call('2> %s',stderr_log);
    if arg.Results.verbose
        add_to_call('| tee %s',stdout_log);
    else
        add_to_call('1> %s',stdout_log);
    end

    disp(maxfilter_call);
    ret = system(maxfilter_call); % Capture log output

    if ~osl_util.isfile(fif_out)
        fprintf(2,'Maxfilter did not generate the output file - something probably went wrong. Check %s\n',stderr_log);
    end
    
    bad_times = [];
    final_offset = 0;

    if any(strcmp(method,{'sss','tsss'}))
        try
            [bad_times,final_offset] = read_bad_times_from_maxfilter(stdout_log);
        catch ME
            fprintf(2,'Error parsing log file! This needs to be investigated, you will need to manually mark bad times\n');
            fprintf(2,'Error was:\n%s\n',ME.message)
        end
    end

end
