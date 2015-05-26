function fif_out = osl_call_maxfilter_remote(S, rmaxf_port, rmaxf_quit)
% Calls Elekta MaxFilter to apply signal space separation (SSS) to input
% .fif file.
%
% Syntax: FIF_OUT = osl_call_maxfilter(S, rmaxf_port, rmaxf_quit)
%
% S needs to contain:
%   -   fif: full name of input fif file (inc. path if not in present
%       working directory, with or without .fif file extention.
%
%   -   fif_out (optional): manually specify name of output file.
%
%   -   logfile (optional): set to 1 to produce log file of SSS output. 
%       Log file name is fif_out_log.txt. Set to 0 for no log file (default).
%
%   -   downsample_factor (optional): integer number >1. Suggested value of 
%       4 for 1000KHz data. If S.downsample_factor is not set then no 
%       downsampling occurs.
%
%   -   spmfile (optional): Include the name of corresponding SPM format 
%       data. Maxfilter
%       will use bad channel information contained in SPM object to ignore bad
%       channels.
%
%   -   nosss (optional): set to 1 to call Maxfilter without applying SSS.
%
%   -   maxfilt_dir (recommended): explicity tell OMT where to find
%       MaxFilter. Defaults to S.maxfilt_dir = '/neuro/bin/util'.
%
%   -   movement_compensation (optional): set to 1 to use cHPI movement
%       compensation or 0 (default) to not use it. Must be combined with SSS.
%
%   -   movecomp_call (optional, default = ''): string specifying a manual 
%       movement compensation call. Setting this field overrides 
%       movement_compensation. Must be combined with SSS.
%
%   -   autobad_off (optional): set to 1 to switch off MaxFilter's automated
%       badchannel detection. Default = 0.
%
%   -   trans_ref_file (optional, default = ''): allow transformation to the 
%       specified reference .fif file using the -trans maxfilter option
%
%   -   headpos (optional, default = 0): set to 1 to output head position 
%       information from the recording to a text file.
%
%   -   cal_file (optional):  full path to a fine calibration file to use
%
%   -   ctc_file (optional):  full path to a cross-talk matrix file to use
%
%
% OUTPUTS:
%   - fif_out: the name of the output file.
% HL 061011


%>SB
rmaxf_remote = 0;
% rmaxf_client = '/net/horus/data/OHBA/legacy/licence-workaround/remote_maxf_client.pl';
rmaxf_client = '/net/aton/data/OHBA/legacy/licence-workaround/remote_maxf_client.pl';



if 1 < nargin
    rmaxf_remote = 1;
    
    if ~isnumeric(rmaxf_port) || ~isscalar(rmaxf_port)
        error('port must be a scalar number');
    end
end

if 2 < nargin
    [rmstat, ~] = system([rmaxf_client ' "port=' num2str(rmaxf_port) ' quit"']);
    displayError(rmstat);
    fif_out = [];
    return;
end
%<SB

if ~isfield(S, 'fif')
    error('No .fif file specified.')
else
    [pathstr,filestr,~] = fileparts(S.fif);
    S.fif = fullfile(pathstr,filestr);
end

if ~isfield(S,'fif_out')
    fif_out=[S.fif '_sss'];
else
    fif_out=S.fif_out;
end

if isfield(S,'cal_file')
    cal_call = [' -cal ' S.cal_file ];
else
    cal_call = '';
end

if isfield(S,'ctc_file')
    ctc_call = [' -ctc ' S.ctc_file ];
else
    cal_call = '';
end

if isfield(S,'downsample_factor')
    if S.downsample_factor==1
        ds_call='';
    else
        ds_call=[' -ds ' num2str(S.downsample_factor)];
    end
else
    ds_call='';
end

if isfield(S,'spmfile')
    D=spm_eeg_load(S.spmfile);
    badchans=[]; chanind=D.badchannels;
    for i=1:length(D.badchannels)
        chanlab=D.chanlabels(chanind(i));
        badchans=[badchans ' ' chanlab{1}(4:end)];
    end
    if ~isempty(badchans)
    badchan_call = [' -bad ' badchans];
    else
        badchan_call = '';
    end
else
    badchan_call='';
end
if isfield(S,'nosss')
    if S.nosss==1
        nosss_call=[' -nosss'];
        if ~isfield(S,'fif_out')
            fif_out=[S.fif '_nosss'];
        end
    else
        nosss_call='';
    end
else
    nosss_call='';
end
if isfield(S,'maxfilt_dir')
    max_dir=[S.maxfilt_dir];
else
    max_dir='/neuro/bin/util';
    warning(['assuming Maxfilter found in ' max_dir]);
end
if isfield(S,'logfile')
    if S.logfile==1
        if rmaxf_remote == 0
            log_call=[' >& ' fif_out '_log.txt'];
        else
            log_call=[fif_out '_log.txt'];
        end
    else
      log_call='';
    end
else
    log_call='';
end
if isfield(S,'movecomp_call') % ICG ADDITION, ALLOW CUSTOM MOVECOMP CALLS
    movecomp_call=S.movecomp_call;
    
    if rmaxf_remote == 1
        warning('remote custom movecomp class not fully supported, please check result');
    end
    
    if isfield(S,'movement_compensation') && S.movement_compensation==0
        warning('Movement compensation call is specifed, even though movemen_compensation == 0!');
    end
elseif isfield(S,'movement_compensation')
    if S.movement_compensation==1 && isempty(nosss_call)
        movecomp_call=[' -movecomp'];
    else
        movecomp_call='';
    end
else
    movecomp_call='';
end
if isfield(S,'trans_ref_file') % ICG ADDITION, ALLOW TRANFORM TO REFERENCE FILE
    if rmaxf_remote == 0
        trans_call= [' -trans ' S.trans_ref_file '.fif'];
    else
        trans_call= [S.trans_ref_file '.fif'];
    end
else
    trans_call='';
end
if isfield(S,'headpos') % ICG ADDITION, ALLOW OUTPUT OF HEAD MOVEMENT INFO
    if S.headpos==1
        if rmaxf_remote == 0
            headpos_call = [' -hp ' fif_out '_headpos.txt'];
        else
            headpos_call = [fif_out '_headpos.txt'];
        end
        
        if isempty(movecomp_call)
            warning('''-hp'': cHPI option(s) not found');
        end
    else
        headpos_call='';
    end
else 
    headpos_call='';
end
if isfield(S,'autobad_off')
    if S.autobad_off==1 && isempty(nosss_call)
        autobad_call=[' -autobad off'];
    else
        autobad_call='';
    end
else
    autobad_call='';
end

if rmaxf_remote == 0
    maxfilter_call=[max_dir '/maxfilter -f ' S.fif '.fif -o ' fif_out '.fif'  nosss_call ds_call badchan_call movecomp_call autobad_call trans_call headpos_call ' -format float -v ' cal_call ctc_call log_call]
    runcmd(maxfilter_call)
    return
end

rmaxf_cmd = [rmaxf_client ' "port=' num2str(rmaxf_port) ' max=' ...
             max_dir '/maxfilter'                               ...
             nosss_call                                         ...
             ds_call                                            ...
             badchan_call                                       ...
             movecomp_call                                      ...
             autobad_call                                       ...
             ' -format float -v'                                ...
             ' #F:' num2str(int32([S.fif '.fif']))              ...
             ' #O:' num2str(int32([fif_out '.fif']))];

if 0 < length(trans_call)
    rmaxf_cmd = [rmaxf_cmd ' #T:' num2str(int32(trans_call))];
end

if 0 < length(headpos_call)
    rmaxf_cmd = [rmaxf_cmd ' #H:' num2str(int32(headpos_call))];
end

if 0 < length(log_call)
    rmaxf_cmd = [rmaxf_cmd ' #L:' num2str(int32(log_call))];
end

rmaxf_cmd   = [rmaxf_cmd '"'];
[rmstat, rmout] = system(rmaxf_cmd);
disp(rmout);
displayError(rmstat);
end

function displayError(status)
switch status
    case 1
        error('local: argument not valid');
    case 2
        error('local: user not identified');
    case 3
        error('local: connection failed');
    case 4
        error('local: temporary folder is not accessible');
    case 5
        error('local: input file(s) not specified/readable');
    case 6
        error('local: sss-output file not specified');
    case 7
        error('local temporary file already exists');

    case 10
        error('remote: request not valid');
    case 11
        error('remote: user not identified');
    case 12
        warning('remote: maxfilter exited with error, please check files');
        
    case 20
        warning('local: server did not acknowledge job, please check files');
        
    case 21
        warning('local: output file(s) not found/writable, please check files');
end
end
